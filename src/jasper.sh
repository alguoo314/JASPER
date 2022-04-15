#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
NUM_THREADS=1
CMD=jasper.py
BATCH_SIZE=0
PASSES=2
KMER=25
JF_SIZE=0
QUERY="random.fa"
QUERY_FN="random.fa"
READS="random.fastq"
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
kill -9 0
exit 1
}

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

function usage {
    echo "Usage: jasper.sh [options]"
    echo "Options:"
    echo "Options (default value in (), *required):"
    echo "-b, --batch=uint64               Desired batch size for the query (default value based on number of threads and assembly size). For the efficiency of jellyfish database loading, the max number of batches is limited to 8."
    echo  "-t, --threads=uint32             Number of threads (1)"
    echo "-a --assembly                    *Path to the assembly file"
    echo "-j --jf                          Path to the jellyfish database file. Required if --reads is not provided"
    echo "-r --reads                       Path to the file(s) containing the reads to construct a jellyfish database. If two or more files are provided, please enclose the list with single-quotes, e.g. -r '/path_to/file1.fa /path_to/file2.fa'. Required if --jf is not provided"
    echo "-k, --kmer=uint64                k-mer size (25)"
    echo "-p, --num_passes=utint16         The number of iterations of running jasper for fixing (2). A number smaller than 6 is usually more than sufficient" 
    echo "-h, --help                       This message"
    echo "-v, --verbose                    Output information (False)"
    
}

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -b|--batch)
            export BATCH_SIZE="$2"
            shift
            ;;
        -t|--threads)
            export NUM_THREADS="$2"
            shift
            ;;
        -a|--assembly)
            export QUERY="$2"
	    export QUERY_FN=`basename $QUERY`
            shift
            ;;
        -j|--jf)
            JF_DB="$2";
            shift
            ;;
	-r|--reads)
	    export READS="$2"
	    export JF_SIZE=`stat -c%s $READS |awk '{n+=$1}END{print int(n/100)}'`
	    shift
	    ;;
        -p|--num_passes)
            export PASSES="$2";
            shift
            ;;
	-k|--kmer)
            export KMER="$2";
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done



#calculate the threshold
if [ ! -s $PYTHONPATH/jellyfish.py ];then
error_exit "jellyfish.py not found on the python library path $PYTHONPATH, or variable \$PYTHONPATH is not set; please refer to https://github.com/gmarcais/Jellyfish for instructions on setting the variable."
fi

if [ ! -s $MYPATH/jasper.py ];then
error_exit "jasper.py not found in $MYPATH. jasper.py must be in the directory as jasper.sh.  Please make reinstall JASPER."
fi

if [ ! -s $QUERY ];then
error_exit "The query file does not exist. Please supply a valid fasta file to be polished with -a option."
fi

#Create database if the reads file is given instead of a database file                                       
if [ -z ${JF_DB+x} ];then
    if [ -n ${READS+x} ];then
	for filename in $READS
	do
	    if  [ ! -s $filename ];then
		error_exit "The reads file  $filename does not exist. Please supply a series of valid reads files separated by space and wrapped in one pair of quotation marks."
	     fi
	done
        JF_DB="mer_counts.jf"
        if [ -s "mer_counts.jf" ];then
          log "Using existing jellyfish database mer_counts.jf"
        else
          log "Creating jellyfish database mer_counts.jf"
          zcat -f $READS | jellyfish count -C -s $JF_SIZE -m 25 -o $JF_DB -t $NUM_THREADS /dev/stdin
        fi
    else
        error_exit "Either a jf database or files of polishing reads must be provided in the argument."
    fi
fi


if [ $NUM_THREADS -lt 9 ];then
    NEW_NUM_THREADS=$NUM_THREADS
else
    NEW_NUM_THREADS=8
fi
if ! [[ $BATCH_SIZE =~ ^[0-9]+$ ]];then
    log "BATCH SIZE supplied is not a positive integer. Calculating BATCH SIZE from QUERY SIZE"
    BATCH_SIZE=0
fi

BS=`grep -v '^>' $QUERY | tr -d '\n' |wc |awk '{print int($3/'$NEW_NUM_THREADS'*.9)}'`
if [ $BS -gt $BATCH_SIZE ];then
    BATCH_SIZE=$BS
fi    
log "Using BATCH SIZE $BATCH_SIZE"

if ! [[ $(($PASSES-1)) =~ ^[0-9]+$ ]];then
error_exit "The number of passes supplied by -p must be a positive integer"
fi

if ! [[ $(($KMER-1)) =~ ^[0-9]+$ ]];then
error_exit "The k-mer size supplied by -k must be a positive integer"
fi

if [ ! -e jasper.threshold.success ];then
log "Determining lower threshold for bad kmers"
jellyfish histo -t $NUM_THREADS $JF_DB > jfhisto.csv && \
jellyfish.py  jfhisto.csv > threshold.txt && \
rm jfhisto.csv && \
rm -f jasper.correct.success && \
touch jasper.threshold.success || error_exit "Computing threshold failed"
fi

LAST_IT=$(($PASSES-1))

#create batches
if [ ! -e jasper.split.success ];then 
log "Splitting query into batches for parallel execution"
rm -f $QUERY_FN.batch.*.fa && \
perl -ane 'BEGIN{$seq="";$bs=int('$BATCH_SIZE');}{if($F[0] =~ /^>/){if(not($seq eq "")){for($ci=0;$ci<length($seq);$ci+=$bs){print "$ctg:$ci\n",substr($seq,$ci,$bs),"\n";}}$ctg=$F[0];$seq=""}else{$seq.=$F[0]}}END{if(not($seq eq "")){for($ci=0;$ci<length($seq);$ci+=$bs){print "$ctg:$ci\n",substr($seq,$ci,$bs),"\n";}}}' $QUERY | \
perl -ane 'BEGIN{$batch_index=0;$output=0;open(FILE,">'$QUERY_FN'.batch.".$batch_index.".fa");}{if($F[0]=~/^>/){if($output>int('$BATCH_SIZE')){close(FILE);$batch_index++;open(FILE,">'$QUERY_FN'.batch.".$batch_index.".fa");$output=0;}}else{$output+=length($F[0]);}print FILE join(" ",@F),"\n";}' && \
rm -f jasper.correct.success && \
touch jasper.split.success || error_exit "Splitting files failed, do you have enough dick space?"
fi

if [ ! -e jasper.correct.success ];then
log "Polishing"
cat $JF_DB > /dev/null && \
echo "#!/bin/bash" >run_jasper.sh && \
echo "$CMD --db $JF_DB --query \$1 --ksize $KMER -p $PASSES --fix --fout \$1.fix.csv -ff \$1.fixed.fa.tmp --test -thre `head -n 1 threshold.txt| awk '{print $1}'` 1>jasper.out 2>jasper.err && mv _iter${LAST_IT}_\$1.fixed.fa.tmp _iter${LAST_IT}_\$1.fixed.fa" >> run_jasper.sh && \
chmod 0755 run_jasper.sh && \
ls $QUERY_FN.batch.*.fa | xargs -P $NEW_NUM_THREADS -I{} ./run_jasper.sh {} && \
rm -f run_jasper.sh && \
rm -f jasper.join.success && \
touch jasper.correct.success || error_exit "Polishing failed"
fi 

if [ ! -e jasper.join.success ];then
log "Joining"
cat _iter${LAST_IT}_$QUERY_FN.batch.*.fa.fixed.fa | perl -ane 'BEGIN{$seq="";$bs=int('$BATCH_SIZE');$bs=1 if($bs<=0);}{if($F[0] =~ /^>/){if(not($seq eq "")){$h{$ctg}=$seq;$seq=""}$ctg=$F[0]}else{$seq.=$F[0]}}END{$h{$ctg}=$seq;foreach $c(keys %h){if($c =~ /\:0$/){@f=split(/:/,$c);$ctg=join(":",@f[0..($#f-1)]);print "$ctg\n";$b=0;while(defined($h{$ctg.":$b"})){print $h{$ctg.":$b"};$b+=$bs;}print "\n";}}}' > $QUERY_FN.fixed.fasta.tmp && mv $QUERY_FN.fixed.fasta.tmp $QUERY_FN.fixed.fasta && \
rm -f _iter*_$QUERY_FN.batch.*.fa.fixed.fa _iter*_$QUERY_FN.batch.*.fa.fixed.fa.tmp && \
#cat _iter*_$QUERY_FN.batch.*.fa.fix.csv > $QUERY_FN.fixes.csv.tmp && mv $QUERY_FN.fixes.csv.tmp $QUERY_FN.fixes.csv && \
awk 'NR==1 || FNR>1' _iter*_$QUERY_FN.batch.*.fa.fix.csv  > $QUERY_FN.fixes.csv.tmp && mv $QUERY_FN.fixes.csv.tmp $QUERY_FN.fixes.csv && \
rm -f _iter*_$QUERY_FN.batch.*.fa.fix.csv && \
rm -f $QUERY_FN.batch.*.fa && \
touch jasper.join.success || error_exit "Joining failed"
fi

#outputting Q value
total_error_kmers_before=$(awk '{total+=$1}END{print total}' 0qValCalcHelper.csv)
total_kmers_before=$(awk '{total+=$2}END{print total}' 0qValCalcHelper.csv)
total_error_kmers_after=$(awk '{total+=$1}END{print total}' ${PASSES}qValCalcHelper.csv)
total_kmers_after=$(awk '{total+=$2}END{print total}' ${PASSES}qValCalcHelper.csv)
pgood_before=$(echo "scale=10; 1-$total_error_kmers_before/$total_kmers_before" | bc)
error_rate_before=$(echo "scale=50; 1 - e(l($pgood_before)*(1/$KMER))" | bc -l)
if (( $(echo "$error_rate_before > 0" | bc -l) )); then
    Q_before=$(echo "scale=5; -10*l($error_rate_before) / l(10)" | bc -l)
else
    Q_before="Inf"

fi
log "Before Polishing: Q value = $Q_before"

pgood_after=$(echo "scale=10; 1-$total_error_kmers_after/$total_kmers_after" | bc)
error_rate_after=$(echo "scale=50; 1 - e(l($pgood_after)*(1/$KMER))" | bc -l)
if (( $(echo "$error_rate_after > 0" | bc -l) )); then
    Q_after=$(echo "scale=5; -10*l($error_rate_after) / l(10)" | bc -l)
else
    Q_after="Inf"
fi

log "After Polishing: Q value = $Q_after"
rm -f *qValCalcHelper.csv

log "Polished sequence is in $QUERY_FN.fixed.fasta"
