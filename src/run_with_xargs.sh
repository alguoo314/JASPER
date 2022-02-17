#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
NUM_THREADS=1
CMD=jasper.py
BATCH_SIZE=0

if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
kill 0
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
    echo "Usage: bash run_with_xargs.sh [options]"
    echo "Options:"
    echo "Options (default value in (), *required):"
    echo "-b, --batch=uint64               Desired batch size for the query (default value based on number of threads and assembly size)"
    echo  "-t, --threads=uint32             Number of threads (1)"
    echo "-a --assembly                    *Path to the assembly file"
    echo "-j --jf                          Path to the jellyfish database file. Required if --reads is not provided"
    echo "-r --reads                       Path to the file(s) containing the reads to construct a jellyfish database. If two or more files are provided, please enclose the list with  a pair of quotation marks. Required if --jf is not provided"
    echo "-k, --kmer=uint64                *k-mer size"
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
            shift
            ;;
        -j|--jf)
            JF_DB="$2";
            shift
            ;;
	-r|--reads)
	    export READS="$2"
	    export JF_SIZE=`stat -c%s $READS |awk '{n+=$1}END{print int(n/4)}'`
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
error_exit "jellyfish.py not found at python library path $PYTHONPATH, or path is not set; please refer to https://github.com/gmarcais/Jellyfish for instructions on setting the variable."
fi

if [ $BATCH_SIZE -lt 1 ];then
  BATCH_SIZE=`grep -v '^>' $QUERY | tr -d '\n' |wc |awk '{print int($3/'$NUM_THREADS'*.7)}'`
  log "Using BATCH SIZE $BATCH_SIZE"
fi



#Create database if the reads file is given instead of a database file             
if [ -z ${JF_DB+x} ];then
    if [ -n ${READS+x} ];then
	stringarray=($READS)
	firstfile=${stringarray[0]}
	filename_with_ext=${firstfile##*/}
	filename=${filename_with_ext%.*}
	JF_DB="$filename.jf"
	zcat -f $READS | jellyfish count -C -s $JF_SIZE -m 25 -o $JF_DB -t $NUM_THREADS
    else
	error_exit "Either a jf database or files of reads must be provided in the argument."
    fi
fi





if [ ! -e jasper.threshold.success ];then
log "Determining the lower threshold for bad kmers"
jellyfish histo -t $NUM_THREADS $JF_DB > jfhisto.csv && \
jellyfish.py  jfhisto.csv > threshold.txt && \
rm jfhisto.csv && touch jasper.threshold.success
fi



#create batches
if [ ! -e jasper.split.success ];then 
log "Splitting query into batches for parallel execution"
rm -f $QUERY.batch.*.fa && \
perl -ane 'BEGIN{$seq="";$bs=int('$BATCH_SIZE');}{if($F[0] =~ /^>/){if(not($seq eq "")){for($ci=0;$ci<length($seq);$ci+=$bs){print "$ctg:$ci\n",substr($seq,$ci,$bs),"\n";}}$ctg=$F[0];$seq=""}else{$seq.=$F[0]}}END{if(not($seq eq "")){for($ci=0;$ci<length($seq);$ci+=$bs){print "$ctg:$ci\n",substr($seq,$ci,$bs),"\n";}}}' $QUERY | \
perl -ane 'BEGIN{$batch_index=0;$output=0;open(FILE,">'$QUERY'.batch.".$batch_index.".fa");}{if($F[0]=~/^>/){if($output>int('$BATCH_SIZE')){close(FILE);$batch_index++;open(FILE,">'$QUERY'.batch.".$batch_index.".fa");$output=0;}}else{$output+=length($F[0]);}print FILE join(" ",@F),"\n";}' && \
touch jasper.split.success
fi

if [ ! -e jasper.correct.success ];then
log "Polishing"
cat $JF_DB > /dev/null && \
echo "#!/bin/bash" >run.sh && \
echo "$CMD --db $JF_DB --query \$1 --ksize 25 -p $PASSES --fix --fout \$1.fix.csv -ff \$1.fixed.fa.tmp -thre `head -n 1 threshold.txt| awk '{print $1}'` 1>jasper.out 2>jasper.err && mv _iter1_\$1.fixed.fa.tmp _iter1_\$1.fixed.fa" >>run.sh && \
chmod 0755 run.sh && \
ls $QUERY.batch.*.fa | xargs -P $NUM_THREADS -I{} ./run.sh {} && \
rm -f run.sh && \
touch jasper.correct.success
fi 

if [ ! -e jasper.join.success ];then
log "Joining"
cat _iter1_$QUERY.batch.*.fa.fixed.fa | perl -ane 'BEGIN{$seq="";$bs=int('$BATCH_SIZE');$bs=1 if($bs<=0);}{if($F[0] =~ /^>/){if(not($seq eq "")){$h{$ctg}=$seq;$seq=""}$ctg=$F[0]}else{$seq.=$F[0]}}END{$h{$ctg}=$seq;foreach $c(keys %h){if($c =~ /\:0$/){@f=split(/:/,$c);$ctg=join(":",@f[0..($#f-1)]);print "$ctg\n";$b=0;while(defined($h{$ctg.":$b"})){print $h{$ctg.":$b"};$b+=$bs;}print "\n";}}}' > $QUERY.fixed.fasta.tmp && mv  $QUERY.fixed.fasta.tmp  $QUERY.fixed.fasta && \
rm -f _iter?_$QUERY.batch.*.fa.fixed.fa $QUERY.batch.*.fa && \
log "Polished sequence is in $QUERY.fixed.fasta" && \
touch jasper.join.success
fi

#cat $QUERY.fix.*.csv > $QUERY.fix.csv.tmp && mv $QUERY.fix.csv.tmp $QUERY.fix.csv && \
#cat $QUERY.test.*.csv > $QUERY.test.csv.tmp && mv $QUERY.test.csv.tmp $QUERY.test.csv && \
#rm -f  $QUERY.fix.*.csv $QUERY.test.*.csv $QUERY.fixed.*.fasta && \
#grep ^Q  jasper.*.out > jasper.out 

