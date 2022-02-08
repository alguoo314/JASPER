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
echo "Usage:"
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


#create batches
if [ ! -e jasper.split.success ];then 
log "Splitting query into batches for parallel execution"
if [ $BATCH_SIZE -lt 1 ];then
  BATCH_SIZE=`grep -v '^>' $QUERY |wc | tr -d '\n' |awk '{print int($3/'$NUM_THREADS'*.7)}'`
fi
log "Using BATCH SIZE $BATCH_SIZE"
rm -f $QUERY.batch.*.fa && \
perl -ane 'BEGIN{$seq="";$bs=int('$BATCH_SIZE');}{if($F[0] =~ /^>/){if(not($seq eq "")){for($ci=0;$ci<length($seq);$ci+=$bs){print "$ctg:$ci\n",substr($seq,$ci,$bs),"\n";}}$ctg=$F[0];$seq=""}else{$seq.=$F[0]}}END{if(not($seq eq "")){for($ci=0;$ci<length($seq);$ci+=$bs){print "$ctg:$ci\n",substr($seq,$ci,$bs),"\n";}}}' $QUERY | \
perl -ane 'BEGIN{$batch_index=0;$output=0;open(FILE,">'$QUERY'.batch.".$batch_index.".fa");}{if($F[0]=~/^>/){if($output>int('$BATCH_SIZE')){close(FILE);$batch_index++;open(FILE,">'$QUERY'.batch.".$batch_index.".fa");$output=0;}}else{$output+=length($F[0]);}print FILE join(" ",@F),"\n";}' && \
touch jasper.split.success
fi

if [ ! -e jasper.correct.success ];then
log "Polishing"
cat $JF_DB > /dev/null && \
echo "#!/bin/bash" >run.sh && \
echo "$CMD --db $JF_DB --query \$1 --ksize 25 --fix --fout \$1.fix.csv -ff \$1.fixed.fa.tmp -thre 10 -rep_thre 10000000000 1>jasper.out 2>jasper.err && mv _iter1_\$1.fixed.fa.tmp _iter1_\$1.fixed.fa" >>run.sh && \
chmod 0755 run.sh && \
ls $QUERY.batch.*.fa | xargs -P $NUM_THREADS -I{} ./run.sh {} && \
rm -f run.sh && \
touch jasper.correct.success
fi 

if [ ! -e jasper.join.success ];then
log "Joining"
cat _iter1_$QUERY.batch.*.fa.fixed.fa | perl -ane 'BEGIN{$seq="";$bs=int('$BATCH_SIZE');}{if($F[0] =~ /^>/){if(not($seq eq "")){$h{$ctg}=$seq;$seq=""}$ctg=$F[0]}else{$seq.=$F[0]}}END{$h{$ctg}=$seq;foreach $c(keys %h){if($c =~ /\:0$/){@f=split(/:/,$c);$ctg=join(":",@f[0..($#f-1)]);print "$ctg\n";$b=0;while(defined($h{$ctg.":$b"})){print $h{$ctg.":$b"};$b+=$bs;}print "\n";}}}' > $QUERY.fixed.fasta.tmp && mv  $QUERY.fixed.fasta.tmp  $QUERY.fixed.fasta && \
rm -f _iter?_$QUERY.batch.*.fa.fixed.fa $QUERY.batch.*.fa && \
log "Polished sequence is in $QUERY.fixed.fasta" && \
touch jasper.join.success
fi

#cat $QUERY.fix.*.csv > $QUERY.fix.csv.tmp && mv $QUERY.fix.csv.tmp $QUERY.fix.csv && \
#cat $QUERY.test.*.csv > $QUERY.test.csv.tmp && mv $QUERY.test.csv.tmp $QUERY.test.csv && \
#rm -f  $QUERY.fix.*.csv $QUERY.test.*.csv $QUERY.fixed.*.fasta && \
#grep ^Q  jasper.*.out > jasper.out 

