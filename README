# My_Project
## Evaluate genome assemblies with Jellyfish and k-mers

(short description) 

## Dependency
* Python 3
* Jellyfish (https://github.com/gmarcais/Jellyfish)



## Installation

### Stable Release

#### Direct installation 
1. Get a working Jellyfish and perform Python binding
See instructions in Jellyfish's manual (https://github.com/gmarcais/Jellyfish)

2. Download the release version of My_Project

Blah blah blah


## Run
Options:
  --db, the path to the .jf  database file. REQUIRED if --reads is not given.
  --fix, output the index of fixed bases and output the new sequence.
  -ff, --fixedfasta, the path to output the fixed assembly sequences. Default = fixed_seq.fasta. No use if --fix is not provided
  --fout, the path to output the index of the fixed bases (0-based). Default = fout.csv. No use if --fix is not provided
  --help,	display this help message
	-k, --ksize, the kmer size. Type = int. REQUIRED
  -q, --query, the path to the .fasta query file. REQUIRED
  --reads, the path to the .fasta file(s）containing the reads to build the jellyfish database. REQUIRED if --db is not given.
	-rep_thre, the threshold  of occurrance for a kmer at a repeitive region. Type = int. If not provided, will be be calculated from the .jf file.
	--test, output the indexes of bad kmers, total num of bad kmers, and an estimation for Q value.
  -thre, --threshold, the threshold for a bad kmer. Type = int. If not provided, will be be calculated from the .jf file.
  --tout, the path to output the locations of bad kmers (1-based). Default = tout.csv. No use if --test is not provided

*Note: 
1. One and only one between the contigs and database argument should be given.
2. If --contigs is not given, --db, --rep_thre, and --thre must all be given. The two thresholds can be manually calculated by the user from the Jellyfish histogram file:
jellyfish histo your.jf > your.histo
Look at the second column in the .histo file. Find the row of the first local minimum in the second column, then 0.5*the value in the first column in this row is the bad kmer threshold. 
Then locate the glocal maximum in the second column among the rows below the bad kmer threshold row, then the 2* the value in the first column in this row is the repeititve region threshold 

Below is showing examples how to run My_Program (Each # represent an integer)
If I have the reads to build the .jf file and would like the program to determine the two threshold. Also want both testing and fixing outputs
```shell
python My_Program.py --reads read1.fastq (read2.fastq, etc.)  --query assembly.fasta --ksize # --test --fix --fout fix.csv  --tout test.csv -ff fixed.fasta
```
If I have the reads to build the .jf file and want a specific repetitve region threshold (or bad kmer threshold, just providing the --threshold argument instead of -repthre)
```shell
python My_Program.py --reads read1.fastq (read2.fastq, etc.) -rep_thre # --query assembly.fasta --ksize # --test --fix --fout fix.csv  --tout test.csv -ff fixed.fasta
```
If I have the .jf file already along with the two thresholds. Also want both testing and fixing outputs
```shell
python My_Program.py --db mydb.jf  --threshold # -rep_thre # --query assembly.fasta --ksize # --test --fix --fout fix.csv  --tout test.csv -ff fixed.fasta 
```
