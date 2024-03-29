# JASPER software for polishing genome assemblies and creating personalized reference genomes

JASPER (Jellyfish based Assembly Sequence Polisher for Error Reduction) is an efficient polishing tool for draft genomes.  It uses accurate reads (PacBio HiFi or Illumina) to evaluate consensus quality and correct consensus errors in genome assemblies.  JASPER is substantially faster than polishing methods based on sequence alignment, and more accurate than currently available k-mer based methods.  The efficiency and scalability of JASPER allows one to use it to create personalized reference genomes for specific populations very efficiently, even for large sequenced populations, by polishing the reference genome, such as GRCh38 or chm13v2.0 for human, with Illumina reads sequenced from many individuals from the population. 

Please see this manuscript for more details: Guo A, Salzberg SL, Zimin AV. JASPER: A fast genome polishing tool that improves accuracy of genome assemblies. PLoS Comput Biol. 2023 Mar 31;19(3):e1011032. doi: 10.1371/journal.pcbi.1011032. PMID: 37000853; PMCID: PMC10096238.

### Note: This version of JASPER does not install Jellyfish, a required dependency. For convenience we created an integrated version that includes all required dependencies.  You can install it from here: https://github.com/alekseyzimin/JASPER_release/releases
To install, download the release tarball JASPER-v1.0.4.tar.gz and then run:
```
tar xzf JASPER-v1.0.4.tar.gz
cd JASPER-v1.0.4
./install.sh
```
This will compule jellyfish and configure and install JASPER.  You can then run JASPER from /path_to/JASPER-v1.0.4/bin/jasper.sh. 

## Dependencies
* Python 3
* Jellyfish version 2 or above (https://github.com/gmarcais/Jellyfish)
* Biopython (https://biopython.org/)

## Installation
To install, please download the latest release tarball from the Releases section, and use "tar -xzf" to unpack the archive.  Then cd to the resulting folder and run ./configure --prefix=$PWD && make install.  For example:
```shell
tar -xzf jasper-1.0.3.tar.gz
cd jasper-1.0.3
./configure --prefix=$PWD
make install
```
JASPER will be available as $PWD/bin/jasper.sh.  Upon successful install, you can type $PWD/bin/jasper.sh -h to get usage information.
    
JASPER uses Python binding of Jellyfish. To configure Jellyfish Python binding, download and install the latest release of Jellyfish from https://github.com/gmarcais/Jellyfish, and then configure and compile with:

```shell
./configure --enable-python-binding 
make -j 4
sudo make install
```
By default, Jellyfish installs in /usr/local and the bindings are installed in the proper system location. The "python" command must alias to "python3" or "python3.x" in the system for Jellyfish binding to install properly. If the "python" is aliased to "python3.8" then Jellyfish binding will install into python3.8/site-packages/. Please type ``` which python ``` to ensure that "python" points to a python3.x version before running the command below. Furthermore, if you do not have root access, you can pass the --prefix switch to have the bindings install in the given directory. For example:

```shell
./configure --prefix=$HOME --enable-python-binding
make -j 4
make install
```
This will install the python binding in $HOME/lib/python3.8/site-packages (adjust based on your Python version). 

Then, an environment variable PYTHONPATH needs to be set. For example:

```shell
export PYTHONPATH=$HOME/lib/python3.8/site-packages 
```
    
## Running JASPER
To run JASPER, execute <PATH>/bin/jasper.sh with the following options:\
(default value in (), *required):\
```
-b, --batch=uint64               Desired batch size for the query (default value based on number of threads and assembly size). For the efficiency of jellyfish database loading, the max number of batches is limited to 8.
-t, --threads=uint32             Number of threads (1)
-a --assembly                    *Path to the assembly file
-j --jf                          Path to the jellyfish database file. Required if --reads is not provided
-r --reads                       Path to the file(s) containing the reads to construct a jellyfish database. If two or more files are provided, please enclose the list with single-quotes, e.g. -r '/path_to/file1.fastq /path_to/file2.fastq'. Both fasta and fastq formats are acceptable. Required if --jf is not provided.
-k, --kmer=uint64                k-mer size (37)
-p, --num_passes=utint16         The number of iterations of running jasper for fixing (2). A number smaller than 4 is usually more than sufficient
-h, --help                       This message
-v, --verbose                    Output information (False)
-d. --debug                      Debug mode. If supplied, all the _iter*batch*csv and _iter*batch*fa.temp files will be kept for debugging 
```
When JASPER runs with --reads, it creates mer_counts.jf Jellyfish k-mer count database from the reads. On subsequent runs in the same folder, if mer_counts.jf exists, it re-uses the mer_counts.jf.  Here is an example of how to run JASPER:

```shell
/PATH/bin/jasper.sh --reads '/path/read1.fastq /path/read2.fastq' -a assembly.fasta -k 25 -t 16 -p 4 1>jasper.out 2>&1
```
This command will polish assembly.fasta using reads /path/read1.fastq and /path/read2.fastq with 16 threads and with 4 passes of polishing. jasper.out will contain the diagnostic output including the QV value before and after polishing, and the polished assembly will be output as assembly.fasta.fixed.fasta.  
