# Jasper

JASPER (Jellyfish based Assembly Sequence Polisher for Error Reduction) is an efficient polishing tool for draft genomes.  It uses accurate reads (PacBio HiFi or Illumina) to evaluate consensus quality and correct consensus errors in genome assemblies.  JASPER is substantially faster than polishing methods based on sequence alignment, and more accurate than currently available k-mer based methods.  The efficiency and scalability of JASPER allows one to use it to create personalized reference genomes for specific populations very efficiently, even for large sequenced populations, by polishing the reference genome, such as GRCh38 or chm13v2.0 for human, with Illumina reads sequenced from many individuals from the population. 

## Dependencies
* Python 3
* Jellyfish (https://github.com/gmarcais/Jellyfish)

## Installation
To install, please download the latest release tarball from the Releases section, and use "tar -xzf" to unpack the archive.  Then cd to the resulting folder and run ./configure --prefix=$PWD && make install.  For example:
```shell
tar -xzf jasper-1.0.0.tar.gz
cd jasper-1.0.0
./configure --prefix=$PWD
make install
```
JASPER will be available as $PWD/bin/jasper.sh.  Upon successful install, you can type $PWD/bin/jasper.sh -h to get usage information.
    
JASPER uses python binding of Jellyfish. To configure Jellyfish Python binding, download Jellyfish, and then configure and compile with:

```shell
./configure --enable-python-binding 
make -j 4
sudo make install
```
By default, Jellyfish is installed in /usr/local and the bindings are installed in the proper system location. If you do not have root acess, you can pass the --prefix switch is passed, to have the bindings install in the given directory. For example:

```shell
./configure --prefix=$HOME --enable-python-binding
make -j 4
make install
```
This will install the python binding in $HOME/lib/python2.7/site-packages (adjust based on your Python version).
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
-r --reads                       Path to the file(s) containing the reads to construct a jellyfish database. If two or more files are provided, please enclose the list with single-quotes, e.g. -r '/path_to/file1.fa /path_to/file2.fa'. Required if --jf is not provided
-k, --kmer=uint64                k-mer size (25)
-p, --num_passes=utint16         The number of iterations of running jasper for fixing (2). A number smaller than 6 is usually more than sufficient
-h, --help                       This message
-v, --verbose                    Output information (False)
-d. --debug                      Debug mode. If supplied, all the _iter*batch*csv and _iter*batch*fa.temp files will be kept for debugging 
```
When JASPER runs with --reads, it creates mer_counts.jf Jellyfish k-mer count database from the reads. On subsequent runs in the same folder, if mer_counts.jf exists, it re-uses the mer_counts.jf.  Here is an example of how to run JASPER:

```shell
/PATH/bin/jasper.py --reads '/path/read1.fastq /path/read2.fastq' -a assembly.fasta -k 25 -t 16 -p 4 1>jasper.out 2>&1
```
This command will polish assembly.fasta using reads /path/read1.fastq and /path/read2.fastq with 16 threads and with 4 passes of polishing. jasper.out will contain the diagnostic output including the QV value before and after polishing, and the polished assembly will be output as assembly.fasta.fixed.fasta.  
