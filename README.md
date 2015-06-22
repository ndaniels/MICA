ABOUT
=====
MICA (Metagenomic Inquiry Compressive Acceleration) is a family of programs for performing compressively-accelerated
metagenomic sequence searches based on BLASTX and DIAMOND.
MICA also includes compressively accelerated versionf of the BLASTP family of 
tools (including PSI-BLAST and DELTA-BLAST), as well as a compression tool (mica-compress)
for creating searchable, compressed databases based on an input FASTA file.

If you use MICA, please cite:

* Daniels N, Gallant A, Peng J, Cowen L, Baym M, Berger B
"Compressive Genomics for Protein Databases." Bioinformatics 29.13 (2013): i283-i290.

* Yu YW, Daniels N, Danko DC, Berger B 
"Entropy-scaling search of massive biological data." (2015) Submitted.

* Buchfink B, Chao X, Huson, D. "Fast and sensitive protein alignment using DIAMOND." Nature methods 12.1 (2015): 59-60.

MICA is licensed under the GNU public license version 2.0. If you would
like to license MICA in an environment where the GNU public license is
unacceptable (such as inclusion in a non-GPL software package) commercial
MICA licensing is available through MIT office of Technology Transfer.
Contact bab@mit.edu for more information.
Contact ndaniels@csail.mit.edu for issues involving the code.


QUICK EXAMPLE
=============
Assuming you have [Go](https://golang.org/), [DIAMOND](https://github.com/bbuchfink/diamond/) and [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) installed, here is a quick example of how to 
perform a compressively accelerated MICA search using a compressed database 
that has already been created.

    # Install MICA
    go get github.com/ndaniels/MICA/...

    # Download and extract the database. It is large and could take a while.
    # Make sure to check for a newer version!
    wget http://giant.csail.mit.edu/gems/nr-20140917-mica.tgz
    tar zxf nr-20140917-mica.tgz

    # Search.
    mica-xsearch --dmnd-fine=result.txt nr-20140917-mica query.fasta

There are more examples covering more use cases further down.


INSTALLATION
============
The easiest way to install is to download binaries compiled for your operating
system. No other dependencies are required (sans BLAST+ and DIAMOND, which 
should already be in your PATH).
They can be downloaded here: http://gems.csail.mit.edu/

Compiling from source is also easy; compiling MICA only requires that git 
and Go are installed. If Go is not already available via your package manager,
it can be installed from source by following the directions here:
http://golang.org/doc/install

Once Go is installed, you'll need to set your GOPATH, which is where MICA 
(and other Go packages) will be installed. We recommend running

    mkdir $HOME/go

And adding the following to your `~/.profile` or equivalent:

    export GOPATH="$HOME/go"
    export PATH="$PATH:$GOPATH/bin"

Finally, run the following command to download, compile and install CaBLASTP:

    go get github.com/ndaniels/MICA/...

The MICA executables should be installed in `$GOPATH/bin`.

MICA has been tested against Go 1.4.2.


EXECUTABLES
===========
There are seven binary executables in the MICA suite, also available as 
binaries for users without Go installed. They are: 

    mica-xsearch      A compressively accelerated translated search (like BLASTX),
                        which can use DIAMOND or BLASTX for fine search.
    mica-psearch      A compressively accelerated protein search (like BLASTP),
                        which can use DIAMOND or BLASTP for fine search.
    
    mica-compress     Compresses FASTA input files (such as nr.fasta or
                          nr.gz) into a compressed database for quick searching.

    mica-decompress   A rarely-needed inverse of mica-compress.

    mica-search       A compressively accelerated version of BLASTP.

    mica-psisearch    A compressively accelerated version of PSI-BLAST.

    mica-deltasearch  A compressively accelerated version of DELTA-BLAST.

Every executable can be run with the `--help` flag to get a list of command 
line options.


PREREQUISITES
=============
MICA boosts BLAST+ and DIAMOND protein-database search, and as such it is not 
completely self-contained. It relies on BLAST+ and DIAMOND.

To use MICA, you must already have BLAST+ 2.2 or later installed, so that
the BLAST binaries are in your PATH. DELTA-BLAST requires BLAST+ 2.2.26 or 
later and we recommend 2.2.27. DELTA-BLAST also requires an RPS database 
configured per NCBI's instructions.
You must also have DIAMOND installed (tested with DIAMOND 0.7.9) so that the
`diamond` binary is in your PATH.

We provide binaries for Mac OS X (64-bit intel, tested on OS X 10.10.3 and
built with Go 1.4.2) and Linux (64-bit intel/AMD, tested on Linux kernel 3.13.0 
and Go 1.4.2). With Go installed, MICA might work on Microsoft Windows but 
is untested and unsupported.

You do not need the Go compiler installed to use the binary distributions of
MICA.


ADDITIONAL FILES
================
As compression is compute-intensive, we provide an already-compressed database
based on NCBI's NR from September 17, 2014, which we will update quarterly.
Since the MICA compressed database format is actually a directory 
structure, we provide it as a .tar.gz file, so should be unarchived with 
`tar zxf nr-20140917-mica.tgz`.

The result will be a directory, 'nr-20140917-mica', which contains the 
various files necessary for MICA to run.

Should you wish to create your own compressed database, you would use the
mica-compress binary. The database we provide was created with:

    mica-compress --match-seq-id-threshold 60 --ext-seed-size 0
                      --ext-seq-id-threshold 50 --max-seeds 20 -p 40
                      nr-20140917-mica nr.fasta

Several of the command-line arguments are tuning parameters that affect the
run-time performance of compression.

The --max-seeds argument caps the size of the seeds table to, in this case, 20
gigabytes. Compressing large databases can require a great deal of RAM. A
significantly smaller cap will harm compression.

The --ext-seed-size argument allows for larger k-mer seeds without the memory
overhead associated with the larger size, by greedily requiring the additional
residues to be exact matches.

The --match-seq-id-threshold argument sets the sequence identity percentage
required for a match during compression.

The --ext-seq-id-threshold argument sets the sequence identity percentage
required for a single instance of extension during compression.

The -p argument simply sets the number of processor cores used during 
compression, and bears no relevance to the resulting compressed database.

In this case, the input file is `nr.fasta`, and the output name for the
compressed database is `nr-20140917-mica`.
Note that the compressed database is actually a directory that will be created
by `mica-compress`.


USAGE
=====
Run mica-compress -help, mica-xsearch -help, mica-psearch -help, 
mica-deltasearch -help, mica-search -help, 
or mica-psisearch -help for detailed help as to command-line arguments.


EXAMPLES
========
To perform a compressively accelerated DIAMOND search, you might do:

    mica-xsearch --dmnd-fine=result.txt
                         /path/to/mica_database /path/to/query.fasta

where:

    result.txt is the local file path to output the DIAMOND results

    /path/to/mica_database is the local file path to your MICA 
    compressed database (it will be the path to nr-20140917-mica if you are 
    using the provided September, 2014 database)

    /path/to/query.fasta is simply the local file path to the FASTA file you 
    wish to use as a query.

To perform a compressively accelerated BLASTX search, you might do:

    mica-xsearch /path/to/mica_database /path/to/query.fasta
      --blast-args -evalue 1e-7 -outfmt 6 > result.txt

where:

    /path/to/mica_database is the local file path to the MICA 
    compressed database,

    /path/to/query.fasta is the local file path to the FASTA file you wish to 
    use as a query,
    
    1e-7 is the BLAST e-value you wish to use,
    
    -outfmt 6 is the standard BLAST argument to produce tabular output, and
    
    result.txt is where you wish the results to appear.
    
    

Arguments the user wishes to pass to the underlying BLAST program, if BLAST is 
used for fine search, such as
adjusting the output format or the E-value threshold, may be passed via the
`--blast-args` flag.

For example, to specify XML output, one might run:

    mica-xsearch /path/to/mica_database /path/to/query.fasta
                    --blast-args -outfmt 5

Where `-outfmt 5` is, as indicated in the NCBI blastp user guide, the 
command-line argument for XML output.


REPORTING BUGS
==============
If you find any bugs or have any problems using CaBLASTP, please submit a bug
report on our issue tracker:

    https://github.com/ndaniels/MICA/issues

