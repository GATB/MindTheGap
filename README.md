# MindTheGap 

| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-macos-10.9.5-gcc-4.2.1/)

Travis CI : [![Build Status](https://travis-ci.org/GATB/MindTheGap.svg?branch=master)](https://travis-ci.org/GATB/MindTheGap)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is MindTheGap ?

MindTheGap  performs detection and assembly of **DNA insertion variants** in NGS read datasets with respect to a reference genome. It is designed to call insertions of any size, whether they are novel or duplicated, homozygous or heterozygous in the donor genome. It takes as input a set of reads and a reference genome. It outputs two sets of FASTA sequences: one is the set of breakpoints of detected insertion sites, the other is the set of assembled insertions for each breakpoint.

**New !** MindTheGap can also be used as a **genome assembly finishing tool**: it can fill the gaps between a set of input contigs without any a priori on their relative order and orientation. It outputs the results in a gfa file. 

MindTheGap is a [Genscale](http://team.inria.fr/genscale/) tool, built upon the [GATB](http://gatb.inria.fr/) C++ library, and developed by:
* Claire Lemaitre
* Cervin Guyomar
* Wesley Delage
* Guillaume Rizk
* Former developers: Rayan Chikhi, Pierre Marijon. 

# Installation instructions

## Requirements

CMake 3.1+; see http://www.cmake.org/cmake/resources/software.html

C++/11 capable compiler (e.g. gcc 4.7+, clang 3.5+, Apple/clang 6.0+)

## Getting the latest source code with git

    # get a local copy of MindTheGap source code
    git clone --recursive https://github.com/GATB/MindTheGap.git
    
    # compile the code
    cd MindTheGap
    sh INSTALL
    # the binary file is located in directory build/bin/
    ./build/bin/MindTheGap -help

Note: when updating your local repository with `git pull`, if you see that thirdparty/gatb-core has changed, you have to run also : `git submodule update`. 

## Installing a stable release

Retrieve a binary archive file from one of the official MindTheGap releases (see "Releases" tab on the Github web page); file name is `MindTheGap-vX.Y.Z-bin-Linux.tar.gz` (for Linux) or `MindTheGap-vX.Y.Z-bin-Darwin.tar.gz` (for MacOs).

    tar -zxf MindTheGap-vX.Y.Z-bin-Darwin.tar.gz
    cd MindTheGap-vX.Y.Z-bin-Darwin
    chmod u+x bin/MindTheGap
    
    # run a simple example
    ./bin/MindTheGap find -in data/reads_r1.fastq,data/reads_r2.fastq -ref data/reference.fasta -out example
    ./bin/MindTheGap fill -graph example.h5 -bkpt example.breakpoints -out example

In case the software does not run appropriately on your system, you should consider to install it from its source code. Retrieve the source archive file `MindTheGap-vX.Y.Z-Source.tar.gz`.

    tar -zxf MindTheGap-vX.Y.Z-Source.tar.gz
    cd MindTheGap-vX.Y.Z-Source
    sh INSTALL
    # the binary file is located in directory build/bin/
    ./build/bin/MindTheGap -help

## Using docker or conda

Pull the docker image of the latest release of MindTheGap:

    docker pull clemaitr/mindthegap

MindTheGap is also distributed as a [Bioconda package](https://anaconda.org/bioconda/mindthegap):

``` 
conda install -c bioconda mindthegap
```



# USER MANUAL	 

## Description

MindTheGap is a software that performs integrated detection and assembly of **genomic insertion variants** in NGS read datasets with respect to a reference genome. It is designed to call insertions of any size, whether they are novel or duplicated, homozygous or heterozygous in the donor genome. 

Alternatively and since release 2.1.0, MindTheGap can also be used as a **genome assembly finishing tool**. It is integrated as an essential step in the targeted assembly tool [MinYS](https://github.com/cguyomar/MinYS) (MineYourSymbiont in metagenomics datasets).

**Insertion variant detection**

It takes as input a set of reads and a reference genome. It outputs two sets of FASTA sequences: one is the set of breakpoints of detected insertion sites, the other is the set of assembled insertions for each breakpoint. For each breakpoint, MindTheGap either returns a single insertion sequence (when there is no assembly ambiguity), or a set of candidate insertion sequences (due to ambiguities) or nothing at all (when the insertion is too complex to be assembled).

For a detailed user manual specific to insertion variants see [doc/MindTheGap_insertion_caller.md](doc/MindTheGap_insertion_caller.md).

**Genome assembly gap-filling** (New feature !)

When given a set of reads and a set of contigs as input, MindTheGap tries to fill the gaps between all pairs of contigs by de novo local assembly without any a priori on their relative order and orientation. It outputs the results in gfa file. 

For a detailed user manual specific to contig gap-filling see [doc/MindTheGap_assembly.md](doc/MindTheGap_assembly.md).

**Performances**

MindTheGap performs de novo assembly using the [GATB](http://gatb.inria.fr) C++ library and inspired from algorithms from Minia. Hence, the computational resources required to run MindTheGap are significantly lower than that of other assemblers (for instance it uses less than 6GB of main memory for analyzing a full human NGS dataset).


For more details on the method and some recent results, see the [web page](http://gatb.inria.fr/software/mind-the-gap/).
	
## Usage and examples

MindTheGap is composed of two main modules : breakpoint detection (`find` module) and the local assembly of insertions or gaps (`fill` module). Both steps are implemented in a single executable, MindTheGap, and can be run independently by specifying the module name as follows :

    MindTheGap <module> [module options] 

1. **Basic command lines**

        #Find module:
        MindTheGap find (-in <reads.fq> | -graph <graph.h5>) -ref <reference.fa> [options]
        #To get help:
        MindTheGap find -help
	    
        #Fill module:
        MindTheGap fill (-in <reads.fq> | -graph <graph.h5>) (-bkpt <breakpoints.fa> | -contig <contigs.fa>) [options]
        #To get help:
        MindTheGap fill -help

2. **Examples**

   These examples can be run with the small datasets in directory `data/`
   
	**Example for insertion variant calling:**
   
	    #find
	    bin/MindTheGap find -in data/reads_r1.fastq,data/reads_r2.fastq -ref data/reference.fasta -out example
	    # 3 files are generated: 
	    #   example.h5 (de bruijn graph), 
	    #   example.othervariants.vcf (SNPs and deletion variants), 
	    #   example.breakpoints (breakpoints of insertion variants).
	    
	    #fill
	    bin/MindTheGap fill -graph example.h5 -bkpt example.breakpoints -out example
	    # 3 files are generated:
	    #   example.insertions.fasta (insertion sequences)
	    #   example.insertions.vcf (insertion variants)
	    #   example.info.txt (log file)
	
	**Example for gap-filling between contigs:**
	
	```
	build/bin/MindTheGap fill -in data/contig-reads.fasta.gz -contig data/contigs.fasta -abundance-min 3 -out contig_example
	# 4 files are generated
	#   contig_example.h5 (de bruijn graph)
	#   contig_example.insertions.fasta (gap-filling sequences)
	#   contig_example.gfa (genome graph)
	#   contig_example.info.txt (log file)
	```
	
	The usage of the `fill` module is a little bit different depending on the type of gap-filling : assembling insertion variants (using the `-bkpt`option with a breakpoint file) or gap-filling between contigs (using the `-contig` option with a contig fasta file). 

## Details

1. **Input sequencing read data**
	
	For both modules, read dataset(s) are first indexed in a De Bruijn graph. The input format of read dataset(s) is either the read files themselves (option `-in`), or the already computed de bruijn graph in hdf5 format (.h5) (option `-graph`).   
	NOTE: options `-in` and `-graph` are mutually exclusive, and one of these is mandatory.
	
	If the input is composed of several read files, they can be provided as a list of file paths separated by a comma or as a "file of file" (fof), that is a text file containing on each line the path to each read file. All read files will be treated as if concatenated in a single sample. The read file format can be fasta, fastq or gzipped. 
	
3. **de Bruijn graph creation options**

   In addition to input read set(s), the de Bruijn graph creation uses two main parameters, `-kmer-size` and `-abundance-min`: 

   * `-kmer-size`: the k-mer size [default '31']. By default, the largest kmer-size allowed is 128. To use k>128, you will need to re-compile MindTheGap with the two commands in the build directory: `cmake -DKSIZE_LIST="32 64 96 256" ..` and then `make`. To go back to default, replace 256 by 128. Note that increasing the range between two consecutive kmer-sizes in the list can have an impact on the size of the output h5 files (but none on the results).
* `-abundance-min`: the minimal abundance threshold, k-mers having less than this number of occurrences are discarded from the graph [default 'auto', ie. automatically inferred from the dataset]. 
   * `-abundance-max`: the maximal abundance threshold, k-mers having more than this number of occurrences are discarded from the graph [default '2147483647' ie. no limit].

6. **Computational resources options**

    Additional options are related to computational runtime and memory:
    
    * `-nb-cores`: number of cores to be used for computation [default '0', ie. all available cores will be used].
* `-max-memory`: max RAM memory for the graph creation (in MBytes)  [default '2000']. Increasing the memory will speed up the graph creation phase.
    * `-max-disk`: max usable disk space for the graph creation (in MBytes)  [default '0', ie. automatically set]. Kmers are counted by writting temporary files on the disk, to speed up the counting you can increase the usable disk space.
    
4. **MindTheGap Output**

    All the output files are prefixed either by a default name: "MindTheGap_Expe-[date:YY:MM:DD-HH:mm]" or by a user defined prefix (option `-out` of MindTheGap)
    Both MindTheGap modules generate the graph file if reads were given as input: 
    * a graph file (`.h5`). This is a binary file, to obtain information stored in it, you can use the utility program dbginfo located in your bin directory or in ext/gatb-core/bin/.

    `MindTheGap find` generates the following output files:
    * a breakpoint file (`.breakpoints`) in fasta format. 
    * a variant file (`.othervariants.vcf`) in vcf format. It contains SNPs and deletion events.

    `MindTheGap fill` generates the following output files:
    * a sequence file (`.insertions.fasta`) in fasta format. It contains the inserted sequences or contig gap-fills that were successfully assembled. 

    * an insertion variant file (`.insertions.vcf`) in vcf format, in the case of insertion variant detection. 

    * an assembly graph file (`.gfa`) in GFA format, in the case of contig gap-filling. It contains the original contigs and the obtained gap-fill sequences (nodes of the graph), together with their overlapping relationships (arcs of the graph).

    * a log file (`.info.txt`), a tabular file with some information about the filling process for each breakpoint/grap-fill. 

      

Other optional parameters and details on input and output file formats are given in [doc/MindTheGap_insertion_caller.md](doc/MindTheGap_insertion_caller.md) and [doc/MindTheGap_assembly.md](doc/MindTheGap_assembly.md), depending on the usage.



## Utility programs

Either in your `bin/` directory or in `ext/gatb-core/bin/`, you can find additional utility programs :
* `dbginfo` : to get information about a graph stored in a .h5 file
* `dbgh5` : to build a graph from read set(s) and obtain a .h5 file
* `h5dump` : to extract data stored in a .h5 file



## Reference

MindTheGap: integrated detection and assembly of short and long insertions. Guillaume Rizk, Anaïs Gouin, Rayan Chikhi and Claire Lemaitre. Bioinformatics 2014 30(24):3451-3457. http://bioinformatics.oxfordjournals.org/content/30/24/3451

[Web page](https://gatb.inria.fr/software/mind-the-gap/) with some updated results.


# Contact

To contact a developer, request help, or for any feedback on MindTheGap, please use the issue form of github: https://github.com/GATB/MindTheGap/issues

You can see all issues concerning MindTheGap [here](https://github.com/GATB/MindTheGap/issues) and GATB [here](https://www.biostars.org/t/GATB/).

If you do not have any github account, you can also send an email to claire dot lemaitre at inria dot fr
