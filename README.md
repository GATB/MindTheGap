# MindTheGap 

| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/MindTheGap/job/tool-mindthegap-build-macos-10.9.5-gcc-4.2.1/)

Travis CI : [![Build Status](https://travis-ci.org/GATB/MindTheGap.svg?branch=master)](https://travis-ci.org/GATB/MindTheGap)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is MindTheGap ?

MindTheGap  performs detection and assembly of **DNA insertion variants** in NGS read datasets with respect to a reference genome. It is designed to call insertions of any size, whether they are novel or duplicated, homozygous or heterozygous in the donor genome. It takes as input a set of reads and a reference genome. It outputs two sets of FASTA sequences: one is the set of breakpoints of detected insertion sites, the other is the set of assembled insertions for each breakpoint.

**New !** MindTheGap can also be used as a **genome assembly finishing tool**: it can fill the gaps between a set of input contigs without any a priori on their relative order and orientation. It outputs the results in gfa file. 

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

## Using docker

Pull the docker image of the latest release of MindTheGap:

    docker pull clemaitr/mindthegap

# USER MANUAL	 

## Description

MindTheGap is a software that performs integrated detection and assembly of **genomic insertion variants** in NGS read datasets with respect to a reference genome. It is designed to call insertions of any size, whether they are novel or duplicated, homozygous or heterozygous in the donor genome. 

Alternatively and since release 2.1.0, MindTheGap can also be used as a **genome assembly finishing tool**.

### Insertion variant detection

It takes as input a set of reads and a reference genome. It outputs two sets of FASTA sequences: one is the set of breakpoints of detected insertion sites, the other is the set of assembled insertions for each breakpoint. For each breakpoint, MindTheGap either returns a single insertion sequence (when there is no assembly ambiguity), or a set of candidate insertion sequences (due to ambiguities) or nothing at all (when the insertion is too complex to be assembled).

Since version 2.0.0, MindTheGap can detect other types of variants, not only insertion events. These are homozygous SNPs and homozygous deletions of any size. They are detected by the find module of MindTheGap and are output separately in a VCF file. Importantly, even if the user is not interested in these types of variants, it is worth to detect them  since it can improve the recall of the insertion event detection algorithm : it is now possible to find insertion events that are located at less than k nucleotides from an other such variant.

### Genome assembly gap-filling

New feature !

When given a set of reads and a set of contigs as input, MindTheGap tries to fill the gaps between all pairs of contigs by de novo local assembly without any a priori on their relative order and orientation. It outputs the results in gfa file. 

### Performances

MindTheGap performs de novo assembly using the [GATB](http://gatb.inria.fr) C++ library and inspired from algorithms from Minia. Hence, the computational resources required to run MindTheGap are significantly lower than that of other assemblers (for instance it uses less than 6GB of main memory for analyzing a full human NGS dataset).


For more details on the method and some recent results, see the [web page](http://gatb.inria.fr/software/mind-the-gap/).
	
## Usage

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

2. **Input read data**

   For both modules, read dataset(s) are first indexed in a De Bruijn graph. The input format of read dataset(s) is either the read files themselves, or the already computed de bruijn graph in hdf5 format (.h5). In the first case, the option is `-in` and the user can provide the de Bruijn graph building options, in the second case the option is -graph and only options for the detection or assembly are to be given.  
   NOTE: options `-in` and `-graph` are mutually exclusive, and one of these is mandatory.
	
   If the input is composed of several read files, they can be provided as a list of file paths separated by a comma or as a "file of file" (fof), that is a text file containing on each line the path to each read file. All read files will treated as if concatenated in a single sample. The read file format can be fasta, fastq or gzipped. 
		
3. **de Bruijn graph creation options**

   In addition to input read set(s), the de Bruijn graph creation uses two main parameters, `-kmer-size` and `-abundance-min`: 

   * `-kmer-size`: the k-mer size [default '31']. By default, the largest kmer-size allowed is 128. To use k>128, you will need to re-compile MindTheGap with the two commands in the build directory: `cmake -DKSIZE_LIST="32 64 96 256" ..` and then `make`. To go back to default, replace 256 by 128. Note that increasing the range between two consecutive kmer-sizes in the list can have an impact on the size of the output h5 files (but none on the results).

   * `-abundance-min`: the minimal abundance threshold, k-mers having less than this number of occurrences are discarded from the graph [default 'auto', ie. automatically inferred from the dataset]. 

   * `-abundance-max`: the maximal abundance threshold, k-mers having more than this number of occurrences are discarded from the graph [default '2147483647' ie. no limit].
	
4. **Find module specific options**
    
    In addition to the read or graph files, the find module has one mandatory option `-ref` and several optional options:
    * `-ref`: the path to the reference genome file (in fasta format).
    * `-max-rep`: maximal repeat size allowed for fuzzy sites  [default '5']. 
    * `-snp-min-val`: minimal number of kmers to validate a SNP [default '5']. A SNP is validated if it is validated by at least this number of consecutive overlapping kmers. This corresponds also to the minimal distance between two consecutive SNPs to be able to detect both of them.
    * `-het-max-occ`: maximal number of occurrences of a (k-1)mer in the reference genome allowed for heterozyguous insertion breakpoints  [default '1']. In order to detect an heterozyguous insertion breakpoints, both flanking k-1-mers, at each side of the insertion site, must have strictly less than this number of occurrences in the reference genome. This prevents false positive predictions inside repeated regions. Warning : increasing this parameter may lead to numerous false positives (genomic approximate repeats).
    * `-no-[type]`: to disable the detection of certain types of variants.
    * `-[type]-only`: to detect only certain types of variants.
    
    NOTE: MindTheGap can find mainly homozygous variants, except for insertion variants for which it can also find heterozygous variants. Therefore -homo-only and -hete-only only apply to insertion variants.

5. **Fill module specific options**
    
    In addition to the read or graph files, the fill module has one other mandatory option, either `-bkpt` or `-contig` depending on the type of gap-filling : assembling insertion variants or gap-filling between contigs respectively: 	
    * `-bkpt`: the breakpoint file path. This is one of the output of the Find module and contains for each detected insertion site its left and right kmers from and to which the local assembly will be performed (see section E for details about the format).
	* `-contig`: the contig file path in fasta format. Note that only contigs larger than 2*kmerSize will be used.
	
	The fill module has several optional options:
    * `-max-nodes`: maximum number of nodes in contig graph  [default '100']. This arguments limits the computational time, this is especially useful for complex genomes.
    * `-max-length`: maximum length of insertions (nt)  [default '10000']. This arguments limits the computational time, this is especially useful for complex genomes.
	* `-overlap`: size of sequence overlap between input contigs in `-contig` mode [default '0']. To be specified if it is larger than the kmer size used for gap-filling (expert usage).

6. **MindTheGap Output**
    
    All the output files are prefixed either by a default name: "MindTheGap_Expe-[date:YY:MM:DD-HH:mm]" or by a user defined prefix (option `-out` of MindTheGap)
    Both MindTheGap modules generate the graph file if reads were given as input: 
    * a graph file (`.h5`). This is a binary file, to obtain information stored in it, you can use the utility program dbginfo located in your bin directory or in ext/gatb-core/bin/.
    
    `MindTheGap find` generates the following output files:
    * a breakpoint file (`.breakpoints`) in fasta format. It contains the breakpoint sequences of each detected insertion site. Each insertion site corresponds to 2 consecutive entries in the fasta file : sequences are the left and right side flanking kmers.
    * a variant file (`.othervariants.vcf`) in vcf format. It contains SNPs and deletion events.
    
    `MindTheGap fill` generates the following output files:
    * a sequence file (`.insertions.fasta`) in fasta format. It contains the inserted sequences or contig gap-fills that were successfully assembled. In the case of insertion variants, the location of each insertion on the reference genome can be found in its fasta header. In the case of contig gap-fills, the fasta header contains the source and target contigs with their relative orientation ("_Rc" for reversed). In both cases, the fasta header includes also information about each gap-fill such as its length, quality score and median kmer abundance.
    * an insertion variant file (`.insertions.vcf`) in vcf format, in the case of insertion variant detection. This file contains all information of assembled insertion variants as in the `.insertions.fasta` file but in a different format. Here, insertion site positions are 1-based and left-normalized according to the VCF format specifications (contrary to positions indicated in the `.breakpoints` and `insertions.fasta` files which are right-normalized). Normalization occurs when multiple positions are possible for a single variation due to a small repeat. 
	* an assembly graph file (`.gfa`) in GFA format, in the case of contig gap-filling. It contains the original contigs and the obtained gap-fill sequences (nodes of the graph), together with their overlapping relationships (arcs of the graph).
    * a log file (`.info.txt`), a tabular file with some information about the filling process for each breakpoint/grap-fill. 

	
7. **Computational resources options**
    
    Additional options are related to computational runtime and memory:
    * `-nb-cores`: number of cores to be used for computation [default '0', ie. all available cores will be used].
    * `-max-memory`: max RAM memory for the graph creation (in MBytes)  [default '2000']. Increasing the memory will speed up the graph creation phase.
    * `-max-disk`: max usable disk space for the graph creation (in MBytes)  [default '0', ie. automatically set]. Kmers are counted by writting temporary files on the disk, to speed up the counting you can increase the usable disk space.



## Details on output formats

1. Breakpoint format
    
    A breakpoint file is output by MindTheGap find and is required by MindTheGap fill. This is a plain text file in fasta format, where each insertion site (or gap to fill) corresponds to two consecutive fasta entries: left kmer and right kmer. Sequences are kmer (with k being the same value used in the de bruijn graph). MindTheGap fill will try to find a path in the de bruijn graph from the left kmer to the right one. 
    In the breakpoint file output by MindTheGap find, one can find useful information in the fasta headers, such as the genomic position of the insertion site and its genotype (detected by the homozygous or heterozygous algorithm). A typical header is as follows: 
        
        >bkpt5_chr1_pos_39114_fuzzy_0_HOM left_kmer
        #bkpt5 : this is the id of the insertion event. 
        #chr1_pos_39114 : the position of the insertion site is on chr1 at position 39114 (position just before the insertion, 1-based). 
        #fuzzy_0 : is the size of the repeated sequence at the breakpoint (here 0 means it is a clean insertion site).
        #HOM : it was detected by the homozygous algorithm.

    Note: in the case of a small repeat at the breakpoint site (fuzzy>0), the exact position can not be known inside the repeat, the reported position is always the right-most.

    Note #2: sometimes the header can contain the word `REPEATED` next to `left kmer` or `right kmer`. This concerns also repeated sequences but must not be confused with the `fuzzy` field. Fuzzy indicates if a small repeat (typically <5 bp) is exactly repeated at the breakpoint site and at an extremity of the inserted sequence (this can happen very often by chance and may have no biological meaning). Whereas "REPEATED" indicates that this insertion site is probably located in a repeated region of the reference genome, with the repeat size being >=(k-1). These breakpoints have more probability to be false positives.   
	
2. VCF variant format

    For both `.othervariants.vcf` and `insertions.vcf` files, the format follows the VCF specifications version 4.1 (see https://samtools.github.io/hts-specs/VCFv4.1.pdf). Positions are 1-based.
	
	For insertion variants, positions are left-normalized. This happens when there is a small (typically <5bp)) repeated sequence between the breakpoint site and one extremity of the inserted sequence. In this case, multiple positions are possible. Here, the leftmost position is reported and the number of possible positions is indicated in the INFO field (NPOS id). This latter value is the size of the repeated sequence + 1 (= fuzzy +1). Note that in this case the REF field not only contains the nucleotide before the insertion but also the repeated sequence (the REF field size is therefore equal to NPOS) and the ALT field contains the two copies of the repeated sequence (at both extremities). Example:
	
		chr4    618791     bkpt20  TAGG    TAGGTGTATTTAGCTCCGAGG   .       PASS    TYPE=INS;LEN=17;NPOS=4;AVK=22.71;MDK=23.00      GT      1/1
	 	#AGG is a repeat of size 3, there are 4 (NPOS) possible positions for an insertion of 17 nt from positions 618791 to 618794 on chr4
	
3. Assembled insertion format
    
    MindTheGap fill outputs a file in fasta format containing the obtained inserted sequences. Breakpoint kmers are not included in the output sequences. For each insertion breakpoint for which the filling succeeded, one can find in this file either one or several sequences with the following header:
    
        >bkpt5_chr1_pos_39114_fuzzy_0_HOM_len_59_qual_50_avg_cov_21.69_median_cov_17.00
        #same info as in the breakpoint file
        #len_59: the length in bp of the inserted sequence, here 59 bp
        #qual_50: quality of 50 (quality scores range from 0 to 50, 50 being the best quality) 
        #avg_cov_21.69: average abundance of the filled sequence (average of all its kmer abundances)
        #median_cov_17.00: median abundance of the filled sequence (median of all its kmer abundances)

    If more than one sequence are assembled for a given breakpoint, the header is as follows:
    
        >bkpt5_chr1_pos_39114_fuzzy_0_HOM_len_57_qual_15_avg_cov_21.69_median_cov_17.00 solution 2/3
        #this is the second sequence out of 3

	**Contig gap-fill header specificies**:
	
	When used with option `-contig`, the fasta header is a little bit different: it contains notably two contig identifiers (their fasta headers in the original contig file) with optionnally a suffix "_Rc" if it is reversed.
	
	 	>contig3_len_3652;contig18_len_19822_Rc;len_117_qual_50_median_cov_1350
		#contig3_len_3652: header of the source contig, contig3 in the original input file contigs.fa
		#contig18_len_19822: header of the target contig, contig18 in the original input file contigs.fa
		#_Rc: absent for the source contig and present for the target contig, this means that the end of contig3 is gap-filled with the end of contig18 (that is with the beginning of the reverse complement of contig18).
		#len_117_qual_50_median_cov_1350: information about the assembled gap-fill sequence
	 
    **Quality scores**:
    
    Each insertion is assigned a quality score ranging from 0 (low quality) to 50 (highest quality). This quality score reflects mainly repeat-associated criteria:
    * `qual=5`: if one of the breakpoint kmer could not be found exactly but with 2 errors (mismatches)
    * `qual=10`: if one of the breakpoint kmer could not be found exactly but with 1 error (mismatch)
    * `qual=15`: if multiple sequences can be assembled for a given breakpoint (note that to output multiple sequences, they must differ from each other significantly, ie. <90% id)
    * `qual=25`: if one of the breakpoint kmer is repeated in the reference genome (REPEATED field in the breakpoint file)
    * `qual=50`: otherwise.

    **Info file**:

    For each breakpoint, some informations about the filling process are given in the file `.info.txt`, whether it has been successfully filled or not. This can help understand why some breakpoints could not be filled. Here are the description of the columns:
    * column 1 : breakpoint name       
    * column 2-4 : number of nodes in the contig graph, total nt assembled, number of nodes containing the right breakpoint kmer
    * (optionnally) column 5-7 : same informations as in column 2-4 but for the filling process in the reverse direction from right to left kmer, activated only if the filling failed in the forward direction
    * last 2 columns : number of alternative filled sequences before comparison, number of output filled sequences (can be reduced if some pairs of alternative sequences are more than 90% identical).
 

## Full example

This example can be run with the small dataset in directory `data/`, for instance from the build directory:

    #find
    bin/MindTheGap find -in ../data/reads_r1.fastq,../data/reads_r2.fastq -ref ../data/reference.fasta -out example
    # 3 files are generated: 
    #   example.h5 (graph), 
    #   example.othervariants.vcf (SNPs and deletion variants), 
    #   example.breakpoints (breakpoints of insertion variants).

    #fill
    bin/MindTheGap fill -graph example.h5 -bkpt example.breakpoints -out example
    # 2 files are generated:
    #   example.insertions.fasta
    #   example.insertions.vcf

## Utility programs

Either in your bin/ directory or in ext/gatb-core/bin/, you can find additional utility programs :
* dbginfo : to get information about a graph stored in a .h5 file
* dbgh5 : to build a graph from read set(s) and obtain a .h5 file
* h5dump : to extract data stored in a .h5 file
	
## Reference

MindTheGap: integrated detection and assembly of short and long insertions. Guillaume Rizk, Anaïs Gouin, Rayan Chikhi and Claire Lemaitre. Bioinformatics 2014 30(24):3451-3457. http://bioinformatics.oxfordjournals.org/content/30/24/3451

[Web page](https://gatb.inria.fr/software/mind-the-gap/) with some updated results.
 

# Contact

To contact a developer, request help, or for any feedback on MindTheGap, please use the issue form of github: https://github.com/GATB/MindTheGap/issues

You can see all issues concerning MindTheGap [here](https://github.com/GATB/MindTheGap/issues) and GATB [here](https://www.biostars.org/t/GATB/).

If you do not have any github account, you can also send an email to claire dot lemaitre at inria dot fr
