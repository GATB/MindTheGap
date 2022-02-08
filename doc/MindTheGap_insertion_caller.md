# MindTheGap  for genomic insertion variant calling

MindTheGap  performs detection and assembly of **DNA insertion variants** in short read datasets with respect to a reference genome. It is designed to call insertions of any size, whether they are novel or duplicated, homozygous or heterozygous in the donor genome.

It takes as input a set of reads and a reference genome. It produces two sets of FASTA sequences: one is the set of breakpoints of detected insertion sites, the other is the set of assembled insertions for each breakpoint. For each breakpoint, MindTheGap either returns a single insertion sequence (when there is no assembly ambiguity), or a set of candidate insertion sequences (due to ambiguities) or nothing at all (when the insertion is too complex to be assembled).  Its final output is a VCF file, giving for each insertion variant, its insertion site location on the reference genome, one or several candidate insertion sequences, and its genotype in the sample.

Since version 2.0.0, MindTheGap can detect other types of variants, not only insertion events. These are homozygous SNPs and homozygous deletions of any size. They are detected by the find module of MindTheGap and are output separately in a VCF file. Importantly, even if the user is not interested in these types of variants, it is worth to detect them  since it can improve the recall of the insertion event detection algorithm : it is now possible to find insertion events that are located at less than k nucleotides from an other such variant.


For more details on the method and some recent results, see the [web page](http://gatb.inria.fr/software/mind-the-gap/).
	
## Usage

MindTheGap is composed of two main modules : breakpoint detection (`find` module) and the local assembly of insertions (`fill` module). Both steps are implemented in a single executable, MindTheGap, and can be run independently by specifying the module name as follows :

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

2. **Common options**

   Common options for input read files, de Bruijn graph construction and computational resource settings are detailed in the main [README.md](../README.md).
   
4. **Find module specific options**
  
    In addition to the read or graph files, the `find` module has one mandatory option `-ref` and several optional options:
    
    * `-ref`: the path to the reference genome file (in fasta format).
    * `-homo-only`: only homozygous insertions are reported (default: not activated).
    * `-max-rep`: maximal repeat size allowed for fuzzy sites  [default '5']. 
    * `-het-max-occ`: maximal number of occurrences of a (k-1)mer in the reference genome allowed for heterozyguous insertion breakpoints  [default '1']. In order to detect an heterozyguous insertion breakpoints, both flanking k-1-mers, at each side of the insertion site, must have strictly less than this number of occurrences in the reference genome. This prevents false positive predictions inside repeated regions. Warning : increasing this parameter may lead to numerous false positives (genomic approximate repeats).
    * `-bed`: the path to a bed file defining genomic regions, to limit the find algorithm to particular regions of the genome. This can be usefull for exome data.
    
5. **Fill module specific options**
  
    In addition to the read or graph files, the `fill` module has one other mandatory option, `-bkpt`:
     	
    * `-bkpt`: the breakpoint file path. This is one of the output of the `find` module and contains for each detected insertion site its left and right kmers from and to which the local assembly will be performed (see section E for details about the format).
	
	The fill module has several optional options:
	
	* `-max-nodes`: maximum number of nodes in contig graph for each insertion assembly [default '100']. This arguments limits the computational time, this is especially useful for complex genomes.
    * `-max-length`: maximum number of assembled nucleotides in the contig graph (nt)  [default '10000']. This arguments limits the computational time, this is especially useful for complex genomes.
    * `-filter`: if set, insertions with multiple solutions are not output in the final vcf file (default : not activated).
	* `-fwd-only`: if set, inserted sequences are searched in only one direction : from the left kmer to the right kmer. Default behavior : not set, ie. when no solution is found in the forward direction, a solution is searched in the opposite direction (revcomp(right) --> revcomp(left)).
	
6. **MindTheGap Output**
  
    All the output files are prefixed either by a default name: "MindTheGap_Expe-[date:YY:MM:DD-HH:mm]" or by a user defined prefix (option `-out` of MindTheGap)
    Both MindTheGap modules generate the graph file if reads were given as input: 
    
    * a graph file (`.h5`). This is a binary file, to obtain information stored in it, you can use the utility program dbginfo located in your bin directory or in ext/gatb-core/bin/.
    
    `MindTheGap find` generates the following output files:
    
    * a breakpoint file (`.breakpoints`) in fasta format. It contains the breakpoint sequences of each detected insertion site. Each insertion site corresponds to 2 consecutive entries in the fasta file : sequences are the left and right side flanking kmers.
    * a variant file (`.othervariants.vcf`) in vcf format. It contains SNPs, deletion and very small insertions (1-2 bp).
    
    `MindTheGap fill` generates the following output files:
    
    * a sequence file (`.insertions.fasta`) in fasta format (for insertions >2 bp). It contains the inserted sequences or contig gap-fills that were successfully assembled. In the case of insertion variants, the location of each insertion on the reference genome can be found in its fasta header. The fasta header includes also information about each gap-fill such as its length, quality score and median kmer abundance.
    * an insertion variant file (`.insertions.vcf`) in vcf format (for insertions >2 bp). This file contains all information of assembled insertion variants as in the `.insertions.fasta` file but in a different format. Here, insertion site positions are 1-based and left-normalized according to the VCF format specifications (contrary to positions indicated in the `.breakpoints` and `insertions.fasta` files which are right-normalized). Normalization occurs when multiple positions are possible for a single variation due to a small repeat. 
	* a log file (`.info.txt`), a tabular file with some information about the filling process for each breakpoint/grap-fill. 
    



## Details on output formats

1. Breakpoint format
  
    A breakpoint file is output by MindTheGap find and is required by MindTheGap fill. This is a plain text file in fasta format, where each insertion site (or gap to fill) corresponds to two consecutive fasta entries: left kmer and right kmer. Sequences are kmer (with k being the same value used in the de bruijn graph). MindTheGap fill will try to find a path in the de bruijn graph from the left kmer to the right one. 
    In the breakpoint file output by MindTheGap find, one can find useful information in the fasta headers, such as the genomic position of the insertion site and its genotype (detected by the homozygous or heterozygous algorithm). A typical header is as follows: 
        
        >bkpt5_chr1_pos_39114_fuzzy_0_HOM left_kmer
        #bkpt5 : this is the id of the insertion event. 
        #chr1_pos_39114 : the position of the insertion site is on chr1 at position 39114 (position just before the insertion, 1-based). 
        #fuzzy_0 : is the size of the repeated sequence at the breakpoint (here 0 means it is a clean insertion site).
        #HOM : it was detected by the homozygous algorithm.

    Note: in the case of a small repeat at the breakpoint site (fuzzy>0), the exact position can not be known inside the repeat, the reported position here is always the right-most (contrary to the VCF output of the `fill` module which is left-normalized).

    Note #2: some times the header can contain the word `REPEATED` next to `left kmer` or `right kmer`. This concerns also repeated sequences but must not be confused with the `fuzzy` field. Fuzzy indicates if a small repeat (typically <5 bp) is exactly repeated at the breakpoint site and at an extremity of the inserted sequence (this can happen very often by chance and may have no biological meaning). Whereas "REPEATED" indicates that this insertion site is probably located in a repeated region of the reference genome, with the repeat size being >=(k-1). These breakpoints have more probability to be false positives.   
	
2. VCF variant format

    For both `.othervariants.vcf` and `insertions.vcf` files, the format follows the VCF specifications version 4.1 (see https://samtools.github.io/hts-specs/VCFv4.1.pdf). Positions are 1-based.
	
	For insertion variants (`insertions.vcf` file only), positions are **left-normalized**. This happens when there is a small (typically <5bp)) repeated sequence between the breakpoint site and one extremity of the inserted sequence. In this case, multiple positions are possible. Here, the leftmost position is reported and the number of possible positions is indicated in the INFO field (NPOS id). This latter value is at least the size of the repeated sequence + 1 (>= fuzzy +1). 
	
		chr4    618791     bkpt20  T    TAGGTGTATTTAGCTCCG   .       PASS    TYPE=INS;LEN=17;QUAL=50;NSOL=1;NPOS=4;AVK=22.71;MDK=23.00      GT      1/1
		#there are 4 (NPOS) possible positions for an insertion of 17 nt from positions 618791 to 618794 on chr4, therefore the repeat is of size 3 and is AGG.
	
	FILTER field: can be `PASS`or `LOWQUAL` (for insertions with multiple solutions)
	
	INFO fields:  
	
	* `TYPE`: variant type, INS for insertion
	* `LEN`: insertion size in bp
	* `QUAL`: quality of the insertion (quality scores range from 0 to 50, 50 being the best quality, see the different quality scores [below](#quality))
	* `NSOL`: number of alternative sequences that were assembled at this position (note that to output multiple sequences, they must differ from each other significantly, ie. <90% identity)
	* `NPOS`: number of possible positions where the insertion event can occur giving the same ALT sequence (see the left-normalization paragraph).
	* `AVK`: average abundance of the inserted sequence (average value of the abundances of all its overlapping kmers)
	* `MDK`: median abundance of the inserted sequence (median value of the abundances of all its overlapping kmers)
	
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

	
	
4. Assembled insertion quality scores:
<a name="quality"></a>
  
    Each insertion is assigned a quality score ranging from 0 (low quality) to 50 (highest quality). This quality score reflects mainly repeat-associated criteria:
    
    * `qual=5`: if one of the breakpoint kmer could not be found exactly but with 2 errors (mismatches)
    * `qual=10`: if one of the breakpoint kmer could not be found exactly but with 1 error (mismatch)
    * `qual=15`: if multiple sequences can be assembled for a given breakpoint (note that to output multiple sequences, they must differ from each other significantly, ie. <90% identity)
    * `qual=25`: if one of the breakpoint kmer is repeated in the reference genome (REPEATED field in the breakpoint file)
    * `qual=50`: otherwise.

5. Gap-filling information file:

    For each gap-fill, some informations about the filling process are given in the file `.info.txt`, whether it has been successfully filled or not. This can help understand why some breakpoints could not be filled. Here are the description of the columns:
    
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
    # 3 files are generated:
    #   example.insertions.fasta
    #   example.insertions.vcf
    #   example.info.txt
