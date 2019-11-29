# Genome assembly gap-filling using *MindTheGap*

## *MindTheGap* contig mode

In addition to the assembly of insertion variants, the `fill` module of MindTheGap can be used as a genome assembly finishing tool, to fill the gaps between a given set of contigs, wihtout any apriori on their relative order and orientation.

The basic usage of this mode is :

```bash
MindTheGap fill (-in <reads.fq> | -graph <graph.h5>) -contig <contigs.fa> [options]
```

It takes as input 2 mandatory files : the sequencing reads or their de bruijn graph if already computed (options `-in` and `-graph` respectively) and the set of contigs in fasta format (option `-contig`). 

### Specific input parameters

Most options are similar to those of the standard mode of MindTheGap, notably for the de Bruijn graph construction or for computational resource settings (see [../README.md](../README.md)). Specific options of the `contig`mode are:

- `-contig`: the contig file path in fasta format. Note that only contigs larger than $3*kmerSize$ will be used. 
  Although MindTheGap has been tested with contigs obtained with the assembler [Minia](https://github.com/GATB/minia), which uses similar assembly heuristics, contigs from any assembler may be used.
- `-overlap`: the maximal potential sequence overlap between input contigs (default = $k$). MindTheGap extract from contig extremities seed and target kmers to perform local assembly between these kmers. To ensure a non-null gap-filled sequence even between contigs that overlap, seed and target kmers are extracted at `overlap` bp from the contig extremities. In case the overlap between your contigs exceeds the 'k' value chosen for *MindTheGap*, this can be specified using the `-overlap` option.

- Local assembly limitations: local assembly may be tuned to allow larger and more complex assemblies between the contigs (than for insertion variants), with the following options:

  - `-max-nodes`: maximum number of nodes in the contig graph for each gap-filling assembly [default '100']. This arguments limits the computational time, but it can be safely set to $300$ or $1000$ in contig mode.
  - `-max-length`: maximum number of assembled nucleotides in the contig graph (nt)  [default '10000']. This arguments limits the computational time, but if gaps are large, it must be increased.

  Increasing these two parameters may improve the results for gapfilling of assemblies much shorter than their expected size.

### Output

In contig mode, *MindTheGap* returns 3 files ; 
1. GFA file : `out.gfa`
    The assembly is returned in a [GFA format graph](https://github.com/GFA-spec/GFA-spec). Both initial contigs and gapfilling sequences are represented by segments. Links indicate sequence overlaps between segments.
    
2. Insertion sequences `out.insertions.fa`
    In addition to the GFA file, gap-filling sequences are reported in a fasta format file (see the header format below).

3. Gap-filling information file :  `out.info.txt`

    For each gap-fill, some informations about the filling process are given in the file `.info.txt`, whether it has been successfully filled or not. 


### Output formats

**Gap-filling sequence header specificies**:

MindTheGap fill outputs a file in fasta format containing the obtained  gap-filling sequences (`.insertions.fasta`). Source and target kmers are not included in the output  sequences. For each pair of contigs for which the filling  succeeded, one can find in this file either one or several sequences with the following header: 

```
>contig3_len_3652;contig18_len_19822_Rc;len_117_qual_50_median_cov_1350
#contig3_len_3652: header of the source contig, contig3 in the original input file contigs.fa
#contig18_len_19822: header of the target contig, contig18 in the original input file contigs.fa
#_Rc: absent for the source contig and present for the target contig, this means that the end of contig3 is gap-filled with the end of contig18 (that is with the beginning of the reverse complement of contig18).
#len_117_qual_50_median_cov_1350: information about the assembled gap-fill sequence, median_ cov giving the median kmer abundance of the sequence in the sample.
```

it contains notably two contig identifiers (their fasta headers in the original contig file) with optionnally a suffix "_Rc" if it is reversed.

**Gap-filling information file**

For each gap-fill, some informations about the filling process are given in the file `.info.txt`, whether it has been successfully filled or not. This can help  understand why some gaps could not be filled. Here are the  description of the columns:

- column 1 : gap-filling name
- column 2-4 : number of nodes in the contig graph, total nt assembled, number of nodes containing the right breakpoint kmer
- (optionnally) column 5-7 : same informations as in column 2-4 but  for the filling process in the reverse direction from right to left  kmer, activated only if the filling failed in the forward direction
- last 2 columns : number of alternative filled sequences before  comparison, number of output filled sequences (can be reduced if some  pairs of alternative sequences are more than 90% identical).

### Dealing and analysing genome graphs (GFA files)

The graph output by MindTheGap can be easily visualized using [Bandage](https://github.com/rrwick/Bandage). 

Note that GFA Graphs supplied by MindTheGap may contain redundant sequence information (for instance this is likely that two contigs are linked in the graph by two gapfillings with reverse-complement sequences). Before further analyses, we recommend to simplify and reduce the redundancy in the graph using the scripts available in (MinYS github repository)[https://github.com/cguyomar/MinYS] :

```
git clone https://github.com/cguyomar/MinYS.git
python3 MinYS/graph_simplification/graph_simplification.py MindTheGap_output.gfa simplified_graph.gfa
```

Other usefull scripts are available in the (MinYS github repository)[https://github.com/cguyomar/MinYS]  to deal with this graph data structure : to enumerate paths, to convert the segments to a fasta file, to filter connected components according to their size...



## *MinYS*: targeted assembly pipeline

*MindTheGap* in contig-mode is an essential step of the targeted assembly tool MinYS which is freely available here: [https://github.com/cguyomar/MinYS](https://github.com/cguyomar/MinYS).

MinYS stands for *Mine Your Symbiont* and was designed for de novo assembly of a bacterial genome of interest, sequenced in a metagenomic context, and with the help of a (potentially distant) reference genome. A typical situation when studying symbiont genomes within their eukaryotic host sequencing. It consists in three steps : 

- Recruiting reads by mapping onto the reference genome,
- Assembly of those reads in *backbone* contigs,
- Gap-filling of these *backbone* contigs with the whole readset.




