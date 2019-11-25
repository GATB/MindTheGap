# Genome assembly gap-filling using *MindTheGap*

## *MindTheGap* contig mode

In addition to the assembly of structural variations, the fill module of MindTheGap can be used as a genome assembly finishing tool.

The basic usage of this mode is :

```bash
MindTheGap fill (-in <reads.fq> | -graph <graph.h5>) -contig <contigs.fa> [options]
```

Options are similar to those of the standard mode of MindTheGap.
`-contig` is a fasta file containing the contigs to use.
Although MindTheGap has been tested with contigs obtained with the assembler [Minia](https://github.com/GATB/minia), which uses similar assembly heuristics, contigs from any assembler may be used.

### Output

In contig mode, *MindTheGap* returns 3 files ; 
1. GFA file
    The assembly is returned in a [GFA format graph](https://github.com/GFA-spec/GFA-spec).
    Both initial contigs and gapfilling sequences are represented by segments. Links indicate an overlap between segments.
    
    GFA Graphs supplied by MindTheGap may contain redundant sequences.
    Before further analyses, it should be simplified using the script available in `pipeline/genome/graph/graph_simplification.py`
    Its usage is `graph_simplification.py MindTheGap_output.gfa simplified_graph.gfa`

    Afterwards, the graph can be further analyzed using [Bandage](https://github.com/rrwick/Bandage) 

2. Info file
    This file is supplied as `out.info.txt`
    It is a tab delimited file, with one by seed kmer used during the gapfilling process (and therefore two by input contigs)
    Column 1 is the seed identifier, made of the contig name and an eventual "_Rc" suffix if it is the reverse complemented seed.
    Column 2 and 3 are the number of nodes and total length of the graph built during the local assembly process
    Column 4 is the number of nodes containing a target kmer.
    Columns 5 and 6 are the number of solutions found, before and after comparison of the sequences under a 90% identity threshold.

3. Insertion sequences
    In addition to the GFA file, gap-filling sequences are reported in the fasta format in the `out.insertions.fa` file.


### Additional input parameters :

#### Contig overlap

In many assembly outputs, contigs ends may overlap.
In particular, contigs from *De Bruijn* based assemblies may overlap from `k`.
In case the overlap between your contigs exceeds the 'k' value chosen for *MindTheGap*, we recommend specifying the overlap using the `-overlap` option.

#### Graph complexity

Local assembly may be tuned by allowing larger and more complex assemblies between the contigs.
Option `-max-length` specifies the maximum length a gapfilling may reach, while option `-max-nodes` is the number of nodes that can be built in the assembly graph.
Increasing these two parameters may improve the results for gapfilling of assemblies much shorter than their expected size.

### Output formats

**Gap-filling sequence header specificies**:

MindTheGap fill outputs a file in fasta format containing the obtained  gap-filling sequences (`.insertions.fasta`). Source and target kmers are not included in the output  sequences. For each pair of contigs for which the filling  succeeded, one can find in this file either one or several sequences with the following header: 

```
>contig3_len_3652;contig18_len_19822_Rc;len_117_qual_50_median_cov_1350
#contig3_len_3652: header of the source contig, contig3 in the original input file contigs.fa
#contig18_len_19822: header of the target contig, contig18 in the original input file contigs.fa
#_Rc: absent for the source contig and present for the target contig, this means that the end of contig3 is gap-filled with the end of contig18 (that is with the beginning of the reverse complement of contig18).
#len_117_qual_50_median_cov_1350: information about the assembled gap-fill sequence
```

it contains notably two contig identifiers (their fasta headers in the original contig file) with optionnally a suffix "_Rc" if it is reversed.



## *MinYS*: targeted assembly pipeline

*MindTheGap* in contig-mode is an essential step of the targeted assembly tool MinYS which is freely available here: [https://github.com/cguyomar/MinYS](https://github.com/cguyomar/MinYS).

MinYS stands for *Mine Your Symbiont* and was designed for de novo assembly of a bacterial genome of interest, sequenced in a metagenomic context, and with the help of a (potentially distant) reference genome. A typical situation when studying symbiont genomes within their eukaryotic host sequencing. It consists in three steps : 

- Recruiting reads by mapping onto the reference genome,
- Assembly of those reads in *backbone* contigs,
- Gap-filling of these *backbone* contigs with the whole readset.




