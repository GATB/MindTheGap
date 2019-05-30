# Reference guided genome assembly using *MindTheGap*

## *MindTheGap* contig mode

In addition to the assembly of structural variations, the fill module of MindTheGap can be used as a genome assembly finishing tool.
The basic usage of this mode is :

```bash
MindTheGap fill (-in <reads.fq> | -graph <graph.h5>) -contig <contigs.fa>) [options]
```

Options are similar to those of the standard mode of MindTheGap.
`-contig` is a fasta file containing the contigs to use.
Although MindTheGap has been tested with contigs from Minia, which uses similar assembly heuristics, contigs from any assembler may be used.

### Output

In contig mode, *MindTheGap* returns 3 files ; 
1. GFA output
    The assembly is returned in a [GFA format graph](https://github.com/GFA-spec/GFA-spec).
    Both initial contigs and gapfilling sequences are represented by segments. Links indicate an overlap between segments.
    
    GFA Graph supplied by MindTheGap may contain redundant sequences.
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

3. Insertions sequences
    In addition to the GFA file, insertion sequences are reported in the fasta format in the `out.insertions.fa` file.


### Additional parameters :

#### Contig overlap

In many assembly outputs, contigs ends may overlap.
In particular, contigs from *De Bruijn* based assemblies may overlap from `k`.
In case the overlap between your contigs exceeds the 'k' value chosen for *MindTheGap*, we recommand specifying the overlap using the `-overlap` option.

#### Graph complexity

Local assembly may be tuned by allowing larger and more complex assemblies between the contigs.
Option `-max-length` specifies the maximum length a gapfilling may reach, while option `-max-nodes` is the number of nodes that can be built in the assembly graph.
Increasing these two parameters may improve the results for gapfilling of assemblies much shorter than their expected size.

## *MindTheGap* assembly pipeline

Along with *MindTheGap* is distributed a pipeline enabling reference guided genome assembly.
It consists in three steps : 
- Recruiting reads by mapping onto the reference genome.
- Assembly of those reads.
- Gapfilling of the contigs with the whole readset.

This pipeline is available in `pipeline/mtg_pipeline.py`.

### Requirements

- MindTheGap
- [BWA](http://bio-bwa.sourceforge.net/) (read mapping)
- [Minia](https://github.com/GATB/minia) (contig assembly)
- Biopython (graph simplification)
- [Bandage](https://github.com/rrwick/Bandage) (Optionnal, for assembly graph visualization) 

### Usage

```
[main options]:
  -in                   (1 arg) :    input reads file
  -1                    (1 arg) :    input reads first file
  -2                    (1 arg) :    input reads second file
  -fof                  (1 arg) :    input file of read files
  -out                  (1 arg) :    output directory for result files [Default: ./mtg_results]

[mapping options]:
  -ref                  (1 arg) :    bwa index

[assembly options]:
  -minia-bin            (1 arg) :    path to Minia binary
  -assembly-kmer-size   (1 arg) :    kmer size used for Minia assembly (should be given even if bypassing minia assembly step, usefull knowledge for gap-filling) [Default: 31]
  -assembly-abundance-min 
                        (1 arg) :    Minimal abundance of kmers used for assembly [Default: auto]
  -min-contig-size      (1 arg) :    minimal size for a contig to be used in gapfilling [Default: 0]

[gapfilling options]:
  -mtg-dir              (1 arg) :    path to MindTheGap build directory
  -gapfilling-kmer-size 
                        (1 arg) :    kmer size used for gapfilling [Default: 31]
  -gapfilling-abundance-min 
                        (1 arg) :    Minimal abundance of kmers used for gapfilling [Default: auto]
  -max-nodes            (1 arg) :    Maximum number of nodes in contig graph [Default: 100]
  -max-length           (1 arg) :    Maximum length of gapfilling (nt) [Default: 10000]

[continue options]:
  -contigs              (1 arg) :    Contigs in fasta format - override mapping and assembly
  -graph                (1 arg) :    Graph in h5 format - override graph creation

[core options]:
  -nb-cores             (1 arg) :    number of cores [Default: 0]
```
- If *minia* of *MindTheGap* are not in $PATH, a path to the minia binary of MindTheGap build directory has to be supplied using `-minia-bin` or `-mtg-dir`
- `-contigs` and `-graph` may be used to bypass the mapping/assembly step, or the graph creation. 
    In the first case, `-assembly-kmer-size` should be supplied as the overlap between contigs.




