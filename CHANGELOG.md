# Change Log

--------------------------------------------------------------------------------
## [Unreleased]

--------------------------------------------------------------------------------
## [2.2.2] - 2020-06-19

* A bug fix: updating gatb-core version, notably this fixes a bug in the `fill` module: nodes at extremities of contigs of size exactly `k` were not marked correctly, potentially leading to duplicated contigs in the contig graph. This could prevent exploring some parts of the graph, and if graph exploration parameters where set too large (`-max-nodes` and `-max-length`), it could lead in some rare cases to extreme running times and/or memory consumptions. This should no longer happen now.

--------------------------------------------------------------------------------
## [2.2.1] - 2019-11-29

* Some bug fixes:
     * updating gatb-core version, notably this fixes a potentially important bug, where the de bruijn graph was erroneous for large datasets (such as human re-sequencing ones)
     * bug fix in the fill module, node marking in contig graph construction was not working properly leading to obtain too many solutions, in the case of multiple solutions
* Some improvements:
    * optimization of the algorithm to find paths in the  fill module (should be faster);
    * new options for insertion variant detection:
        * `-bed` (find module): to limit the search of insertion breakpoint in specific regions.
        * `-filter` (fill module): to remove insertions with multiple sequence solutions from the final vcf file (since they most often are a sign of false positive)

Note that the targeted assembly pipeline, including the gfa graph simplification scripts, is no longer included in this repository. This is now a proper tool, called MinYS (for MineYourSymbiont), which is distributed independently of MindTheGap and has now its own github repository : [https://github.com/cguyomar/MinYS](https://github.com/cguyomar/MinYS)
    

--------------------------------------------------------------------------------
## [2.2.0] - 2018-07-06

* A nice novel feature: insertion variants are now output in vcf format! and with left-normalization (ie. if several equivalent positions are possible for a given insertion event, the left-most is output and the size of the ambiguity is indicated).
* Some improvements and bug fixes:
	* faster graph loading in fill module;
	* if multiple inserted sequence solutions, better handling of very similar ones, the number of output solutions can be reduced;
	* better handling of N stretches in the reference genome, resulting in less False Positive calls in find module;
	* better recall for very small heterozygous insertion variants (bug fix when the insertion is size smaller or equal than the ambiguity size).
	* a CI simple test for the Fill module with option `-contig`.
	

--------------------------------------------------------------------------------
## [2.1.0] - 2018-06-13

A nice novel feature:

MindTheGap can now also be used as a genome assembly finishing tool: it can fill the gaps between a set of input contigs without any a priori on their relative order and orientation. This new feature is available in the Fill module with option `-contig`.

Some bug and compilation fixes, by updating the gatb-core version to 1.4.1 and more.


--------------------------------------------------------------------------------
## [2.0.2] - 2017-07-06

Some new features:
* the Fill module is now parallelized and can use several cores
* additional information is output by the fill module:
	* the abundance of each filled sequence is now computed and written in the fasta file
	* a log file is output giving details about each gap-filling process

Bug fix:
* some gap-filled sequences were incorrect (this happened only for multiple filled sequences in rare cases)

--------------------------------------------------------------------------------
## [2.0.1] - 2016-07-21

This is a bug-fix release :
* fixed a compilation issue with old version of clang compilers (prior to clang 4.3 on mac), by updating the gatb-core version to 1.2.2.

--------------------------------------------------------------------------------
## [2.0.0] - 2016-06-29

*   Initial release after refactoring the whole code of MindTheGap to use the **GATB library**.
    Some of the benefits:
    * faster (thanks to GATB improvements in kmer counting!);
    * no longer need to recompile for changing the `k` parameter;
    * automatic estimation of the paramater `abundance-min`;
    * more user-friendly usage, with more readable help, progress bars, input-output summaries, etc.
    * compatibility with other GATB tools: input-output graph in `h5` format.
*   **New features** (with respect to the published version, August 2014):
    * detection of homozygous SNPs and deletions (output in a separate VCF file); this should also improve the recall of insertion event detection;
    * a quality score is now associated to each insertion prediction, this enables to filter out some predictions and to obtain a high-confidence subset.

Have a look at a comparison between the published and the 2.0.0 versions on simulated data [here](https://gatb.inria.fr/mindthegap-insertion-event-detection/).
