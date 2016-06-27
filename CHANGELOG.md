# Change Log

--------------------------------------------------------------------------------
## [Unreleased]

--------------------------------------------------------------------------------
## [2.0.0] - 2016-06-28 

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
