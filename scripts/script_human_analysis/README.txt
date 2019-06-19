Two scripts are available to improve performance of MindTheGap for human genome analysis :

Inser_snp_in_ref.py :
It allows user to integrate SNP called from GATK HaplotypeCaller in a reference genome.
Three paramaters are required : -s GATK.vcf, -g reference_genome.fa, -o altered_genome.fa

Context_genome.py :
It allows user to filter potential false positive.
The script will check k-mer connectivity around each breakpoints.
By default, if  more than 20% of the last 50 k-mers contain unusual connectivity (number of branching k-mer for a k-mer is greater than 2) the breakpoints is not kept.
Four parameters are required :
-g MindTheGap_file.h5
-p Reference_genome.fa
-b Breakpoint_file.breakpoints
-o Breakpoints_filtered.breakpoints

Use -m to set a specific threshold of connectivity (0 to 1)

Example of running pipeline :
python3.5 /MindTheGap/script/python3/Inser_snp_in_ref.py -g genome.fa -s GATKHC.vcf -o altered_genome.fa
/MindTheGap/build/bin/MindTheGap find -ref altered_genome.fa -in part1.fastq.gz,part2.fastq.gz  -abundance-min auto -out OUTPUT_FIND
python3.5 /MindTheGap/script/python3/Context_genome_WG.py -g OUTPUT_FIND.h5 -p altered_genome.fa -b OUTPUT_FIND.breakpoints -o OUTPUT_FIND_filter.breakpoints
/MindTheGap/build/bin/MindTheGap fill -graph OUTPUT_FIND.h5 -bkpt OUTPUT_FIND_filter.breakpoints -out OUTPUT_FIND_filter -filter
