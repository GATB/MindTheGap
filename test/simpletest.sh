#!/bin/bash

run_test()
{
    # param : reads_file ref_file true_result prefix
    ../build/MindTheGap find -in $1 -ref $2 -kmer-size 31 -out output/$4_find $5 1> output/$4_find.out 2> output/$4_find.err

    ../build/MindTheGap fill -bkpt output/$4_find.breakpoints -graph output/$4_find.h5  -out output/$4_fill 1> output/$4_fill.out 2> output/$4_fill.err

    diff --ignore-matching-lines=">" output/$4_fill.insertions $3 1> /dev/null 2>&1

    var=$?
    if [ $var -eq 0 ]
    then
	echo  PASSED
    else
	echo  FAILED
    fi
}

mkdir -p output

echo  "Testing canonical k-1 insert site"

run_test reads/master.fasta references/deleted.fasta truths/insertion.fasta k-1 "-insert-only -no-backup -homo-only"

echo "Testing solo SNP"

run_test reads/master.fasta references/sSNP.fasta truths/sSNP.fasta sSNP "-snp-only -no-backup"

echo "Testing duo SNP"

run_test reads/master.fasta references/dSNP.fasta truths/dSNP.fasta dSNP "-snp-only -no-backup"

echo "Testing heterozygote"

run_test reads/deleted.fasta,reads/master.fasta references/deleted.fasta truths/insertion.fasta hete "-no-backup -max-rep 2"

echo "Testing witht N in reference :"
echo "    N in stretch :"

run_test reads/master.fasta references/n_in_stretch.fasta truths/n_in_stretch.fasta n_in_stretch

echo "    N in before gap :"

run_test reads/master.fasta references/n_before_gap.fasta truths/n_before_gap.fasta n_before_gap


echo "    N after gap :"

run_test reads/master.fasta references/n_after_gap.fasta truths/n_after_gap.fasta n_after_gap

