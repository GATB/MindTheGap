#!/bin/bash

run_test()
{
    # param : reads_file ref_file true_result prefix
    ../build/MindTheGap find -in $1 -ref $2 -kmer-size 31 -abundance-min 4 -out output/$4_find 1> output/$4_find.out 2> output/$4_find.err

    ../build/MindTheGap fill -bkpt output/$4_find.breakpoints -graph output/$4_find.h5  -out output/$4_fill 1> output/$4_fill.out 2> output/$4_fill.err

    diff --ignore-matching-lines=">" output/$4_fill.insertions $3

    var=$?
    if [ $var -eq 0 ]
    then
	echo  PASSED
    else
	echo  FAILED
    fi
}

echo  "Testing canonical k-1 insert site"

run_test reads/master.fasta references/deleted.fasta truths/insertion.fasta k-1

echo "Testing solo SNP"

run_test reads/master.fasta references/sSNP.fasta truths/insertion.fasta sSNP

echo "Testing heterozygote"

run_test reads/deleted.fasta,reads/master.fasta references/deleted.fasta truths/insertion.fasta hete
