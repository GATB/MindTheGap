#!/bin/bash

run_test()
{
    # param : reads_file ref_file true_result prefix
    ../build/MindTheGap find -in $1 -ref $2 -kmer-size 31 -abundance-min 4 -out $4_find   &> /dev/null

    ../build/MindTheGap fill -bkpt $4_find.breakpoints -graph $4_find.h5  -out $4_fill &> /dev/null

    sed  '/^>/d' $4_fill.insertions > $4.res
    diff $4.res $3

    var=$?
    if [ $var -eq 0 ]
    then
	echo  PASSED
    else
	echo  FAILED
    fi
}

#insert should be TCGTTCGGTGACCGGCCATAGGA
echo  "Testing canonical k-1 insert site"

run_test reads1.fasta.gz ref1.fasta true_res1 k-1

#SNP at pos 40 A in G before the insert
echo "Testing SNP before k-1 insert site"

run_test reads_bSNP.fasta.gz ref1.fasta true_resbSNP bSNP

#SNP at pos 172 A in C
echo "Testing solo SNP"

run_test reads_sSNP.fasta.gz ref1.fasta true_ressSNP sSNP
