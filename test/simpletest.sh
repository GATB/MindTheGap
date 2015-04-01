#!/bin/bash



#insert should be TCGTTCGGTGACCGGCCATAGGA


echo  "Testing canonical k-1 insert site"


../build/MindTheGap find -in reads1.fasta.gz -ref ref1.fasta -kmer-size 31 -abundance-min 4 -out insert_simple   &> /dev/null

../build/MindTheGap fill -bkpt insert_simple.breakpoints -graph insert_simple.h5
#../build/Mindthegap fill insert_simple.breakpoints -p insert_simple  &> /dev/null


sed  '/^>/d' insert_simple.insertions.fa > tres
diff tres true_res1

var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
else
    echo  FAILED
fi

