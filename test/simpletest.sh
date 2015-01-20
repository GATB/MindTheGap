#!/bin/bash



#insert should be TCGTTCGGTGACCGGCCATAGGA


echo  "Testing canonical k-1 insert site"

../build/Mindthegap index reads.fasta.gz  -g 1000 -p insert_simple -k 31 -t 4  &> /dev/null


../build/Mindthegap find ref.fasta -p  insert_simple   &> /dev/null


../build/Mindthegap fill insert_simple.breakpoints -p insert_simple  &> /dev/null


sed  '/^>/d' insert_simple.insertions.fa > tres
diff tres true_res1

var=$?
if [ $var -eq 0 ] 
then
    echo  PASSED
else
    echo  FAILED
fi

