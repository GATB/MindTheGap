#!/bin/bash

# look for MindTheGap binary. In devel mode, it's in ../build/bin directory.
# In production mode, it's in ../bin directory.
if [ -f "../bin/MindTheGap" ]
then
 bindir="../bin"
elif [ -f "../build/bin/MindTheGap" ]
then
 bindir="../build/bin"
else
 echo "could not find a compiled Bloocoo binary"
 exit 1
fi

testOK="true"

run_test()
{
    # param : reads_file ref_file true_result prefix
    $bindir/MindTheGap find -in $1 -ref $2 -kmer-size 31 -out output/$4_find $5 1> output/$4_find.out 2> output/$4_find.err

    $bindir/MindTheGap fill -bkpt output/$4_find.breakpoints -graph output/$4_find.h5  -out output/$4_fill 1> output/$4_fill.out 2> output/$4_fill.err

    diff --ignore-matching-lines=">" output/$4_fill.insertions $3 1> /dev/null 2>&1

    var=$?
    if [ $var -eq 0 ]
    then
	     eval $6="passed"
    else
	     testOK="false" ; eval $6="FAILED"
    fi
}

mkdir -p output
output=""

output=$output"clean-insert : "
run_test reads/master.fasta references/deleted.fasta truths/insertion.fasta k-1 "-insert-only" retvalue
output=${output}${retvalue}

#output=$output"\n1-SNP : "
#run_test reads/master.fasta references/sSNP.fasta truths/sSNP.fasta sSNP "-snp-only" retvalue
#output=${output}${retvalue}

#output=$output"\n3-SNP*2 : "
#run_test reads/master.fasta references/multiSNP.fasta truths/multiSNP.fasta multiSNP "-snp-only" retvalue
#output=${output}${retvalue}

output=$output"\nsnp-before-clean-insert : "
run_test reads/master.fasta references/deleted_before_SNP.fasta truths/insertion_before_SNP.fasta k-1_before_SNP "-no-deletion -homo-only" retvalue
output=${output}${retvalue}

#output=$output"\nsnp-begin-fuzzy : "
#run_test reads/beginfuzzySNP.fasta references/beginfuzzySNP.fasta truths/beginfuzzySNP.fasta beginfuzzySNP "-snp-only" retvalue
#output=${output}${retvalue}

output=$output"\nhetero-insert : "
run_test reads/deleted.fasta,reads/master.fasta references/deleted.fasta truths/insertion.fasta hete "-hete-only -max-rep 2" retvalue
output=${output}${retvalue}

output=$output"\ndeletion : "
run_test reads/deleted.fasta references/master.fasta truths/deletion.fasta deletion "-deletion-only" retvalue
output=${output}${retvalue}

output=$output"\nfuzzy-deletion : "
run_test reads/deletionfuzzy.fasta references/deletionfuzzy.fasta truths/deletionfuzzy.fasta deletionfuzzy "-deletion-only" retvalue
output=${output}${retvalue}

output=$output"\nn-in-solid-stretch : "
run_test reads/master.fasta references/n_in_stretch.fasta truths/n_in_stretch.fasta n_in_stretch "-insert-only" retvalue
output=${output}${retvalue}

output=$output"\nn-in-before-clean-insert : "
run_test reads/master.fasta references/n_before_gap.fasta truths/n_before_gap.fasta n_before_gap "-insert-only" retvalue
output=${output}${retvalue}

output=$output"\nn-after-clean-insert : "
run_test reads/master.fasta references/n_after_gap.fasta truths/n_after_gap.fasta n_after_gap "-insert-only" retvalue
output=${output}${retvalue}

echo -e $output | column -t
echo "test 3: $testOK"

if [ $testOK == "false" ]; then
  exit 1
fi
