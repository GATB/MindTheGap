#! /bin/bash

# look for MindTheGap binary. In devel mode, it's in ../build/bin directory.
# In production mode, it's in ../bin directory.
if [ -f "../bin/MindTheGap" ]
then
 bindir="../bin"
elif [ -f "../build/bin/MindTheGap" ]
then
 bindir="../build/bin"
else
 echo "could not find a compiled MindTheGap binary"
 exit 1
fi

RETVAL=0
testDir="test-output"
outputPrefix=$testDir/full-test
goldPrefix="full_test/gold"

outputPrefix2=$testDir/contig-test
goldPrefix2="contig_test/gold"


## First cleaning test dir
if [ -d  $testDir ]; then
  rm -rf $testDir
fi
mkdir $testDir



################################################################################
# we launch the find module
################################################################################
${bindir}/MindTheGap find -in ../data/reads_r1.fastq,../data/reads_r2.fastq -ref ../data/reference.fasta -out $outputPrefix >$outputPrefix.out -nb-cores 1 2> /dev/null

################################################################################
# we check the results 
################################################################################

# Checking the .othervariants.vcf :
sh compare_vcf.sh $outputPrefix.othervariants.vcf $goldPrefix.othervariants.vcf 1> /dev/null 2>&1
var=$?

if [ $var -eq 0 ]
then
echo "full-test find vcf         : PASS"
else
echo "full-test find vcf         : FAILED"
RETVAL=1
fi

# Checking the .breakpoints :
#diff --ignore-matching-lines=">" $outputPrefix.breakpoints $goldPrefix.breakpoints 1> /dev/null 2>&1
#var=$?

tmp1=$outputPrefix.breakpoints.tmp
tmp2=$testDir/tmp2

grep -v "^>" $outputPrefix.breakpoints > $tmp1
grep -v "^>" $goldPrefix.breakpoints > $tmp2


diff $tmp1 $tmp2 1> /dev/null 2>&1
var=$?


if [ $var -eq 0 ]
then
echo "full-test find breakpoints : PASS"
else
echo "full-test find breakpoints : FAILED"
RETVAL=1
fi

################################################################################
# we launch the find module with bed option
################################################################################
${bindir}/MindTheGap find -in ../data/reads_r1.fastq,../data/reads_r2.fastq -ref ../data/reference.fasta -bed ${goldPrefix}.bed -out ${outputPrefix}_bed >${outputPrefix}_bed.out -nb-cores 1 2> /dev/null

################################################################################
# we check the results 
################################################################################

# Checking the .othervariants.vcf :
sh compare_vcf.sh ${outputPrefix}_bed.othervariants.vcf ${goldPrefix}_bed.othervariants.vcf 1> /dev/null 2>&1
var=$?

if [ $var -eq 0 ]
then
echo "full-test find with bed option vcf         : PASS"
else
echo "full-test find with bed option vcf         : FAILED"
RETVAL=1
fi

# Checking the .breakpoints :
#diff --ignore-matching-lines=">" $outputPrefix.breakpoints $goldPrefix.breakpoints 1> /dev/null 2>&1
#var=$?

tmp1=${outputPrefix}_bed.breakpoints.tmp
tmp2=$testDir/tmp2

grep -v "^>" ${outputPrefix}_bed.breakpoints > $tmp1
grep -v "^>" ${goldPrefix}_bed.breakpoints > $tmp2


diff $tmp1 $tmp2 1> /dev/null 2>&1
var=$?


if [ $var -eq 0 ]
then
echo "full-test find breakpoints with bed otpion : PASS"
else
echo "full-test find breakpoints with bed otpion : FAILED"
RETVAL=1
fi

################################################################################
# we launch the fill module
################################################################################
${bindir}/MindTheGap fill -graph $outputPrefix.h5 -bkpt $outputPrefix.breakpoints -out $outputPrefix -nb-cores 1 >>$outputPrefix.out 2> /dev/null

################################################################################
# we check the results 
################################################################################
tmp1=$outputPrefix.insertions.fasta.tmp
tmp2=$testDir/tmp2

grep -v "^>" $outputPrefix.insertions.fasta > $tmp1
grep -v "^>" $goldPrefix.insertions.fasta > $tmp2


diff $tmp1 $tmp2 1> /dev/null 2>&1
var=$?

if [ $var -eq 0 ]
then
echo "full-test fill fasta       : PASS"
else
echo "full-test fill fasta       : FAILED"
RETVAL=1
fi

# Checking the .insertions.vcf :
sh compare_vcf.sh $outputPrefix.insertions.vcf $goldPrefix.insertions.vcf 1> /dev/null 2>&1
var=$?

if [ $var -eq 0 ]
then
echo "full-test fill vcf         : PASS"
else
echo "full-test fill vcf         : FAILED"
RETVAL=1
fi


################################################################################
# we launch the fill module in contig mode
################################################################################
${bindir}/MindTheGap fill -in ../data/contig-reads.fasta.gz -contig ../data/contigs.fasta -abundance-min 3 -out $outputPrefix2 -nb-cores 1 >>$outputPrefix2.out 2> /dev/null

################################################################################
# we check the results
################################################################################
tmp1=$outputPrefix2.insertions.fasta.tmp
tmp2=$testDir/tmp2

grep -v "^>" $outputPrefix2.insertions.fasta > $tmp1
grep -v "^>" $goldPrefix2.insertions.fasta > $tmp2


diff $tmp1 $tmp2 1> /dev/null 2>&1
var=$?

if [ $var -eq 0 ]
then
echo "contig-test fill fasta       : PASS"
else
echo "contig-test fill fasta       : FAILED"
RETVAL=1
fi

# Checking the .gfa :
diff $outputPrefix2.gfa $goldPrefix2.gfa 1> /dev/null 2>&1
var=$?

if [ $var -eq 0 ]
then
echo "contig-test fill gfa         : PASS"
else
echo "contig-test fill gfa         : FAILED"
RETVAL=1
fi

################################################################################
# we launch the graph simplifications script   --> no longer in this repo (see MinYS)
################################################################################
#../pipeline/genome_graph/graph_simplification.py ../pipeline/genome_graph/data/simple4.gfa $outputPrefix.simplified.gfa 1> /dev/null 2>&1

################################################################################
# we check the results 
################################################################################

#nblink=$(grep "^[^L]" $outputPrefix.simplified.gfa  | wc -l)
#nbsegment=$(grep "^[^S]" $outputPrefix.simplified.gfa  | wc -l)

#if [ $nblink -eq 4 ] && [ $nbsegment -eq 4 ]
#then
#echo "graph simplification         : PASS"
#else
#echo "graph simplification         : FAILED"
#RETVAL=1
#fi



################################################################################
# clean up
################################################################################
rm -rf  $testDir

# for Jenkins CI platform, we need an exit code: PASS (0) vs. FAILED (1)
exit $RETVAL
