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
# we launch the find module
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
# clean up
################################################################################
rm -rf  $testDir

# for Jenkins CI platform, we need an exit code: PASS (0) vs. FAILED (1)
exit $RETVAL
