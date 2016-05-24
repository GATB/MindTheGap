#! /bin/bash

## Script to compare vcf, in a smarter way than a simple diff : comparing only the fields chrom position ref alt (must be identical)

## 2 arguments = 2 vcf files

vcf1=$1
vcf2=$2

tmp1=$vcf1.temp
tmp2=$vcf2.temp

grep -v "^#" $vcf1 | cut -f1,2,4,5 | sort > $tmp1 
grep -v "^#" $vcf2 | cut -f1,2,4,5 | sort > $tmp2 

diff $tmp1 $tmp2
