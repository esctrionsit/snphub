#!/usr/bin/env bash

Output_path=$1
gff3_file=$2

cd $Output_path

gawk -F"[\t=;]" -vOFS="\t" '$3=="gene"{print $1,$4,$5,$10}' ${gff3_file} > geneinfo.txt
