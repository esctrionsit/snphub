#!/usr/bin/env bash
set -euxo pipefail

Output_path=$1
vcf_folder=$2

cd $vcf_folder

echo "Transfering VCF files into BCF..."

bcftools concat '*.vcf.gz' -Ob -o ${Output_path}/unannotated.bcf.gz
