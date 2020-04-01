#!/usr/bin/env bash
set -euxo pipefail

Output_path=$1
vcf_folder=$2
gff_name=$3
path_bcftools=$4

tool_path=$(pwd)/Downloads/snpEff/snpEff.jar

cd $vcf_folder

vcf_real_path=$(pwd)

cd -

cd $Output_path

echo "Indexing VCF files..."

$path_bcftools index *.vcf.gz

echo "Converting VCF files into BCF..."

$path_bcftools concat ${vcf_real_path}/*.vcf.gz -Ob -o unannotated.bcf.gz

echo "Indexing BCF file..."

$path_bcftools index unannotated.bcf.gz

echo "Annotating BCF file..."

$path_bcftools view unannotated.bcf.gz --threads 4 \
    | java -jar $tool_path -t SnphubBuilding - \
    | $path_bcftools view --threads 4 -o output.ann.bcf.gz -Ob

echo "Indexing annotated BCF file..."

$path_bcftools index output.ann.bcf.gz

