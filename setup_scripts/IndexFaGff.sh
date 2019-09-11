#!/usr/bin/env bash
set -euxo pipefail

Output_path=$1
santools_path=$2
fasta_name=$3
tabix_path=$4
gff3_name=$5

cd $Output_path

echo "Indexing fasta..."
$santools_path faidx $fasta_name

echo "Indexing gff3 file..."
sort -k1,1 -k4,4n $gff3_name > ${gff3_name}.sorted

rm -f $gff3_name

mv ${gff3_name}.sorted ${gff3_name}

bgzip < ${gff3_name} > ${gff3_name}.gz

$tabix_path -C -p gff ${gff3_name}.gz
