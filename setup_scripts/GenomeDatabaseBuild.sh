Output_path=$1
gff3_name=$2
fasta_name=$3

cd $Output_path

Output_real_path=$(pwd)

cd -

mkdir Downloads

cd Downloads

wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip

unzip snpEff_latest_core.zip

cd snpEff

rm -f snpEff.config

cp ../../setup_scripts/snpEff.config ./

mkdir data && cd data

mkdir SnphubBuilding

mkdir genomes

cp ${Output_real_path}/${gff3_name} SnphubBuilding/

mv SnphubBuilding/$gff3_name SnphubBuilding/genes.gff

cp ${Output_real_path}/${fasta_name} genomes/

mv genomes/${fasta_name} genomes/SnphubBuilding.fa

cd ../

java -jar snpEff.jar build -gff3 -v SnphubBuilding
