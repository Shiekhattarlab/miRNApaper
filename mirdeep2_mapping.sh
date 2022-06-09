#!/bin/bash

FOLDER='/path/to/folder'

mkdir $FOLDER/mapper_files
mkdir $FOLDER/mirdeep2_files

cd $FOLDER/fastq_files
mkdir unzip_files
for f in *cut.fastq.gz;
do
cp "$f" unzip_files
cd unzip_files
gunzip "$f"
cd ..
done;

cd $FOLDER/fastq_files/unzip_files
for f in *.fastq;
do
cd $FOLDER/mapper_files
mkdir "${f%.fastq}"
cd "${f%.fastq}"
nohup mapper.pl $FOLDER/fastq_files/unzip_files/"$f" \
-e -h -i -j -l 18 -m -p /ninas_files/ref/hg19_genome.fa \
-s "${f%.fastq}"_reads_collapsed.fa \
-t "${f%.fastq}"_reads_vs_ref.arf \
-v -o 4

cd $FOLDER/mirdeep2_files
mkdir "${f%.fastq}"
cd "${f%.fastq}"
nohup miRDeep2.pl $FOLDER/mapper_files/"${f%.fastq}"/"${f%.fastq}"_reads_collapsed.fa \
/ninas_files/ref/hg19_genome.fa \
$FOLDER/mapper_files/"${f%.fastq}"/"${f%.fastq}"_reads_vs_ref.arf \
/ninas_files/ref/mature_ref.fa /ninas_files/ref/mature_other.fa /ninas_files/ref/hairpin_ref.fa -t hsa \
2>"${f%.fastq}"_report.log

cd $FOLDER/fastq_files/unzip_files
done;

#rename csv and pdf files
cd $FOLDER/mirdeep2_files
for d in *;
do
cd "$d"
mv miRNAs_expressed_all_samples_*.csv miRNAs_expressed_all_samples_"$d".csv
mv pdfs_* pdfs_"$d"
cd ..
done;

