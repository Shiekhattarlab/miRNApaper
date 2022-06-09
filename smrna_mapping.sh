

FOLDER='/path/to/folder'
cd $FOLDER

mkdir fastq_files bam_files bigwig_files


# copy and rename files to $FOLDER/fastq_files

cat 2020-12-08-Nina-smRNA-01/Files/2020-12-08-Nina-smRNA-01_S* > $FOLDER/fastq_files/filename_r1.fastq.gz
cat 2020-12-08-Nina-smRNA-01/Files/2020-12-08-Nina-smRNA-02_S* > $FOLDER/fastq_files/filename_r2.fastq.gz


cd $FOLDER/fastq_files
for f in *.fastq.gz;
do
cutadapt -m 17 -O 1 --match-read-wildcards -u 3 -a AAAAAAAAAA -o "${f%.fastq.gz}"_cut.fastq.gz \
"$f" > "${f%.fastq.gz}".metrics
done;

# fastqc quality control
# before mapping
cd $FOLDER/fastq_files
mkdir fastqc_before_mapping

for f in *cut.fastq.gz;
do
fastqc "$f" -o $FOLDER/fastq_files/fastqc_before_mapping > \
$FOLDER/fastq_files/fastqc_before_mapping/"${f%.fastq.gz}".control_fastqc
done;

cd fastqc_before_mapping
multiqc . -n multiqc_report_fastqc_before_mapping.html


# store mapped files in respective folder
cd $FOLDER/bam_files
mkdir rep_map hg19_map

# map against repetitive regions
cd $FOLDER/fastq_files
for f in *_cut.fastq.gz;
do
STAR --runMode alignReads  --runThreadN 10  --genomeDir /path/to/references/hg19/STAR/hg19_repeat \
--genomeLoad LoadAndRemove  --readFilesIn "$f" \
--outSAMunmapped Within  --outFilterMultimapNmax 30  --outFilterMultimapScoreRange 1  \
--outFileNamePrefix $FOLDER/bam_files/"${f%.fastq.gz}"_repeat.bam  \
--outSAMattributes All  --readFilesCommand zcat  --outStd BAM_Unsorted  --outSAMtype BAM Unsorted  --outFilterType BySJout
--outReadsUnmapped Fastx  \
--outFilterScoreMin 10  --outSAMattrRGline ID:foo  --alignEndsType EndToEnd \
> $FOLDER/bam_files/"${f%.fastq.gz}"_repeat.bam
done;

# map unmapped reads against hg19 genome
cd $FOLDER/bam_files
for f in *_repeat.bam;
do
STAR  --runMode alignReads  --runThreadN 10  --genomeDir /path/to/references/hg19/STAR/hg19_genome \
--genomeLoad LoadAndRemove  --readFilesIn "$f"Unmapped.out.mate1 \
--outSAMunmapped Within  --outFilterMultimapNmax 1  --outFilterMultimapScoreRange 1  \
--outFileNamePrefix "${f%_repeat.bam}"_rmRep.bam  \
--outSAMattributes All  --outStd BAM_Unsorted  --outSAMtype BAM Unsorted  --outFilterType BySJout  --outReadsUnmapped Fastx  \
--outFilterScoreMin 10  --outSAMattrRGline ID:foo  --alignEndsType EndToEnd \
> "${f%_repeat.bam}"_rmRep.bam
done;


# sort and create bigwig files
cd $FOLDER/bam_files

for f in *_rmRep.bam;
do
samtools sort -@ 10 "$f" -o "${f%_rmRep.bam}"_rmRep_sorted.bam
samtools index "${f%_rmRep.bam}"_rmRep_sorted.bam "${f%_rmRep.bam}"_rmRep_sorted.bam.bai

bamCoverage --bam "${f%_rmRep.bam}"_rmRep_sorted.bam \
-o $FOLDER/bigwig_files/"${f%_rmRep.bam}"_bs1_cpm_reverse.bw -p 5 --normalizeUsing CPM --filterRNAstrand forward -bs 1 --sca
leFactor -1
bamCoverage --bam "${f%_rmRep.bam}"_rmRep_sorted.bam \
-o $FOLDER/bigwig_files/"${f%_rmRep.bam}"_bs1_cpm_forward.bw -p 5 --normalizeUsing CPM  --filterRNAstrand reverse -bs 1
done;



# fastqc after mapping
cd $FOLDER/bam_files
mkdir fastqc_after_mapping

for f in *rmRep_sorted.bam;
do
fastqc "$f" -o $FOLDER/bam_files/fastqc_after_mapping > \
$FOLDER/bam_files/fastqc_after_mapping/"${f%_rmRep_sorted.bam}".after_fastqc
done;

cd fastqc_after_mapping
multiqc . -n multiqc_report_filename_fastqc_after_mapping.html
















# map unmapped output against repetitive genome
cd $FOLDER/bam_files/dme_map
for f in *_dme.bam;
do
STAR --runMode alignReads  --runThreadN 50  --genomeDir /media/labdisk/data/references/hg19/STAR/hg19_repeat \
--genomeLoad LoadAndRemove  --readFilesIn "$f"Unmapped.out.mate1 \
--outSAMunmapped Within  --outFilterMultimapNmax 30  --outFilterMultimapScoreRange 1  \
--outFileNamePrefix $FOLDER/bam_files/rep_map/"${f%.bam}"_hg19_repeat.bam  \
--outSAMattributes All  --outStd BAM_Unsorted  --outSAMtype BAM Unsorted  --outFilterType BySJout  --outReadsUnmapped Fastx  \
--outFilterScoreMin 10  --outSAMattrRGline ID:foo  --alignEndsType EndToEnd \
> $FOLDER/bam_files/rep_map/"${f%.bam}"_hg19_repeat.bam
done;

#map unmapped output against hg19 genome
cd $FOLDER/bam_files/rep_map
for f in *_repeat.bam;
do
STAR  --runMode alignReads  --runThreadN 50  --genomeDir /media/labdisk/data/references/hg19/STAR/hg19_chipseq \
--genomeLoad LoadAndRemove  --readFilesIn "$f"Unmapped.out.mate1 \
--outSAMunmapped Within  --outFilterMultimapNmax 1  --outFilterMultimapScoreRange 1  \
--outFileNamePrefix $FOLDER/bam_files/hg19_map/"${f%_repeat.bam}"_rmRep.bam  \
--outSAMattributes All  --outStd BAM_Unsorted  --outSAMtype BAM Unsorted  --outFilterType BySJout  --outReadsUnmapped Fastx  \
--outFilterScoreMin 10  --outSAMattrRGline ID:foo  --alignEndsType EndToEnd \
> $FOLDER/bam_files/hg19_map/"${f%_repeat.bam}"_rmRep.bam
done;


cd $FOLDER/bam_files/hg19_map/

for f in *_rmRep.bam;
do
samtools sort -@ 40 "$f" -o "${f%_rmRep.bam}"_rmRep_sorted.bam
samtools index "${f%_rmRep.bam}"_rmRep_sorted.bam "${f%_rmRep.bam}"_rmRep_sorted.bam.bai

bamCoverage --bam "${f%_rmRep.bam}"_rmRep_sorted.bam \
-o $FOLDER/bigwig_files/"${f%_rmRep.bam}"_bs1_cpm_reverse.bw -p 50 --normalizeUsing CPM --filterRNAstrand forward -bs 1 --scaleFactor -1
bamCoverage --bam "${f%_rmRep.bam}"_rmRep_sorted.bam \
-o $FOLDER/bigwig_files/"${f%_rmRep.bam}"_bs1_cpm_forward.bw -p 50 --normalizeUsing CPM  --filterRNAstrand reverse -bs 1
done;


# fastqc after mapping
cd $FOLDER/bam_files/hg19_map
mkdir fastqc_after_mapping

for f in *rmRep_sorted.bam;
do
fastqc "$f" -o $FOLDER/bam_files/hg19_map/fastqc_after_mapping > \
$FOLDER/bam_files/hg19_map/fastqc_after_mapping/"${f%_rmRep_sorted.bam}".after_fastqc
done;

cd fastqc_after_mapping
multiqc . -n  multiqc_report_20201210_hela_nuclear_vs_cytoplasm_fastqc_after_mapping.html.html 

##########
# end script 20201211_hela_nuclear_vs_cytoplasm_mapping.sh


#####
#scp ndk25@scccepigenetics.cgcent.miami.edu:$FOLDER/fastq_files/fastqc_before_mapping/multiqc_report_*html \
#/Users/ndk25/Desktop/Nina/data/smrna/20201211_smrna_nuclear_vs_cytoplasm

#scp ndk25@scccepigenetics.cgcent.miami.edu:$FOLDER/bam_files/hg19_map/fastqc_after_mapping/multiqc_report_*html \
#/Users/ndk25/Desktop/Nina/data/smrna/20201211_smrna_nuclear_vs_cytoplasm


# copy bigwig files to amazon using aws conda environment on server
conda activate aws # configuration of environment from Helena's mail 20200811

# need to create respective folder on amazon server, first 

cd $FOLDER/bigwig_files

for f in *.bw;
do
aws s3 cp "$f" s3://ninaki/20201211_smrna_nuclear_vs_cytoplasm/ &
done;




#rename 'rep_map2020' 2020 *dme_hg19_repeat.bam

#############################################################################################

                            # mirdeep on nuclear vs cytoplasm #

#############################################################################################

docker attach miRNA_counts

FOLDER='/ninas_files/20201211_smrna_nuclear_vs_cytoplasm'

cd $FOLDER


# vim loop_mirdeep_nuclear_vs_cytoplasm.sh
#!/bin/bash

# run in docker container (docker attach miRNA_counts)

# nohup bash loop_mirdeep_nuclear_vs_cytoplasm.sh > nohup_loop_mirdeep_nuclear_vs_cytoplasm.out &
# echo $! > save_pid_mirdeep.txt

FOLDER='/ninas_files/20201211_smrna_nuclear_vs_cytoplasm'

mkdir $FOLDER/mapper_files
mkdir $FOLDER/mirdeep2_files


# use unmapped output after dme_map to avoid bias from Drosophila RNA

cd $FOLDER/bam_files/dme_map

for f in *Unmapped.out.mate1;
do
cd $FOLDER/mapper_files
mkdir "${f%.bamUnmapped.out.mate1}"
cd "${f%.bamUnmapped.out.mate1}"
nohup mapper.pl $FOLDER/bam_files/dme_map/"$f" \
-e -h -i -j -l 18 -m -p /ninas_files/ref/hg19_genome.fa \
-s "${f%.bamUnmapped.out.mate1}"_reads_collapsed.fa \
-t "${f%.bamUnmapped.out.mate1}"_reads_vs_ref.arf \
-v -o 4

cd $FOLDER/mirdeep2_files
mkdir "${f%.bamUnmapped.out.mate1}"
cd "${f%.bamUnmapped.out.mate1}"
nohup miRDeep2.pl $FOLDER/mapper_files/"${f%.bamUnmapped.out.mate1}"/"${f%.bamUnmapped.out.mate1}"_reads_collapsed.fa \
/ninas_files/ref/hg19_genome.fa \
$FOLDER/mapper_files/"${f%.bamUnmapped.out.mate1}"/"${f%.bamUnmapped.out.mate1}"_reads_vs_ref.arf \
/ninas_files/ref/mature_ref.fa /ninas_files/ref/mature_other.fa /ninas_files/ref/hairpin_ref.fa -t hsa \
2>"${f%.bamUnmapped.out.mate1}"_report.log

cd $FOLDER/bam_files/dme_map
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





