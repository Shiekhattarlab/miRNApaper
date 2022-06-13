# bin/bash!

# generate STAR reference

cd /path/to/folder
mkdir STAR_primir

rsem-prepare-reference \
-gtf primir_final_annotation_v87_v99_StringTie_hg19_manualClean.gtf \
-star \
-p 30 \
hg19_genome.fa \
STAR_primir/primir_rsem

# trim adaptors using trimmomatic
trimmomatic SE -threads 40 -phred33 sample_r1.fastq.gz  sample_r1_trim.fastq.gz \
ILLUMINACLIP:/share/apps/trimmomatic/0.32/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# map trimmed files to premirnas
cd /path/to/fastq_files
mkdir /path/to/bam_files

for f in *trim.fastq.gz;
do
STAR  --runThreadN 20  --genomeDir /path/to/STAR_primir \
--readFilesIn "$f" \
--outFileNamePrefix /path/to/bam_files/"${
f%.fastq.gz}"_premir.bam \
--readFilesCommand gunzip -c  --quantMode TranscriptomeSAM GeneCounts \
--outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate
done;


cd /path/to/bam_files
mkdir rsem_final

for f in *.toTranscriptome.out.bam;
do
rsem-calculate-expression -p 20 \
--bam "$f" \
--forward-prob 0 \
/path/to/STAR_primir/primir_rsem \
rsem_final/"${f%.bamAligned.toTranscriptome.out.bam}"
done;