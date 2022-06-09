# Strategy for final miRNA analysis, tracks and deseq2

# Subsample all files to 30 million reads
# Merge shGFP replicates (and shINTS11) and also subsample to 30 million reads



##########################################################################################

								# SUBSAMPLES AND MERGE #
								
##########################################################################################

FOLDER='/path/to/folder'
cd $FOLDER
mkdir original_fastq subsampled_fastq


# link to original fastqs
ln -s \
/path/to/original/fastq_files/*r{1,2}.fastq.gz original_fastq


# subsampling to 30 million
for f in *.fastq.gz;
do
cp "$f" unzip_files
cd unzip_files
gunzip "$f"
cd .. 
done;

# from http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/subsampling_reads.pdf
cd $FOLDER/original_files/unzip_files
for f in *.fastq;
do
cat "$f" | awk '{ printf("%s",$0); n++; if(n%4==0) {
printf("\n");} else { printf("\t");} }' |
awk -v k=30000000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | 
awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4}' > $FOLDER/subsampled_fastq/unzip_files/"${f%.fastq}"_30mio.fastq
done;

## after subsampling double check number of reads:
cd $FOLDER/subsampled_fastq/unzip_file
for f in *mio.fastq;
do
echo "$f"
echo $(cat "$f"|wc -l)/4|bc
done;

# zip files
for f in *30mio.fastq;
do
echo "$f"
gzip "$f"  
done;


# proceed with adapter trimming, mapping against repetitive regions, genome mapping