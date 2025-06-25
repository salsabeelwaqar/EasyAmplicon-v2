
# Conda software installation directory, check with `conda env list`, e.g., /anaconda3  
## Set up Miniconda environment
soft=~/miniconda3  

# Set working directory (wd), e.g., /d/Amplicon  
wd="/mnt/d/Amplicon"

# Set database directory (db) to /mnt/d/EasyMicrobiome2
db="/mnt/d/EasyMicrobiome2"

# Change to the working directory
mkdir -p $wd && cd $wd

# Optionally, check the directories
echo "Working directory: $wd"
echo "Database directory: $db"

## 1. Initialize Files

# Create required directories: seq (for sequencing data), result (for metadata and outputs), temp (temporary storage)
mkdir -p seq result temp

# ├── pipeline.sh
# ├── result
# │   └── metadata.txt
# ├── seq
# │   ├── HS_1.fq.gz
# │   ├── ...
# │   └── ZYmo_2.fq.gz
# └── temp

### 1.1. Metadata Preparation

# Save sample metadata in "result/metadata.txt"
# Use "csvtk" to count the number of rows (samples, excluding the header) and columns in a table. Use the -t option to set the column delimiter as a tab (default is comma);
csvtk -t stat result/metadata_raw.txt
# The metadata should have at least 3 columns. The first column should be the sample ID (SampleID), and the last column should be the description (Description)
# You can use the cat command to view a file, and the -A option to display symbols. The '|' symbol is used as a pipe to connect commands. The head command is used to display the file header, and the -n3 option limits the output to the first 3 lines.
cat -A result/metadata_raw.txt | head -n3
# Windows users end with ^M, run the sed command to remove, and then check with cat -A
sed 's/\r//' result/metadata_raw.txt > result/metadata.txt
cat -A result/metadata.txt | head -n3

### 1.2. Sequencing Data Preparation
The sequencing results returned by the company are usually single FastQ compressed files in .gz format for one sample
# The file name and the sample name must be corresponding: manually modify it when it is inconsistent.
# If the sequencing data is a compressed file of .gz, sometimes you need to use gunzip to decompress and use, vsearch can usually read the compressed file directly
# gunzip seq/*.gz
# zless view compressed files by page, space page turn, q exit; head looks at the first 10 rows by default, and -n specifies the rows
# Verify downloaded sequencing files
The sample downloads a single file and renames it
# mkdir -p seq
ls -sh seq/
## View sequence file content
zless seq/HS_1.fastq | head -n4
# If needed, specify character range per line
zless seq/HS_1.fastq | head | cut -c 1-60

#Check sequencing file statistics using seqkit
seqkit stat seq/HS_1.fastq

## Summarize sequencing stats for all samples
seqkit stat seq/*.fastq > result/seqkit.txt
head result/seqkit.txt

### 1.3.Pipeline & databases

# The database must be decompressed the first time it is used, and you can skip this section later
# usearches available 16S/18S/ITS databases: RDP, SILVA and UNITE, local file location ${db}/usearch/
# Usearch database database download page: http://www.drive5.com/usearch/manual/sintax_downloads.html
# Decompress 16S RDP database, gunzip decompress, seqkit stat statistics
# Keep the original compressed file
gunzip -c ${db}/usearch/rdp_16s_v18.fa.gz > ${db}/usearch/rdp_16s_v18.fa
seqkit stat ${db}/usearch/rdp_16s_v18.fa # 21 thousand sequences
## for silva databse 
gunzip -c ${db}/usearch/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz> ${db}/usearch/SILVA_138.2_SSURef_NR99_tax_silva.fa
seqkit stat ${db}/usearch/SILVA_138.2_SSURef_NR99_tax_silva.fa


# To decompress the ITS UNITE database, you need to download it from the official website or the network disk db/amplicon/usearch
# gunzip -c ${db}/usearch/utax_reference_dataset_all_29.11.2022.fasta.gz >${db}/usearch/unite.fa
seqkit stat ${db}/usearch/unite.fa # 326 thousand
# The Greengene database is used for feature annotations: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
# The default decompression will delete the original file, -c specifies the output to the screen, > Write a new file (can be renamed):
gunzip -c ${db}/gg/97_otus.fasta.gz > ${db}/gg/97_otus.fa
seqkit stat ${db}/gg/97_otus.fa
 
# In order to perform quality filtering,first assess sequence length and quality with nanoplot
mkdir -p nanoplot_reports  
##If you're using conda, run:
conda install -c bioconda nanoplot
## If already installed then activate nanoplot (if you already have installed nanoplot_env). 
conda activate nanoplot_env 
# Set directory paths (adjusted for your environment)
INPUT_DIR="/mnt/d/Amplicon/temp"
OUTPUT_DIR="/mnt/d/Amplicon/nanoplot_reports"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# List of merged FASTQ files
SAMPLES=("HS_1" "HS_2" "HS_3" "NW_1" "NW_2" "NW_3" "YW_1" "YW_2" "YW_3" "Zymo_1" "Zymo_2" "Zymo_3")

# Run NanoPlot
for SAMPLE in "${SAMPLES[@]}"; do
echo "Processing $SAMPLE..."
NanoPlot --fastq "$INPUT_DIR/${SAMPLE}.merged.fq" \
--outdir "$OUTPUT_DIR/${SAMPLE}_report/" \
--plots hex dot --loglength --N50
done


## 2. Reads Rename

# Set working directory for clarity
cd /mnt/d/Amplicon

# Run relabeling with a more descriptive output filename
time for i in $(tail -n +2 result/metadata.txt | cut -f1); do
usearch -fastx_relabel "seq/${i}.fastq" -fastqout "temp/${i}.fq" -prefix "${i}."
done

## 2. Trim primers and perform quality control (for Long read data) 
  
##Once you have the FASTQ file (basecalled reads), the next steps depend on the objectives of your analysis. Here’s a general workflow for processing full-length 16S rRNA amplicon sequencing data, such as from Nanopore sequencing:
##1. Quality Control: Remove low-quality reads to ensure high-quality downstream analysis.Perform these method under window subsytem for linux (WSL)
  
# Step 2 Next, quality filtering is conducted by assessing sequence length and quality through 'NANOPLOT'/fastqc result.

## Apply cutadapt for each samples 
conda deactivate

#!/bin/bash
# Define paths
IN_DIR="/mnt/d/Amplicon2/temp"  # <<< Using merged files now
OUT_DIR="/mnt/d/Amplicon2/temp/filtered"
QC_DIR="${OUT_DIR}/qc"

# Create output directories if they don't exist
mkdir -p "$OUT_DIR" "$QC_DIR"

# Loop through each merged FASTQ file
for file in "$IN_DIR"/*.merged.fq; do
# Get sample name from filename
sample=$(basename "$file" .merged.fq)

echo "Processing $sample..."

# Step 1: Cutadapt
cutadapt -g "AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA" \
--error-rate=0.1 \
--action=trim \
--rc \
-m 1000 \
-M 1800 \
-j 20 \
--discard-untrimmed \
-o "${OUT_DIR}/${sample}.filtered.fastq" \
"$file"

# Step 2: FastQC
fastqc "${OUT_DIR}/${sample}.filtered.fastq" -o "$QC_DIR"
done

# Merge all filtered FASTQ files
cat "${OUT_DIR}"/*.filtered.fastq > "${OUT_DIR}/all.filtered.fq"

echo "All samples processed and merged into ${OUT_DIR}/all.filtered.fq"

### 2.3 Copy and rename all.filtered.fq to all.fq in the temp folder 

cp "${OUT_DIR}/all.filtered.fq" "${IN_DIR}/all.fq"
echo "File copied and renamed as ${IN_DIR}/all.fq"

# Count reads in all.fq
echo "Counting reads in all.fq:"
grep -c '^@' "${IN_DIR}/all.fq"

# Optional: Count reads in each filtered file
echo "Counting reads in each filtered.fastq file:"
for file in "${OUT_DIR}"/*.filtered.fastq; do
sample=$(basename "$file" .filtered.fastq)
count=$(grep -c '^@' "$file")
echo "$sample: $count reads"
done

#View file size: 223MB. There are slight differences in results due to different software versions.
ls -lh temp/all.fq
# Check if the sequence name before the '.' is the sample name. Sample names must not contain dots ('.')
# One notable characteristic of sample names with dots (.) is that the generated feature table becomes significantly large, with numerous columns, leading to insufficient memory during subsequent analyses。
# After obtaining the feature table from the subsequent analysis, it's important to take a look to check for any issues. If you encounter memory-related problems, you'll need to go back and troubleshoot。
head -n 6 temp/all.fq|cut -c1-60

## Method 2: Optional  :Apply Porechop directly on your all.fq file located inside the temp folder. If not install just write sudo apt install porchop. 
porechop -i /mnt/d/Amplicon2/temp/all.fq -o /mnt/d/Amplicon2/temp/all_trimmed.fq

Sttep 5 # Be sure to understand the experimental design and primer lengths. If primers have been removed, you can enter 0. Processing 270,000 sequences will take around 14 seconds. As the FASTQ file contains Phred scores up to 50, so fix "--fastq_qmax" to 50. And "--fastq_maxee_rate 0.01" removes reads with too many errors.--fastaout temp/filtered.fa saves only high-quality sequences in FASTA format, which is required for OTU/ASV clustering.
time vsearch --fastx_filter temp/all.fq \
--fastq_stripleft 0 --fastq_stripright 0 \
--fastq_maxee_rate 0.01 \
--fastq_qmax 50 \
--fastaout temp/filtered.fa


#### 2. Remove redundancy and select OTU/ASV 
###  Sequence de-redundancy 
# And add a miniuniqusize of at least 2 or 1/1M to remove low-abundance noise and increase calculation speed. Or 1/1M means that alternatively, a feature should have a frequency of at least 1 per million (1/1M) to be retained. This is useful when dealing with large datasets where low-abundance sequences may be treated as noise.
# -sizeout outputs abundances, and --relabel must be used with sequence prefixes for better organization. Takes 1 second
vsearch --derep_fulllength temp/filtered.fa \
--minuniquesize 1 --sizeout --relabel Uni_ \
--output temp/uniques.fa 
# Uni_1;size=590" - The name of the sequence after dereplication. This sequence appeared 590 times in the sequencing data across all samples
# is the sequence that appears the most。
head -n 2 temp/uniques.fa
##Count the number of sequences:(712392)
grep -c "^>" temp/uniques.fa

###Option 2:with miniuniquesize number 2-10 (try and see the output , after that choose which one is best) 
with minuniquesize= 5 obtained only (285)## very little and with 10 even more less)

adjust to minuniquesize= 2
vsearch --derep_fulllength temp/filtered.fa \
--minuniquesize 2 --sizeout --relabel Uni_ \
--output temp/uniques2.fa 
#Check output summary
ls -lsh temp/uniques2.fa
##Count the number of sequences:(1951)
grep -c "^>" temp/uniques2.fa

# 2.2 : De novo chimera filtering
vsearch --uchime_denovo temp/uniques2.fa \
--nonchimeras temp/uniques_nochim2.fa

##Count the number of sequences:
grep -c "^>" temp/uniques_nochim.fa ##(677701)
grep -c "^>" temp/uniques_nochim2.fa  ##(1471)

## 3. The data is too large to use USEARCH for clustering or denoising. Replace it with vsearch

# When restricted to the limitations of the free version of USEARCH, you can reduce the non-redundant data size by increasing the "minuniquesize" parameter. If you're dealing with over 10,000 OTUs/ASVs and the downstream analysis is taking too long, ensure that the OTU/ASV data is under 5000. Typically, this won't be restricted, and it's also advantageous for faster downstream analyses
# Using vsearch for clustering to generate OTUs is an option, but it lacks automatic de novo chimera removal. With an input of 2155 sequences, the clustering process resulted in 661 output clusters.
# Rename using "relabel" and cluster at a similarity threshold of 97% without masking qmask
# Record the input size ("sizein") and the output frequency ("sizeout"): theoutput is 1487269859 nt in 1024341 seqs, min 1200, max 1662, avg 1452, Clusters: 410369 Size min 1, max 49857, avg 2.5, Singletons: 363336, 35.5% of seqs, 88.5% of clusters.
###Method:2 OTU clustering at 97% identity
vsearch --cluster_size temp/uniques_nochim2.fa \
--relabel ASV_ \
--id 0.97 \
--qmask none \
--sizein \
--sizeout \
--centroids temp/otus_raw2.fa

#Clean sequence headers (remove ";size=...")
sed 's/;size=[0-9]*//' temp/otus_raw2.fa > temp/otus3.fa
 #Count how many ASVs you got 
grep -c "^>" /mnt/d/Amplicon2/temp/otus3.fa
# Preview the output
head -n 4 temp/otus.fa
grep -c "^>" temp/otus.fa (288382)
grep -c "^>" temp/otus3.fa (338)

##Method:2: Denoise sequences to obtain ASVs with single-base accuracy (Yield into small number of ASVs when using minisize 10-2, while stuck ion reducing minisize to 1  due to data size or limitations")
usearch -unoise3 temp/uniques.fa -minsize 2 \
-zotus temp/zotus.fa

usearch -unoise3 temp/uniques_nochim.fa -minsize 2 \
-zotus temp/zotus44.fa

usearch -unoise3 temp/otus.fa -minsize 10 \
-zotus temp/zotus45.fa

# Rename output sequences from "Zotu" to "ASV" for clarity
sed 's/Zotu/ASV_/g' temp/zotus45.fa > temp/otus45.fa

# View the first two lines to confirm
head -n 2 temp/otus2.fa
grep -c "^>" temp/otus2.fa (814) 

grep -c "^>" temp/otus44.fa (98)

grep -c "^>" temp/otus45.fa (98)
#### 4. Aligment with minimap2 and plosihing with meadaka 

## 4.ALignment with minimap2: 
##For Installation, Minimap2 is optimized for x86-64 CPUs.(if not already installed), Skip if already installed.

curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.28_x64-linux/minimap2
# Set the database directory path (adjust as needed)
db="/mnt/d/EasyMicrobiome2/usearch"
# Move to database directory
cd $db
# Create minimap2 index from the reference FASTA
minimap2 -d reference_16S.mmi silva_16s_v138.fa
##Align OTUs/ASVs to Reference with Minimap2
# Back to working directory
cd /mnt/d/Amplicon2
# Align OTUs to reference and create SAM/BAM/Sorted BAM 
## This step will take long time (approx Real time: 6676.350 sec)
minimap2 -ax map-ont /mnt/d/EasyMicrobiome2/usearch/reference_16S.mmi /mnt/d/Amplicon2/temp/otus3.fa > /mnt/d/Amplicon2/temp/aligned3_output.sam

# Convert and sort the alignment
samtools view -Sb temp/aligned3_output.sam > temp/aligned3_output.bam
samtools sort temp/aligned3_output.bam -o temp/aligned3_output_sorted.bam
samtools index temp/aligned3_output_sorted.bam

# Index FASTA for medaka
samtools faidx temp/otus.fa

# View summary statistics
samtools flagstat temp/aligned_output_sorted.bam
# View BAM file header
samtools view -H temp/aligned_output_sorted.bam

#Check overall statistics
samtools flagstat temp/aligned_output_sorted.bam
## In order to view the header of the bam file , run this code 
samtools view -H temp/aligned_output_sorted.bam

##Option2 " Alternative Alignment using VSEARCH
vsearch --usearch_global temp/filtered.fa \
--db /mnt/d/EasyMicrobiome2/usearch/silva_16s_v138.fa \
--strand both \
--id 0.80 \
--blast6out temp/alignment_results.tsv \
--alnout temp/alignment_results.aln \
--samout temp/alignment_results.sam

## 5. ## Polishing Withe Medaka 
#Install Medaka
conda install -c bioconda medaka
#If already installed 
conda activate medaka
##Prepare Your Data
##Make sure you have your aligned BAM file (temp/aligned_output_sorted.bam) and the original Nanopore reads file (temp/otus.fa).
medaka_consensus -i temp/aligned2_output_sorted.bam -d temp/otus2.fa -o temp/medaka_output2 -t 4

##If, for any reason, the program fails to create the consensus.fasta file, and you have a previously created consensus_probs.hdf file, run this code
medaka_consensus -i temp/medaka_output/consensus_probs.hdf -d temp/otus.fa -o temp/medaka_output/ -t 4
##Then convert the consensus.fasta" into "consensus.fa
cd /mnt/d/Amplicon2/temp/medaka_output
mv "consensus.fasta" "consensus.fa"
##If size number present in consesus.fasta, remove using this code
sed 's/;size=[0-9]*//' consensus.fasta > consensus.fa


####6 Perform reference-based chimera detection.
# Removing all OTU_1 sequences seems reasonable since Usearch does not have built-in chimera removal methods。
mkdir -p result/raw

# Reference-based chimera removal using SILVA "silva_reformatted.fasta'

mv "/mnt/d/EasyMicrobiome2/usearch/silva_reformatted.fasta" "/mnt/d/EasyMicrobiome2/usearch/silva_reformatted.fa"
cp usearch/silva_reformatted.fasta usearch/silva_16ss_v138.fa
db="/mnt/d/EasyMicrobiome2"
cd /mnt/d/Amplicon2
vsearch --uchime_ref temp/medaka_output/consensus.fa \
-db ${db}/usearch/silva_16s_v138.fa \
--nonchimeras result/raw/otus.fa \
--chimeras result/raw/otus.chimeras.fa \
--borderline result/raw/otus.borderline.fa
## If you're using Windows and vsearch results have added Windows line endings (^M), you need to remove them. However, if you're using a Mac, you don't need to execute this command
sed -i 's/\r//g' result/raw/otus.fa

##Taxonomic Assignment

##Taxonomic Classification: Assign taxonomy to OTUs or ASVs.

##step 1: Generate a Feature table 

##Method 1:
vsearch --usearch_global temp/filtered.fa \
--db result/raw/otus.fa \
--id 0.97 \
--threads 4 \
--otutabout result/raw/otutab.txt

# For Windows users, you should remove the newline characters (^M) from vsearch results and correct them to the standard Linux format.
sed -i 's/\r//' result/raw/otutab.txt
head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
# To count rows and columns in a table using "csvtk"
# Here, it's important to carefully check the number of columns and whether it matches the number of your samples. If they're not equal, there might be issues with sample naming. Refer to the explanations above for more details.
csvtk -t stat result/raw/otutab.txt


##Method 2 (Not working this code): The output file generate here name as "otutab.txt2" just for noting the clear difference with the output generated above. 
db="/mnt/d/EasyMicrobiome2"
usearch="${db}/win/usearch.exe"
reads="/mnt/d/Amplicon2/temp/filtered.fa"
asvs="/mnt/d/Amplicon2/result/raw/otus.fa"
output="/mnt/d/Amplicon2/result/raw/otutab.txt2"
$usearch -usearch_global $reads \
-db $asvs \
-strand both \
-id 1.0 \
-maxaccepts 8 \
-maxrejects 64 \
-top_hit_only \
-otutabout $output

##step 2

  ## Generate a Feature table Taxonomy Assignment (Use the following command to assign taxonomy to your OTUs using the SILVA database:)

  ##Methods 1: Use this code with latest version of Silva database "SILVA_138.2_SSURef_NR99_tax_silva.fa" formatted version "silva_reformatted.fa"
vsearch --sintax result/raw/otus.fa \
--db ${db}/usearch/silva_16s_v138.fa \
--sintax_cutoff 0.1 \
--tabbedout result/raw/otus.sintax
head result/raw/otus.sintax | cat -A
sed -i 's/\r//' result/raw/otus.sintax

  ## Run the code to identify the taxonomic classification with cutoff '0.8'. 
vsearch --sintax result/raw/otus.fa \
--db ${db}/usearch/silva_16s_v138.fa \
--sintax_cutoff 0.8 \
--tabbedout result/raw/otus2.sintax 
head result/raw/otus2.sintax | cat -A
sed -i 's/\r//' result/raw/otus2.sintax

##Optional: Use RDB silva database 

vsearch --sintax result/raw/otus.fa \
--db ${db}/usearch/rdp_16s_v18.fa \
--sintax_cutoff 0.1 \
--tabbedout result/raw/otus.sintax 
head result/raw/otus.sintax | cat -A
sed -i 's/\r//' result/raw/otus.sintax


Method 2: Taxonomy Assignment.Use the SINTAX method for taxonomy assignment, which aligns with the text you provided.
  # Assign paths to variables (corrected the quotes and formatting)
  # Run this code with silva_16s_v123.fa database and 
  # Define variables
vsearch="${db}/win/vsearch.exe"  # Correct path to the vsearch executable
taxonomy_db="${db}/usearch/silva_reformatted.fa"  # Correct path to the SILVA database
taxonomy_output="/d/Amplicon2/result/raw/otus.sintax"  # Output file path
 # Run vsearch with the corrected parameters
$vsearch --sintax /d/Amplicon2/result/raw/otus.fa \
--db $taxonomy_db \
--strand plus \
--tabbedout $taxonomy_output \
--sintax_cutoff 0.8
 # Check the output
head /d/Amplicon2/result/raw/otus.sintax | cat -A
 # Fix line endings in the output file
sed -i 's/\r//' /d/Amplicon2/result/raw/otus.sintax


# Method1. The number of original feature table rows
wc -l result/raw2/otutab.txt
#R script selects bacterial archaea (eukaryotes), removes chloroplasts, mitochondria and counts the proportions; Output filtered and sorted OTU tables
#The input is the OTU table result/raw/otutab.txt and the species annotation result/raw/otus.sintax
#Output filtered and sorted feature tables result/otutab.txt sum
#Statistical contamination ratio file result/raw/otutab_nonBac.txt and filter details otus.sintax.discard
#For fungal ITS data, use the otutab_filter_nonFungi.R script instead, only screening fungi
# Rscript ${db}/script/otutab_filter_nonBac.R -h # Displays a description of the parameter
Rscript ${db}/script/otutab_filter_nonBac.R \
--input result/raw2/otutab.txt \
--taxonomy result/raw2/otus.sintax \
--output result/otutab.txt \
--stat result/raw2/otutab_nonBac.stat \
--discard result/raw2/otus.sintax.discard
# The number of filtered feature table rows
wc -l result/otutab.txt
#Filter the sequence corresponding to the feature table
cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
usearch -fastx_getseqs result/raw/otus.fa \
-labels result/otutab.id -fastaout result/otus.fa
#Filter feature tables corresponding to sequence annotations
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' \
result/raw2/otus.sintax result/otutab.id > result/otus.sintax


## In oredr to run and save the result of the 'otutab.txt' obtained through miniunique size=1" without using "usearch-unoise with miniuniquesize=2".
mkdir -p result2

 # Method1. The number of original feature table rows
wc -l result/raw/otutab.txt
mkdir -p result2
  #R script selects bacterial archaea (eukaryotes), removes chloroplasts, mitochondria and counts the proportions; Output filtered and sorted OTU tables
  #The input is the OTU table result/raw/otutab.txt and the species annotation result/raw/otus.sintax
  #Output filtered and sorted feature tables result/otutab.txt sum
  #Statistical contamination ratio file result/raw/otutab_nonBac.txt and filter details otus.sintax.discard
  #For fungal ITS data, use the otutab_filter_nonFungi.R script instead, only screening fungi
  # Rscript ${db}/script/otutab_filter_nonBac.R -h # Displays a description of the parameter
Rscript ${db}/script/otutab_filter_nonBac.R \
--input result/raw/otutab.txt \
--taxonomy result/raw/otus.sintax \
--output result2/otutab.txt \
--stat result/raw/otutab_nonBac.stat \
--discard result/raw/otus.sintax.discard
# The number of filtered feature table rows
wc -l result2/otutab.txt
#Filter the sequence corresponding to the feature table
cut -f 1 result2/otutab.txt | tail -n+2 > result2/otutab.id
usearch -fastx_getseqs result/raw/otus.fa \
-labels result2/otutab.id -fastaout result2/otus.fa
#Filter feature tables corresponding to sequence annotations
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' \
result/raw/otus.sintax result2/otutab.id > result2/otus.sintax

# Optional statistical method: Summary OTUs table
usearch -otutab_stats result2/otutab.txt \
-output result2/otutab.stat
cat result2/otutab.stat
#Note the minimum, quantile, or look at the amount of sample detail in result/raw/otutab_nonBac.stat for resampling

### 5.3 Equalization of sampling

# Normlize by subsample
wd=/d/Amplicon2
db=/d/EasyMicrobiome2
PATH=$PATH:${db}/win
cd ${wd}

#Use the vegan package for equal resampling and enter the reads count format Feature table result/otutab.txt
#You can specify the input file, sample size, and random number, and output the equalizer result/otutab_rare.txt and diversity alpha/vegan .txt
mkdir -p result2/alpha
Rscript ${db}/script/otutab_rare.R --input result2/otutab.txt \
--depth 10000 --seed 1 \
--normalize result2/otutab_rare.txt \
--output result2/alpha/vegan.txt
usearch -otutab_stats result2/otutab_rare.txt \
-output result2/otutab_rare.stat
cat result2/otutab_rare.stat

##Analysis with result folder 

mkdir -p result/alpha
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
--depth 10000 --seed 1 \
--normalize result/otutab_rare.txt \
--output result/alpha/vegan.txt
usearch -otutab_stats result/otutab_rare.txt \
-output result/otutab_rare.stat
cat result/otutab_rare.stat








##Taxonomic classification with Kraken2 (Under process)

#if you're using a package manager like conda, you can install Kraken2 directly:
conda install -c bioconda kraken2
# Option2: Download Kraken2 from GitHub
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
make

#Download a pre-built Kraken2 database that includes 16S rRNA sequences (potentially the full SILVA):
kraken2-build --download-db silva_138 --db /mnt/d/EasyMicrobiome2/kraken2_db
#Run Kraken2 Classification:
kraken2 --db /mnt/d/EasyMicrobiome2/kraken2_db \
--output result/taxonomy/kraken2/kraken2_output.txt \
--report result/taxonomy/kraken2/kraken2_report.txt \
D:/Amplicon2/temp/medaka_output/consensus.fasta \
--threads 4
####Remember to adjust the number 8 to match your system's physical core count.
# Option 2: Add the library (SILVA database) to your Kraken2 database
kraken2-build --add-to-library /mnt/d/EasyMicrobiome2/usearch/silva_reformatted.fasta --db /mnt/d/EasyMicrobiome2/kraken2_db

# Step 2: Download the Kraken2 Taxonomy Data:

##This command is generally NOT needed if you are using the --download-db silva_138 option. The pre-built SILVA database already includes the necessary taxonomy information. This step is relevant if you were building a custom database from your own sequences using --add-to-library (which isn't the recommended approach here).	
kraken2-build --download-taxonomy --db D:/EasyMicrobiome2/kraken2_db

# Step 3: Build the Kraken2 database
# If you used --download-db silva_138, you likely do NOT need this step. The downloaded database is usually ready to use. Running --build afterwards might be redundant or could potentially cause issues. Only use --build if you were adding your own sequences with --add-to-library.
kraken2-build --build --db /mnt/d/EasyMicrobiome2/kraken2_db




### 5.3 Equalization of sampling

# Normlize by subsample

#Use the vegan package for equal resampling and enter the reads count format Feature table result/otutab.txt
#You can specify the input file, sample size, and random number, and output the equalizer result/otutab_rare.txt and diversity alpha/vegan .txt

mkdir -p result/alpha-testing
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
--depth 10000 --seed 1 \
--normalize result/otutab_rare.txt \
--output result/alpha-testing/vegan.txt
usearch -otutab_stats result/otutab_rare.txt \
-output result/otutab_rare.stat
cat result/otutab_rare.stat

#While running the above code if all groups are not present in the 'otutab_rare.txt' file its mean that groups being excluded due to insufficient read counts. In order to checking the total reads per sample across all groups with the following command: 
#run the following command:  This will give you the total read count for each sample. If any sample has a count significantly lower than others, it may be the reason it's excluded during the rarefaction process.
awk '{sum=0; for (i=2; i<=NF; i++) sum+=$i; print $1, sum}' result/otutab.txt
#nspect the Rarefied OTU Table:
head -n 10 result/otutab_rare.txt
#If some groups are still missing, it indicates that their total read counts are too low to meet the set depth.
#Adjust the Rarefaction Depth:If the group with the lowest read count (e.g., below 15000) is being excluded, consider lowering the rarefaction depth. For example, you can set the rarefaction depth to a value like 7000 or 10000:
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt --depth 5000 --seed 1 --normalize result/otutab_rare.txt --output result/alpha/vegan.txt
usearch -otutab_stats result/otutab_rare.txt \
-output result/otutab_rare.stat
cat result/otutab_rare.stat







## 6. alpha (α) diversity

### 6.1. calculate alpha (α) diversity

# Use USEARCH to calculate 14 alpha diversity indices (don't use Chao1 if you make a mistake
#details in http://www.drive5.com/usearch/manual/alpha_metrics.html
usearch -alpha_div result/otutab_rare.txt \
-output result/alpha/alpha.txt

### 6.2. calculate rarefaction richness

#Dilution Curve: Take the number of OTUs (Operational Taxonomic Units) from sequences ranging from 1% to 100%, with non-replacement sampling each time.
#Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch -alpha_div_rare result/otutab_rare.txt \
-output result/alpha/alpha_rare.txt \
-method without_replacement
#Preview of Results
head -n2 result/alpha/alpha_rare.txt
#Handling of non-numeric "-" values due to low sample sequencing, see FAQ #8 for details
sed -i "s/-/\t0.0/g" result/alpha/alpha_rare.txt

### 6.3. Filtering for high-abundance bacteria

#Calculate the mean of each feature, and if there are groups, calculate the group means. Group column names should be modified based on the experimental design in metadata.txt.
#The input files are the feature table (result/otutab.txt) and the experimental design metadata (metadata.txt)
#The output is the feature table with group-wise means. A single experiment may have multiple grouping methods
#The "-h" option displays the script's help, providing explanations for the parameters
Rscript ${db}/script/otu_mean.R -h
#The parameter "scale" determines whether to perform scaling, "zoom" normalizes the total sum, "all" outputs mean for all samples, and "type" specifies the calculation type as "mean" or "sum"
Rscript ${db}/script/otu_mean.R --input result/otutab.txt \
--metadata result/metadata.txt \
--group Group --thre 0 \
--scale TRUE --zoom 100 --all TRUE --type mean \
--output result/otutab_mean.txt
# The output includes both the overall mean and the mean values for each group
head -n3 result/otutab_mean.txt

#By filtering with an average abundance threshold of >0.1%, which can be selected as 0.5% or 0.05%, you will obtain the OTU combinations for each group
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=3;i<=NF;i++) a[i]=$i; print "OTU","Group";} \
        else {for(i=3;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
result/otutab_mean.txt > result/alpha/otu_group_exist.txt
head result/alpha/otu_group_exist.txt
cut -f 2 result/alpha/otu_group_exist.txt | sort | uniq -c
# Give it a try: How many OTUs/ASVs does each group have at different abundances?
# Using the following link http://ehbio.com/test/venn/ You can plot and visualize Venn or network diagrams for shared and unique components of each group at different abundances 
# Using the following link http://www.ehbio.com/ImageGP Create Venn, UpSetView, and Sankey diagrams


## 7. Beta(β) diversity

#The results generate multiple files and require a directory (replace otus.fa with otus.nonchimeras.fa )
mkdir -p result/beta/
#Building an evolutionary tree based on OTUs
usearch -cluster_agg result/otus.fa -treeout result/otus.tree
#Generate five distance matrices：bray_curtis, euclidean, jaccard, manhatten, unifrac
usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
-filename_prefix result/beta/
  
  ## While ruuning above code, user may encounter fatal error and alignment issue, like this " Align your sequences using MAFFT (to ensure they are aligned properly before building a tree).
  
  #Use FastTree to build the evolutionary tree from the aligned sequences.
  #Compute distance matrices using QIIME 2 or R packages (e.g., microeco) if needed.
  ##Download and install fastree. FastTree infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences. FastTree can handle alignments with up to a million of sequences in a reasonable amount of time and memory. 
  #Save the file as FastTree.exe in a directory you can easily access, such as D:\EasyAmplicon or a dedicated tools directory like D:\software.Then add FastTree to PATH.Open Open Environment Variables, Under System Variables, find Path and click Edit. Add the directory containing FastTree.exe (e.g., D:\EasyAmplicon) to the list
  ##Save and restart your terminal. Then run the following command to confirm installation (FastTree -help). If the installation is successful, you should see the help message for FastTree

FastTree -nt result/otus.fa > result/otus.tree

## The error Wrong number of characters for ASV_1: expected 1646 but have 1464 instead indicates that there is an inconsistency in sequence lengths in the result/otus.nonchimeras.fa file. FastTree requires all sequences to have the same length, as it expects an aligned sequence file.  
## For alignment user can choose tools like MAFFT or Clustal Omega. 
## while after downloading, installing window 
iqtree -s result/otus.fa -nt AUTO -pre result/otus

  
  
## 8. Compilation of species annotation classifications
  
  #OTU corresponding species annotation in a 2-column format: remove confidence values from sintax, retain species annotations, replace ":" with "_", and remove quotation marks
  cut -f 1,4 result/otus.sintax \
|sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
> result/taxonomy2.txt
head -n3 result/taxonomy2.txt

#OTU corresponding species in an 8-column format: please note that annotations are not uniform
#Generate a species table where blanks in the OTU/ASV are filled with "Unassigned"
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
result/taxonomy2.txt > temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
> result/taxonomy.txt
head -n3 result/taxonomy.txt

#Count phylum, class, order, family, and genus using the rank parameters p, c, o, f, g respectively.
mkdir -p result/tax
for i in p c o f g;do
usearch -sintax_summary result/otus.sintax \
-otutabin result/otutab_rare.txt -rank ${i} \
-output result/tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
#List all files
wc -l result/tax/sum_*.txt
head -n3 result/tax/sum_g.txt

## 9. There is a reference-based quantitative feature table

# Perform alignment against the Greengenes 97% OTUs database for use in PICRUSt/Bugbase functional prediction
mkdir -p result/gg/
  
  #Method 1. Using USEARCH for alignment is faster, but encountering file size limit errors. Choose method 2
  # By default, utilize 1 core for systems with less than 10 cores, and allocate 10 cores for systems with 10 or more cores
  #search -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fa \
  #-otutabout result/gg/otutab.txt -threads 4
  # Alignment rate: 80.0%, 1 core took 11 minutes, 4 cores took 3 minutes, 10 cores took 2 minutes, and memory usage was 743MB
  head -n3 result/gg/otutab.txt

# #Method 2. vsearch alignment is more accurate but slower. However, it exhibits stronger performance with parallelization using 24-96 threads(used this for nanopore data)
vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fa \
--otutabout result/gg/otutab.txt --id 0.97 --threads 12
# sequences in temp/filtered.fa were successfully mapped against the 97_otus.fa database with a 97% identity threshold, and the OTU table was written to result/gg/otutab.txt
head -n3 result/gg/otutab.txt
#Statistics
usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
cat result/gg/otutab.stat


## 10. Space cleanup and data submission.

#Delete intermediate large files
rm -rf temp/*.fq

# Calculate MD5 values for paired-end data to use for data submission
cd seq
md5sum *_1.fq.gz > md5sum1.txt
md5sum *_2.fq.gz > md5sum2.txt
paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > ../result/md5sum.txt
rm md5sum*
  cd ..
cat result/md5sum.txt
  

## 1. Alpha diversity

### 1.1 Box plot of alpha diversity.

# View Help
Rscript ${db}/script/alpha_boxplot.R -h
# Full parameters, diversity index optional richness chao1 ACE shannon simpson invsimpson
Rscript ${db}/script/alpha_boxplot.R --alpha_index richness \
--input result/alpha/vegan.txt --design result/metadata.txt \
--group Group --output result/alpha/ \
--width 89 --height 59
# Use loops to plot 6 commonly used indices
for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
Rscript ${db}/script/alpha_boxplot.R --alpha_index ${i} \
--input result/alpha/vegan.txt --design result/metadata.txt \
--group Group --output result/alpha/ \
--width 89 --height 59
done
mv alpha_boxplot_TukeyHSD.txt result/alpha/
  
  # Alpha diversity histogram + standard deviation
  Rscript ${db}/script/alpha_barplot.R --alpha_index richness \
--input result/alpha/vegan.txt --design result/metadata.txt \
--group Group --output result/alpha/ \
--width 89 --height 59

### 1.2  Dilution curve

Rscript ${db}/script/alpha_rare_curve.R \
--input result/alpha/alpha_rare.txt --design result/metadata.txt \
--group Group --output result/alpha/ \
--width 120 --height 59

### 1.3 Diversity Venn diagram

# Three groups of comparisons: -f input file, -a/b/c/d/g group name, -w/u for width and height inches, -p output file name suffix
bash ${db}/script/sp_vennDiagram.sh \
-f result/alpha/otu_group_exist.txt \
-a HS -b NW -c YW  \
-w 3 -u 3 \
-p HS_NW_YW
# Four sets of comparisons, the diagram and code are shown in the input file directory, and the running directory is the root directory of the current project
bash ${db}/script/sp_vennDiagram.sh \
-f result/alpha/otu_group_exist.txt \
-a NW -b HS -c YW -d Zymo \
-w 3 -u 3 \
-p NW_HS_YW_Zymo
#  plots Venn diagrams online https://www.ehbio.com/test/venn

## 2. Beta diversity

### 2.1 Distance matrix heatmap 

# Taking bray_curtis as an example, -f input file, -h whether the clustering is TRUE/FALSE, -u/v is width and height inches
bash ${db}/script/sp_pheatmap.sh \
-f result/beta/bray_curtis.txt \
-H 'TRUE' -u 6 -v 5
# Add grouping notes such as genotype and location for columns 2,4
cut -f 1-2 result/metadata.txt > temp/group.txt
# -P to add row comment files, -Q to add column comments
bash ${db}/script/sp_pheatmap.sh \
-f result/beta/bray_curtis.txt \
-H 'TRUE' -u 6.9 -v 5.6 \
-P temp/group.txt -Q temp/group.txt
# The distance matrix is similar to correlation, try corrplot or ggcorrplot to plot more styles
# - [Plot the correlation coefficient matrix corrplot](http://mp.weixin.qq.com/s/H4_2_vb2w_njxPziDzV4HQ)
# - [Correlation matrix visualization ggcorrplot](http://mp.weixin.qq.com/s/AEfPqWO3S0mRnDZ_Ws9fnw)

### 2.2 Primary coordinate analysis PCoA

# Input file, choose grouping, output file, image size in mm, and statistical results in beta_pcoa_stat.txt
Rscript ${db}/script/beta_pcoa.R \
--input result/beta/bray_curtis.txt --design result/metadata.txt \
--group Group --label FALSE --width 89 --height 59 \
--output result/beta/bray_curtis.pcoa.pdf
# Add sample labels with the parameter --label TRUE
Rscript ${db}/script/beta_pcoa.R \
--input result/beta/bray_curtis.txt --design result/metadata.txt \
--group Group --label TRUE --width 89 --height 59 \
--output result/beta/bray_curtis.pcoa.label.pdf
mv beta_pcoa_stat.txt result/beta/
  
  ### 2.3 Constrained Principal Coordinates Analysis CPCoA
  
  Rscript ${db}/script/beta_cpcoa.R \
--input result/beta/bray_curtis.txt --design result/metadata.txt \
--group Group --output result/beta/bray_curtis.cpcoa.pdf \
--width 89 --height 59
# Add sample labels with the parameter --label TRUE
Rscript ${db}/script/beta_cpcoa.R \
--input result/beta/bray_curtis.txt --design result/metadata.txt \
--group Group --label TRUE --width 89 --height 59 \
--output result/beta/bray_curtis.cpcoa.label.pdf

## 3. Species composition Taxonomy

### 3.1 Stacked bar chart Stackplot

# Taking the phylum (p) level as an example, the results will include two files: "output.sample.pdf" and "output.group.pdf"
Rscript ${db}/script/tax_stackplot.R \
--input result/tax/sum_p.txt --design result/metadata.txt \
--group Group --color ggplot --legend 7 --width 89 --height 69 \
--output result/tax/sum_p.stackplot
# Modify colors using ggplot manual1 (22), Paired (12), or Set3 (12) color palettes
Rscript ${db}/script/tax_stackplot.R \
--input result/tax/sum_p.txt --design result/metadata.txt \
--group Group --color Paired --legend 12 --width 181 --height 119 \
--output result/tax/sum_p.stackplotPaired

# Batch draw using input that includes phylum (p), class (c), order (o), family (f), and genus (g) for a total of 5 levels
for i in p c o f g; do
Rscript ${db}/script/tax_stackplot.R \
--input result/tax/sum_${i}.txt --design result/metadata.txt \
--group Group --output result/tax/sum_${i}.stackplot \
--legend 8 --width 89 --height 69; done

### 3.2 Chord/circular diagram using Circlize

# Taking class (c) as an example, draw the top 5 groups.
i=c
Rscript ${db}/script/tax_circlize.R \
--input result/tax/sum_${i}.txt --design result/metadata.txt \
--group Group --legend 5
# The results are located in the current directory: "circlize.pdf" (with random colors) and "circlize_legend.pdf" (with specified colors and legend)
# Move and rename to match the classification level
mv circlize.pdf result/tax/sum_${i}.circlize.pdf
mv circlize_legend.pdf result/tax/sum_${i}.circlize_legend.pdf

### 3.3 Treemap visualization (for reference)

# Hierarchical treemap visualization depicting species relationships. Input feature table and species annotation, output treemap
# Specify the number of features and image dimensions. Generating the treemap for 100 ASVs took 12 seconds
Rscript ${db}/script/tax_maptree.R \
--input result/otutab.txt --taxonomy result/taxonomy.txt \
--output result/tax/tax_maptree.pdf \
--topN 100 --width 183 --height 118


# 24、Comparison of differences

## 1. R language difference analysis

### 1.1 Comparison of differences 

mkdir -p result/compare/
  # Input feature table, metadata; Specify the grouping column name, comparison group, and abundance
  # Select the methods wilcox/t.test/edgeR, pvalue, and fdr and output directories
  compare="HS-NW"
Rscript ${db}/script/compare.R \
--input result/otutab.txt --design result/metadata.txt \
--group Group --compare ${compare} --threshold 0.1 \
--method edgeR --pvalue 0.05 --fdr 0.2 \
--output result/compare/
  # Common error: Error in file(file, ifelse(append, "a", "w")): Unable to open link Calls: write.table -> file
  # Solution: The output directory does not exist, just create a directory

### 1.2 Volcano map
  
  #Enter compare. As a result of R, output volcano map with data labels, you can specify the image size
  Rscript ${db}/script/compare_volcano.R \
--input result/compare/${compare}.txt \
--output result/compare/${compare}.volcano.pdf \
--width 89 --height 59

### 1.3 Heat map

# Enter compare. R results, filter the number of columns, specify metadata and grouping, species annotations, figure size inches and font size
bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
-d result/metadata.txt -A Group \
-t result/taxonomy.txt \
-w 8 -h 5 -s 7 \
-o result/compare/${compare}


### 1.4 Map of manhattan

# idifference comparison results, t species annotation, p legend, w width, v height, s size, l legend maximum
# The legend shows no figure, you can increase the height v to 119+, and the AI puzzle is KO-WT.heatmap.emf in the later stage
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_p.txt \
-w 183 -v 109 -s 7 -l 10 \
-o result/compare/${compare}.manhattan.p.pdf

# There are only 6 gates in the picture above, switching to class c and -L class to show the details
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_c.txt \
-w 183 -v 149 -s 7 -l 10 -L Class \
-o result/compare/${compare}.manhattan.c.pdf
# Show the full legend and use AI puzzles
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_c.txt \
-w 183 -v 149 -s 7 -l 10 -L Class \
-o result/compare/${compare}.manhattan.c.legend.pdf


# Input feature table, metadata; Specify the grouping column name, comparison group, and abundance
# Select the methods wilcox/t.test/edgeR, pvalue, and fdr and output directories
compare="NW-YW"
Rscript ${db}/script/compare.R \
--input result/otutab.txt --design result/metadata.txt \
--group Group --compare ${compare} --threshold 0.1 \
--method edgeR --pvalue 0.05 --fdr 0.2 \
--output result/compare/
  # Common error: Error in file(file, ifelse(append, "a", "w")): Unable to open link Calls: write.table -> file
  # Solution: The output directory does not exist, just create a directory
  
  ### 1.2 Volcano map
  
  #Enter compare. As a result of R, output volcano map with data labels, you can specify the image size
  Rscript ${db}/script/compare_volcano.R \
--input result/compare/${compare}.txt \
--output result/compare/${compare}.volcano.pdf \
--width 89 --height 59

### 1.3 Heat map

# Enter compare. R results, filter the number of columns, specify metadata and grouping, species annotations, figure size inches and font size
bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
-d result/metadata.txt -A Group \
-t result/taxonomy.txt \
-w 8 -h 5 -s 7 \
-o result/compare/${compare}


### 1.4 Map of manhattan

# idifference comparison results, t species annotation, p legend, w width, v height, s size, l legend maximum
# The legend shows no figure, you can increase the height v to 119+, and the AI puzzle is KO-WT.heatmap.emf in the later stage
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_p.txt \
-w 183 -v 109 -s 7 -l 10 \
-o result/compare/${compare}.manhattan.p.pdf

# There are only 6 gates in the picture above, switching to class c and -L class to show the details
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_c.txt \
-w 183 -v 149 -s 7 -l 10 -L Class \
-o result/compare/${compare}.manhattan.c.pdf
# Show the full legend and use AI puzzles
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_c.txt \
-w 183 -v 149 -s 7 -l 10 -L Class \
-o result/compare/${compare}.manhattan.c.legend.pdf




# Input feature table, metadata; Specify the grouping column name, comparison group, and abundance
# Select the methods wilcox/t.test/edgeR, pvalue, and fdr and output directories
compare="YW-Zymo"
Rscript ${db}/script/compare.R \
--input result/otutab.txt --design result/metadata.txt \
--group Group --compare ${compare} --threshold 0.1 \
--method edgeR --pvalue 0.05 --fdr 0.2 \
--output result/compare/
  # Common error: Error in file(file, ifelse(append, "a", "w")): Unable to open link Calls: write.table -> file
  # Solution: The output directory does not exist, just create a directory
  
  ### 1.2 Volcano map
  
  #Enter compare. As a result of R, output volcano map with data labels, you can specify the image size
  Rscript ${db}/script/compare_volcano.R \
--input result/compare/${compare}.txt \
--output result/compare/${compare}.volcano.pdf \
--width 89 --height 59

### 1.3 Heat map

# Enter compare. R results, filter the number of columns, specify metadata and grouping, species annotations, figure size inches and font size
bash ${db}/script/compare_heatmap.sh -i result/compare/${compare}.txt -l 7 \
-d result/metadata.txt -A Group \
-t result/taxonomy.txt \
-w 8 -h 5 -s 7 \
-o result/compare/${compare}


### 1.4 Map of manhattan

# idifference comparison results, t species annotation, p legend, w width, v height, s size, l legend maximum
# The legend shows no figure, you can increase the height v to 119+, and the AI puzzle is KO-WT.heatmap.emf in the later stage
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_p.txt \
-w 183 -v 109 -s 7 -l 10 \
-o result/compare/${compare}.manhattan.p.pdf

# There are only 6 gates in the picture above, switching to class c and -L class to show the details
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_c.txt \
-w 183 -v 149 -s 7 -l 10 -L Class \
-o result/compare/${compare}.manhattan.c.pdf
# Show the full legend and use AI puzzles
bash ${db}/script/compare_manhattan.sh -i result/compare/${compare}.txt \
-t result/taxonomy.txt \
-p result/tax/sum_c.txt \
-w 183 -v 149 -s 7 -l 10 -L Class \
-o result/compare/${compare}.manhattan.c.legend.pdf


### 1.5Drawing of individual features

# Differential ASV was screened and displayed in descending order by abundance of KO group, and the top 10 were displayed by ID
awk '$4<0.05' result/compare/HS-NW.txt | sort -k7,7nr | cut -f1 | head
# Differential OTU details display
Rscript ${db}/script/alpha_boxplot.R --alpha_index ASV_2 \
--input result/otutab.txt --design result/metadata.txt \
--transpose TRUE --scale TRUE \
--width 89 --height 59 \
--group Group --output result/compare/feature_ 
# If the ID does not exist, an error will be reported： Error in data.frame(..., check.names = FALSE) : Parameter values mean different number of rows: 0, 18  Calls: alpha_boxplot -> cbind -> cbind -> data.frame
# Specify a column sort: All descending by mean genus abundance
csvtk -t sort -k All:nr result/tax/sum_g.txt | head
# Poor details are shown
Rscript ${db}/script/alpha_boxplot.R --alpha_index Lysobacter \
--input result/tax/sum_g.txt --design result/metadata.txt \
--transpose TRUE \
--width 89 --height 59 \
--group Group --output result/compare/feature_

### 1.5 Ternary Diagram

#For reference examples, see: result\compare\ternary\ternary. RMD documentation
#Alternative Tutorial[246.Application and Drawing Practice of Ternary Diagrams](https://mp.weixin.qq.com/s/3w3ncpwjQaMRtmIOtr2Jvw)




Ternary Plot with ggtern 
install.packages("ggtern")
library(ggtern)








##5. Generate OTU Table

##6. Statistical Analysis (Alpha, beta, taxonomic composition)

##7. Visualization and Reporting







































##Radarplot
Rscript ${db}/script/radar_plot3.R \
-i result/tax/sum_p.txt \
-d result/metadata.txt \
-g Group \
-t Genus \
-n 10 \
-o result/compare/radar_plot.pdf




































##Taxonomic classification with emu (run successfully but found minor issue due to simultaneously running some function in two terminals, trying again )

#Install osfclient: First, if you haven't already, you need to install osfclient (OSF command line tool) which is required for downloading the pre-built databases.
pip install osfclient
#. Define the Database Variable: Set the EMU_PREBUILT_DB variable to silva, as you're using the SILVA 138.1 database.
export EMU_PREBUILT_DB='silva'

#3. Define the Database Directory
export EMU_DATABASE_DIR='/mnt/d/Amplicon2/emu_databases'
cd ${EMU_DATABASE_DIR}
osf -p 56uf7 fetch osfstorage/emu-prebuilt/${EMU_PREBUILT_DB}.tar
tar -xvf ${EMU_PREBUILT_DB}.tar
ls -lh ${EMU_DATABASE_DIR}

export EMU_PREBUILT_DB='silva'
export EMU_DATABASE_DIR='/mnt/d/Amplicon2/emu_databases'

#4. Navigate to the Database Directory
cd ${EMU_DATABASE_DIR}

#Download the Database Using osfclient
osf -p 56uf7 fetch osfstorage/emu-prebuilt/${EMU_PREBUILT_DB}.tar

6. Extract the Tarball
tar -xvf ${EMU_PREBUILT_DB}.tar

#Verify the Contents of the Extracted Directory
ls -lh ${EMU_DATABASE_DIR}

#Use Emu for Taxonomic Classification: Replace input_sequences.fasta with the path to your actual sequence file.
INPUT_SEQUENCE_FILE="/mnt/d/Amplicon2/temp/medaka_output/consensus.fa"
OUTPUT_CLASSIFIED="classified_output"

emu classify -i "${INPUT_SEQUENCE_FILE}" -o "${OUTPUT_CLASSIFIED}" --db "${EMU_DATABASE_DIR}/silva-138"

OR 

emu classify -i /mnt/d/Amplicon2/temp/medaka_output/consensus.fa -o classified_output --db /mnt/d/Amplicon2/emu_databases/silva-138

#Check the Output: The classified_output directory will contain the results of the taxonomic classification. This should include a summary of classifications at different taxonomic levels (e.g., kingdom, phylum, class, etc.).

##Generate Initial Abundance Table:
# Run abundance profiling
emu abundance -i classified_output/classified_result.tsv -o abundance_output

## Collapse Taxonomy to Genus Level:

emu collapse-taxonomy classified_output/classified_result.tsv genus

##Collapse Taxonomy to Phylum Level:
emu collapse-taxonomy classified_output/classified_result.tsv phylum



