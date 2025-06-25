
EasyAmplicon V.2: This pipeline processes long-read amplicon sequencing data (Nanopore/PacBio) from raw reads to polished OTUs/ASVs with taxonomic assignment. 
It includes data preprocessing, quality control, dereplication, clustering/denoising, alignment and polishing, chimera removal, and taxonomic classification.

1. Setup and Initialization
1.1 Environment and Directories
#Install and configure Miniconda for package management.
# Conda software installation directory, check with `conda env list`, e.g., /anaconda3  
soft=~/miniconda3  
# Working directory (wd) 
wd="/mnt/d/Amplicon"
# Database directory (db) 
db="/mnt/d/EasyMicrobiome2"
# Create working directories
mkdir -p $wd/seq $wd/result $wd/temp
cd $wd
# ├── pipeline.sh
# ├── result
# │   └── metadata.txt
# ├── seq
# │   ├── HS_1.fq.gz
# │   ├── ...
# │   └── ZYmo_2.fq.gz
# └── temp

# Verify the directories
echo "Working directory: $wd"
echo "Database directory: $db"

1.2. Metadata Preparation

# Prepare and save sample metadata in "result/metadata.txt"
# Use "csvtk" to count the number of rows (samples, excluding the header) and columns in a table. Use the -t option to set the column delimiter as a tab (default is comma);
csvtk -t stat result/metadata_raw.txt
# The metadata should have at least 3 columns. The first column should be the sample ID (SampleID), and the last column should be the description (Description)
# You can use the cat command to view a file, and the -A option to display symbols. The '|' symbol is used as a pipe to connect commands. The head command is used to display the file header, and the -n3 option limits the output to the first 3 lines.
cat -A result/metadata_raw.txt | head -n3
# Remove Windows line endings if present with "sed" command, and then check with cat -A
sed 's/\r//' result/metadata_raw.txt > result/metadata.txt
cat -A result/metadata.txt | head -n3

1.3. Sequencing Data Preparation
#Place raw sequencing .fastq.gz files in seq/.
#The sequencing results returned by the company are usually single FastQ compressed files in .gz format for one sample, for nanopore data and incase of pacbio, its pairend data. 
# The file name and the sample name must be corresponding and can be manually modify in case of inconsistency.
# If the sequencing data is a compressed file of .gz, sometimes you need to use gunzip to decompress and use, vsearch can usually read the compressed file directly
# gunzip seq/*.gz
# Verify ses folder and downloaded sequencing files
ls -sh seq/
# View sequence file content, specify character range per line
zless seq/HS_1.fastq | head -n4
zless seq/HS_1.fastq | head | cut -c 1-60
#Check sequencing file statistics using seqkit
seqkit stat seq/HS_1.fastq
## Summarize sequencing stats for all samples
seqkit stat seq/*.fastq > result/seqkit.txt
head result/seqkit.txt

1.4.DatabasesPreparation 
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

# To decompress the "ITS" UNITE database, you need to download it from the official website or the network disk db/amplicon/usearch
# gunzip -c ${db}/usearch/utax_reference_dataset_all_29.11.2022.fasta.gz >${db}/usearch/unite.fa
seqkit stat ${db}/usearch/unite.fa # 326 thousand
# The Greengene database is used for feature annotations: wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
# The default decompression will delete the original file, -c specifies the output to the screen, > Write a new file (can be renamed):
gunzip -c ${db}/gg/97_otus.fasta.gz > ${db}/gg/97_otus.fa
seqkit stat ${db}/gg/97_otus.fa

2. Data Preprocessing
2.1. Quality Assessment with NanoPlot
##Install Nanoplot
conda install -c bioconda nanoplot
## If already installed then activate nanoplot_env
conda activate nanoplot_env 
## Run NanoPlot on each sample FASTQ:
mkdir -p nanoplot_reports  
INPUT_DIR="/mnt/d/Amplicon2/Seq"
OUTPUT_DIR="/mnt/d/Amplicon2/nanoplot_reports"
mkdir -p "$OUTPUT_DIR"
SAMPLES=("HS_1" "HS_2" "HS_3" "NW_1" "NW_2" "NW_3" "YW_1" "YW_2" "YW_3" "Zymo_1" "Zymo_2" "Zymo_3")
# Run NanoPlot
for SAMPLE in "${SAMPLES[@]}"; do
echo "Processing $SAMPLE..."
NanoPlot --fastq "$INPUT_DIR/${SAMPLE}.fastq" \
--outdir "$OUTPUT_DIR/${SAMPLE}_report/" \
--plots hex dot --loglength --N50
done

2.2. Reads Rename
#Relabel reads to include sample prefix for traceability:
time for i in $(tail -n +2 result/metadata.txt | cut -f1); do
usearch -fastx_relabel "seq/${i}.fastq" -fastqout "temp/${i}.fq" -prefix "${i}."
done
#Merge all renamed reads:
cat temp/*.fq > temp/all.fq
#View file size: 2.73GB. There are slight differences in results due to different software versions.
ls -lsh temp/all.fq
# Check if the sequence name before the '.' is the sample name. Sample names must not contain dots ('.')
# One notable characteristic of sample names with dots (.) is that the generated feature table becomes significantly large, with numerous columns, leading to insufficient memory during subsequent analyses。
# After obtaining the feature table from the subsequent analysis, it's important to take a look to check for any issues. If you encounter memory-related problems, you'll need to go back and troubleshoot。
head -n 6 temp/all.fq|cut -c1-60

3. Primer Trimming and Quality Filtering
3.1. Trim primers and length filter
#Use cutadapt to trim primers and filter by length
#In result, total reads processed	1,050,569, Reads with adapters	789,024 (75.1%), Reverse-complemented	399,680 (38.0%),Too short (< 1000)	134,153 (12.8%), Too long (> 1800)	3,355 (0.3%), Discarded as untrimmed	129,432 (12.3%), Final retained reads	783,629 (74.6%). 
cutadapt -g "AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA" \
--error-rate=0.1 \
--action=trim \
--rc \
-m 1000 \
-M 1800 \
-j 20 \
--discard-untrimmed \
-o /mnt/d/Amplicon/temp/all.filtered.fastq \
/mnt/d/Amplicon/temp/all.fq
 
3.2. Quality check with FastQC
#Run FastQC on trimmed reads:
mkdir -p /mnt/d/Amplicon/temp/qc
fastqc -o /mnt/d/Amplicon/temp/qc /mnt/d/Amplicon/temp/all.filtered.fastq

3.3. Quality filtering with NanoFilt
#Filter reads by length ≥1000, quality ≥20, headcrop 10 bases, max length 1800:
#Installing method
#conda install -c bioconda nanofilt
#pip install nanofilt
#pip install nanofilt --upgrad

cat /mnt/d/Amplicon/temp/all.filtered.fastq | NanoFilt -l 1000 -q 20 --headcrop 10 --maxlength 1800 > /mnt/d/Amplicon/temp/filtered.fa
sed 's/ rc$//' > /mnt/d/Amplicon/temp/filtered.fa
##Or user can use this code as well: 
cat /mnt/d/Amplicon/temp/all.filtered.fastq | NanoFilt -l 1000 -q 20 --headcrop 10 --maxlength 1800 | sed 's/ rc$//' > /mnt/d/Amplicon/temp/filtered.fa
## For emu pipline: we should have Quality filtering with NanoFilt  and keep FASTQ
cat /mnt/d/Amplicon/temp/all.filtered.fastq | NanoFilt -l 1000 -q 20 --headcrop 10 --maxlength 1800 > /mnt/d/Amplicon/temp/filtered.fastq

##Step 4: FastQC on cleaned reads
fastqc -o /mnt/d/Amplicon/temp/qc /mnt/d/Amplicon/temp/fileterd.fa

4. Dereplication and OTU/ASV Selection
4.1.Dereplication
#Remove redundant sequences with minimum unique size threshold:
# And add a miniuniqusize of at least 2 or 1/1M to remove low-abundance noise and increase calculation speed. Or 1/1M means that alternatively, a feature should have a frequency of at least 1 per million (1/1M) to be retained. This is useful when dealing with large datasets where low-abundance sequences may be treated as noise.
# -sizeout outputs abundances, and --relabel must be used with sequence prefixes for better organization. Takes 1 second.
# size="121" with miniuniquesize 10-, and "1990" with miniuniquesize 2-. Uni_1 is the name of the sequence after dereplication. This sequence appeared 590 times in the sequencing data across all samples.
vsearch --derep_fulllength temp/filtered.fa \
--minuniquesize 2 --sizeout --relabel Uni_ \
--output temp/uniques2.fa

# Uni_1 is the sequence that appears the most。
#Check output summary
ls -lsh temp/uniques2.fa
head -n 2 temp/uniques.fa
##Count the number of sequences:(121)
grep -c "^>" temp/uniques.fa
##Count the number of sequences:(1990)
grep -c "^>" temp/uniques2.fa

4.2. Cluster OTUs or denoise ASVs

#There are two methods: It is recommended to use unoise3 for denoising to obtain ASVs with single-base accuracy. Alternatively, you can use the traditional 97% clustering for OTUs (genus-level accuracy) as an alternative
#Both feature selection methods in usearch come with built-in de novo chimera removal.
#-minsize secondary filtering, control the number of OTU/ASV to 1-5,000, convenient for downstream statistical analysis

#Method 1. ASV Denoise: predict biological sequences and filter chimeras, prefer for high-resolution diversity (species/strain level)
#In 6 seconds, with "minsize 10" 107 sequences were identified as good, , with 1 chimeras removed. While with "minsize 2" obtained 818. 
usearch -unoise3 temp/uniques2.fa -minsize 10 \
-zotus temp/zotus.fa

usearch -unoise3 temp/uniques2.fa -minsize 2 \
-zotus temp/zotus.fa

#Modify sequence names: Change "Zotu" to "ASV" for easier identification
sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
head -n 2 temp/otus.fa
#Count how many ASVs you got (818)
grep -c "^>" /mnt/d/Amplicon/temp/otus.fa

## Method 2: 97% clustering for OTUs, suitable for large datasets/when ASV patterns are not clear/reviewer's request. Prefer for broad taxonomic trends (genus level)
# Rename using "relabel" and cluster at a similarity threshold of 97% without masking qmask
# Record the input size ("sizein") and the output frequency ("sizeout"): the output is Reading file temp/uniques.fa 100%, 174752 nt in 121 seqs, min 1377, max 1498, avg 1444, Sorting by abundance 100%, Counting k-mers 100%, Clustering 100%, Sorting clusters 100%, Writing clusters 100%, Clusters: 46 Size min 10, max 1029, avg 2.6, Singletons: 0, 0.0% of seqs, 0.0% of clusters
###Method:2 OTU clustering at 97% identity
vsearch --cluster_size temp/uniques.fa \
--relabel ASV_ \
--id 0.97 \
--qmask none \
--sizein \
--sizeout \
--centroids temp/otus_raw.fa

#Clean sequence headers (remove ";size=...")
sed 's/;size=[0-9]*//' temp/otus_raw.fa > temp/otus2.fa
#Count how many ASVs you got (46)
grep -c "^>" /mnt/d/Amplicon/temp/otus2.fa
# Preview the output
head -n 4 temp/otus2.fa
#In order to check the formatting of database 
grep "^>" silva_16s_v138.fa | head -n 10
# first sequence line after the header in the file
grep -v "^>" silva_16s_v138.fa | head -n 5
# with format databse where change nucelotide base according to DNA. 
grep "^>" silva_16s_v138_dna.fa | head -n 10
# first sequence line after the header in the file
grep -v "^>" silva_16s_v138_dna.fa | head -n 5

5. Aligment with minimap2 and plosihing with meadaka 
5.1. Allignment with minimap2: 
#For Installation, Minimap2 is optimized for x86-64 CPUs.(if not already installed), Skip if already installed.

curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
  ./minimap2-2.28_x64-linux/minimap2
# Set the database directory path (adjust as needed)
db="/mnt/d/EasyMicrobiome2/usearch"
# Move to database directory
cd $db
minimap2 -ax map-ont temp/otus.fa temp/all.fq > temp/aligned_output.sam
# Set working directory
cd /mnt/d/Amplicon
# Align reads to draft OTUs
minimap2 -ax map-ont temp/otus.fa temp/all.fq > temp/aligned_output.sam
# Convert SAM to BAM, sort, and index
samtools view -Sb temp/aligned_output.sam > temp/aligned_output.bam
samtools sort temp/aligned_output.bam -o temp/aligned_output_sorted.bam
samtools index temp/aligned_output_sorted.bam

5.2. Polishing with RACON 
# Apply RACON to polish OTUs
racon -m 8 -x -6 -g -8 -w 500 -t 14 temp/all.fq temp/aligned_output.sam temp/otus.fa > temp/otus_racon.fa
##Step 3: Prepare for MEDAKA
# Index FASTA (optional for medaka, helpful for inspection)
samtools faidx temp/otus_racon.fa

5.3. Polishing with  Medaka (replace 'r941_min_sup_g507' with correct model for your flowcell + basecaller)
#Install and activate Medaka
conda install -c bioconda medaka
conda activate medaka
#Run Medaka consensus
medaka_consensus -i temp/all.fq \
-d temp/otus_racon.fa \
-o temp/medaka_polished \
-t 14 \
-m r941_min_sup_g507
#Rename consensus file, consensus.fasta to consensus.fa
cd /mnt/d/Amplicon/temp/medaka_polished
mv "consensus.fasta" "consensus.fa"
#If size number present in consesus.fasta, remove using this code
sed 's/;size=[0-9]*//' consensus.fasta > consensus.fa
conda deactivate
####6 Perform reference-based chimera detection.
# Removing all OTU_1 sequences seems reasonable since Usearch does not have built-in chimera removal methods。
cd /mnt/d/Amplicon
mkdir -p result/raw

6. Chimera Detection and Removal
#Reference-based chimera removal using SILVA "silva_16s_v138_dna.fa'.Taking abundance information into account, this corresponds to 146 (17.8%) chimeras, 585 (71.5%) non-chimeras,and 87 (10.6%) borderline sequences in 818 total sequences
db="/mnt/d/EasyMicrobiome2"
cd /mnt/d/Amplicon
vsearch --uchime_ref temp/medaka_polished/consensus.fa \
-db ${db}/usearch/silva_16s_v138_dna.fa \
--nonchimeras result/raw/otus.fa \
--chimeras result/raw/otus.chimeras.fa \
--borderline result/raw/otus.borderline.fa
## If you're using Windows and vsearch results have added Windows line endings (^M), you need to remove them. However, if you're using a Mac, you don't need to execute this command
sed -i 's/\r//g' result/raw/otus.fa

7. Feature Table and Taxonomy (Taxonomic Classification: Assign taxonomy to OTUs or ASVs).
##change terminal here from WSL to gitbash
7.1. Generate a Feature table 
vsearch --usearch_global temp/filtered.fasta \
--db result/raw/otus.fa \
--id 0.97 \
--threads 4 \
--otutabout result/raw/otutab.txt

#In case, If you met this error "Fatal error: FASTA file expected, FASTQ file found (temp/filtered.fa)", than run this code 
vsearch --fastq_filter temp/filtered.fa \
--fastaout temp/filtered.fasta \
--fastq_qmax 50

# For Windows users, you should remove the newline characters (^M) from vsearch results and correct them to the standard Linux format.
sed -i 's/\r//' result/raw/otutab.txt
head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
# To count rows and columns in a table using "csvtk"
# Here, it's important to carefully check the number of columns and whether it matches the number of your samples. If they're not equal, there might be issues with sample naming. Refer to the explanations above for more details.
csvtk -t stat result/raw/otutab.txt

7.2. Generate a Feature table Taxonomy Assignment (Use the following command to assign taxonomy to your OTUs using the SILVA database:)
# Use this code with latest version of Silva database "SILVA_138.2_SSURef_NR99_tax_silva.fa" formatted version "silva_reformatted.fa"
vsearch --sintax result/raw/otus.fa \
--db ${db}/usearch/silva_16s_v138_dna.fa \
--sintax_cutoff 0.8 \
--tabbedout result/raw/otus.sintax
head result/raw/otus.sintax | cat -A
sed -i 's/\r//' result/raw/otus.sintax

# Optional: Use RDB silva database 
vsearch --sintax result/raw/otus.fa \
--db ${db}/usearch/rdp_16s_v18.fa \
--sintax_cutoff 0.8 \
--tabbedout result/raw/otus.sintax 
head result/raw/otus.sintax | cat -A
sed -i 's/\r//' result/raw/otus.sintax

# Optional: To run EMU taxonomy annotation instead of or alongside SILVA sintax
# Assuming EMU is installed and configured
# Run EMU taxonomy annotation
emu taxonomy classify --input-fasta result/raw/otus.fa --output-tsv result/raw/otus_emu_taxonomy.tsv --db emu_db_path

7.3. Number of original feature table rows
wc -l result/raw/otutab.txt
#R script selects bacterial archaea (eukaryotes), removes chloroplasts, mitochondria and counts the proportions; Output filtered and sorted OTU tables
#The input is the OTU table result/raw/otutab.txt and the species annotation result/raw/otus.sintax
#Output filtered and sorted feature tables result/otutab.txt sum
#Statistical contamination ratio file result/raw/otutab_nonBac.txt and filter details otus.sintax.discard
#For fungal ITS data, use the otutab_filter_nonFungi.R script instead, only screening fungi
# Rscript ${db}/script/otutab_filter_nonBac.R -h # Displays a description of the parameter
Rscript ${db}/script/otutab_filter_nonBac.R \
--input result/raw/otutab.txt \
--taxonomy result/raw/otus.sintax \
--output result/otutab.txt \
--stat result/raw/otutab_nonBac.stat \
--discard result/raw/otus.sintax.discard
# The number of filtered feature table rows
wc -l result/otutab.txt
#Filter the sequence corresponding to the feature table
cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
usearch -fastx_getseqs result/raw/otus.fa \
-labels result/otutab.id -fastaout result/otus.fa
#Filter feature tables corresponding to sequence annotations
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' \
result/raw/otus.sintax result/otutab.id > result/otus.sintax
# Optional statistical method: Summary OTUs table
usearch -otutab_stats result/otutab.txt \
-output result/otutab.stat
cat result/otutab.stat
#Note the minimum, quantile, or look at the amount of sample detail in result/raw/otutab_nonBac.stat for resampling

8. Normalization (Equalize Sampling Depth) 
# change terminal here from WSL to gitbash
# Open a new terminal, and Run this code 
wd=/d/Amplicon
db=/d/EasyMicrobiome2
PATH=$PATH:${db}/win
cd ${wd}

#Use the vegan package for equal resampling and enter the reads count format Feature table result/otutab.txt. Adjust the "depth 10000" as per using the sample minimum size, here its 1675, so we are using "0". 
#You can specify the input file, sample size, and random number, and output the equalizer result/otutab_rare.txt and diversity alpha/vegan .txt
mkdir -p result/alpha
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
--depth 0 --seed 1 \
--normalize result/otutab_rare.txt \
--output result/alpha/vegan.txt
usearch -otutab_stats result/otutab_rare.txt \
-output result/otutab_rare.stat
cat result/otutab_rare.stat

### Need change the depth with "1675". 

                             9. alpha (α) diversity

9.1. Calculate alpha (α) diversity
# Use USEARCH to calculate 14 alpha diversity indices (don't use Chao1 if you make a mistake
#details in http://www.drive5.com/usearch/manual/alpha_metrics.html
usearch -alpha_div result/otutab_rare.txt \
-output result/alpha/alpha.txt

9.2. Calculate rarefaction richness

#Dilution Curve: Take the number of OTUs (Operational Taxonomic Units) from sequences ranging from 1% to 100%, with non-replacement sampling each time.
#Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch -alpha_div_rare result/otutab_rare.txt \
-output result/alpha/alpha_rare.txt \
-method without_replacement
#Preview of Results
head -n2 result/alpha/alpha_rare.txt
#Handling of non-numeric "-" values due to low sample sequencing, see FAQ #8 for details
sed -i "s/-/\t0.0/g" result/alpha/alpha_rare.txt

9.3. Filtering for high-abundance bacteria
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


9.4. Beta(β) diversity
# The results generate multiple files and require a directory (replace otus.fa with otus.nonchimeras.fa )
mkdir -p result/beta/
# Building an evolutionary tree based on OTUs
usearch -cluster_agg result/otus.fa -treeout result/otus.tree
# Generate five distance matrices：bray_curtis, euclidean, jaccard, manhatten, unifrac
usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
-filename_prefix result/beta/
  
# While ruuning above code, user may encounter fatal error and alignment issue, like this " Align your sequences using MAFFT (to ensure they are aligned properly before building a tree).
# Use FastTree to build the evolutionary tree from the aligned sequences.
# Compute distance matrices using QIIME 2 or R packages (e.g., microeco) if needed.
# Download and install fastree. FastTree infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences. FastTree can handle alignments with up to a million of sequences in a reasonable amount of time and memory. 
# Save the file as FastTree.exe in a directory you can easily access, such as D:\EasyAmplicon or a dedicated tools directory like D:\software.Then add FastTree to PATH.Open Open Environment Variables, Under System Variables, find Path and click Edit. Add the directory containing FastTree.exe (e.g., D:\EasyAmplicon) to the list
# Save and restart your terminal. Then run the following command to confirm installation (FastTree -help). If the installation is successful, you should see the help message for FastTree
  
FastTree -nt result/otus.fa > result/otus.tree

## The error Wrong number of characters for ASV_1: expected 1646 but have 1464 instead indicates that there is an inconsistency in sequence lengths in the result/otus.nonchimeras.fa file. FastTree requires all sequences to have the same length, as it expects an aligned sequence file.  
## For alignment user can choose tools like MAFFT or Clustal Omega. 
## while after downloading, installing window 
iqtree -s result/otus.fa -nt AUTO -pre result/otus

10. Compilation of species annotation classifications

#OTU corresponding species annotation in a 2-column format: remove confidence values from sintax, retain species annotations, replace ":" with "_", and remove quotation marks
cut -f 1,4 result/otus.sintax \
|sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
> result/taxonomy.txt
head -n3 result/taxonomy.txt

#OTU corresponding species in an 8-column format: please note that annotations are not uniform
#Generate a species table where blanks in the OTU/ASV are filled with "Unassigned"
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
result/taxonomy.txt > temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
> result/taxonomy.txt
head -n3 result/taxonomy.txt

#Count phylum, class, order, family, and genus using the rank parameters p, c, o, f, g respectively.
mkdir -p result/tax
for i in p c o f g s;do
usearch -sintax_summary result/otus.sintax \
-otutabin result/otutab_rare.txt -rank ${i} \
-output result/tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
#List all files
wc -l result/tax/sum_*.txt
head -n3 result/tax/sum_g.txt

11. Reference-based quantitative feature table

# Perform alignment against the Greengenes 97% OTUs database for use in PICRUSt/Bugbase functional prediction
mkdir -p result/gg/
  
# Method 1. Using USEARCH for alignment is faster, but encountering file size limit errors. Choose method 2
# By default, utilize 1 core for systems with less than 10 cores, and allocate 10 cores for systems with 10 or more cores
# search -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fa \
# -otutabout result/gg/otutab.txt -threads 4
# Alignment rate: 80.0%, 1 core took 11 minutes, 4 cores took 3 minutes, 10 cores took 2 minutes, and memory usage was 743MB
 head -n3 result/gg/otutab.txt

# Method 2. vsearch alignment is more accurate but slower. However, it exhibits stronger performance with parallelization using 24-96 threads(used this for nanopore data)
# Matching unique query sequences: 238661 of 720330 (33.13%)
vsearch --usearch_global temp/filtered.fasta --db ${db}/gg/97_otus.fa \
--otutabout result/gg/otutab.txt --id 0.97 --threads 12
# sequences in temp/filtered.fa were successfully mapped against the 97_otus.fa database with a 97% identity threshold, and the OTU table was written to result/gg/otutab.txt
head -n3 result/gg/otutab.txt
#Statistics (Sample sizes: min 13207, lo 16064, med 19653, mean 19888.4, hi 22563, max 30827)
usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
cat result/gg/otutab.stat

12. Space cleanup and data submission.

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















































































