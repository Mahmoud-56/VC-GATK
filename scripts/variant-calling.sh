#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --error=job_error.log

#Script for germline variant calling in a human WGS paired end reads 2 X 100bp


if false
then
# download data
wget -P /projects/bgmp/malm/bioinfo/VC-GATK/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P /projects/bgmp/malm/bioinfo/VC-GATK/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


## Prep files (TO BE GENERATED ONLY ONCE) 

# download reference files
wget -P /projects/bgmp/malm/bioinfo/VC-GATK/ref  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip /projects/bgmp/malm/bioinfo/VC-GATK/ref/hg38.fa.gz

# index ref - .fai file before running haplotype caller
samtools faidx /projects/bgmp/malm/bioinfo/VC-GATK/ref/hg38.fa


# ref dict - .dict file before running haplotype caller
gatk CreateSequenceDictionary R=/projects/bgmp/malm/bioinfo/VC-GATK/ref/hg38.fa O=/projects/bgmp/malm/bioinfo/VC-GATK/ref/hg38.dict


# download known sites files for BQSR from GATK resource bundle
wget -P /projects/bgmp/malm/bioinfo/VC-GATK/ref/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /projects/bgmp/malm/bioinfo/VC-GATK/ref https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx


fi

# ref dict - .dict file before running haplotype caller

gatk CreateSequenceDictionary R=/projects/bgmp/malm/bioinfo/VC-GATK/ref/hg38.fa O=/projects/bgmp/malm/bioinfo/VC-GATK/ref/hg38.dict
## VARIANT CALLING STEPS


# directories
ref="/projects/bgmp/malm/bioinfo/VC-GATK/ref/hg38.fa"
known_sites="/projects/bgmp/malm/bioinfo/VC-GATK/ref/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/projects/bgmp/malm/bioinfo/VC-GATK/aligned_reads"
reads="/projects/bgmp/malm/bioinfo/VC-GATK/reads"
results="/projects/bgmp/malm/bioinfo/VC-GATK/results"
data="/projects/bgmp/malm/bioinfo/VC-GATK/data"


# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

#echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/


# not filtering required since the quality looks good. Checking for adapter content and low quality bases. Low quality bases were not going to be ommitted anyways since the aligned we are using uses soft clipping and we don not want to shorten the length of our reads which might result in incorrect mapping. 



# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}


# BWA alignment (resulting SAM or BAM will be missing the read group section so we are providing it here)
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam



# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} -bqsr ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 



# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics -R /projects/bgmp/malm/bioinfo/VC-GATK/ref/hg38.fa -I /projects/bgmp/malm/bioinfo/VC-GATK/aligned_reads/SRR062634_sorted_dedup_bqsr_reads.bam -O /projects/bgmp/malm/bioinfo/VC-GATK/aligned_reads/alignment_metrics.txt

gatk CollectAlignmentSummaryMetrics -R ${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${aligned_reads}/alignment_metrics.txt

gatk CollectInsertSizeMetrics --INPUT ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam --OUTPUT ${aligned_reads}/insert_size_metrics.txt --Histogram_FILE=${aligned_reads}/insert_size_histogram.pdf

gatk CollectInsertSizeMetrics --INPUT /projects/bgmp/malm/bioinfo/VC-GATK/aligned_reads/SRR062634_sorted_dedup_bqsr_reads.bam --OUTPUT /projects/bgmp/malm/bioinfo/VC-GATK/aligned_reads/insert_size_metrics.txt --Histogram_FILE /projects/bgmp/malm/bioinfo/VC-GATK/aligned_reads/insert_size_histogram.pdf

# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
