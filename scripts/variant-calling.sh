#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 4
#SBATCH --mem=50G


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

#fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
#fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/


# not filtering required since the quality looks good. Checking for adapter content and low quality bases. Low quality bases were not going to be ommitted anyways since the aligned we are using uses soft clipping and we don not want to shorten the length of our reads which might result in incorrect mapping. 



# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}


# BWA alignment (resulting SAM or BAM will be missing the read group section so we are providing it here)
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam


