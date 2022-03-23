#!/usr/bin/env bash

###############################################
# Function "preprocess": basic bash script for preprocessing Fastq-files incl: 
# 	- quality control (FastQC)
#	- adapter trimming (trimmomatic) for NexteraPE-PE
#	- alignment (STAR)
#	- quantification (HTSeq)
#
# REQUIREMENTS:
#  a) packages: fastqc, trimmomatic, star, htseq 
#  b) files to be stored in a file "required" in current workign directory: 
#		- NexteraPE-PE.fa: avaibale via  ~/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa
#		- GTF file for annotation: with name *annotation.gtf
#		- folder "index" with indexed genome/subset (from STAR)
#
# INPUT: two fastq files (r1 and r2)
##############################################

function preprocess {
	# Extracts filename from path and uses the filename until _r1* or _r2* for output file names 
	NAME=$(basename "$1")
	PREFIX=${NAME%%_r*}

	mkdir fastqc trimmed_fastq discarded_fastq alignment_results counts
	
	#1) Quality Control with FastQC
	fastqc -o fastqc $@

	#2) Adapter trimming with trimmomatic: 
	trimmomatic PE $1 $2 \
		./trimmed_fastq/${PREFIX}_r1_trim.fastq ./discarded_fastq/${PREFIX}_r1_un.trim.fastq \
		./trimmed_fastq/${PREFIX}_r2_trim.fastq ./discarded_fastq/${PREFIX}_r2_un.trim.fastq \
		ILLUMINACLIP:./required/NexteraPE-PE.fa:2:40:15

	#3) Alignment with STAR
	STAR --runThreadN 6 --genomeDir ./required/index --readFilesIn ./trimmed_fastq/*r1_trim.fastq ./trimmed_fastq/*r2_trim.fastq --outFileNamePrefix ./alignment_results/${PREFIX}_ --outSAMtype BAM SortedByCoordinate

	#4) Quantification with HTSeq:
	htseq-count --format=bam --order=pos --stranded=no ./alignment_results/*out.bam ./required/*annotation.gtf > ./counts/count_table.tab 
}