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
	# Choose prefixes for alignment output files:
	read -p "Choose file prefixes:" answer
	mkdir fastqc discarded_seq alignment_results counts
	
	#1) Quality Control with FastQC
	fastqc -o fastqc $@

	#2) Adapter trimming with trimmomatic: 
	trimmomatic PE $1 $2 \
		${answer}_r1_trim ./discarded_seq/${answer}_r1_un.trim \
		${answer}_r2_trim ./discarded_seq/${answer}_r2_un.trim \
		ILLUMINACLIP:./required/NexteraPE-PE.fa:2:40:15

	#3) Alignment with STAR
	STAR --runThreadN 6 --genomeDir ./required/index --readFilesIn *r1_trim *r2_trim --outFileNamePrefix ./alignment_results/${answer}_ --outSAMtype BAM SortedByCoordinate

	#4) Quantification with HTSeq:
	htseq-count --format=bam --order=pos --stranded=no ./alignment_results/*out.bam ./required/*annotation.gtf > ./counts/count_table.tab 
}




