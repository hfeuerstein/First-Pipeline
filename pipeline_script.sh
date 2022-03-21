#!/usr/bin/env bash

# Basic bash script for pre-processing fastq files with quality control (FastQC), adapter-trimming (Trimmomatic), alignment (STAR) and quantification (featurecount in R) 
# Input: Fastq file pair: #1 = read1; #2 = read 2
# Trimmomatic: requires NexteraPE-PE.fa in the directory => path: ~/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa  

function preprocess {
	# Makes for each file input (r1 and r2) a new directory and stores the QC files in it 
	for file in $@
	do
		mkdir fastqc_"$file"
		fastqc -o fastqc_"$file" $file
	done
	# trimms read1 and read2 fastq files for NexteraPE adapters and stores the suriving sequences in the current directory and the discarding reads in dir:discarded sequences
	mkdir discarded_sequences 
	trimmomatic PE $1 $2 \
		"$1"_trim ./discarded_sequences/"$1"_un.trim \
		"$2"_trim ./discarded_sequences/"$2"_un.trim \
		ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

	# Alignment: requires indexed genome file "index" in current directory (create index with STAR) 
	mkdir alignment_results
	echo: "Name prefix for outputfiles_"
	read;
	STAR --runThreadN 6 --genomeDir index --readFilesIn $1 $2 --outFileNamePrefix ./alignment_results/${REPLY} --outSAMtype BAM SortedByCoordinate
	samtools index ./alignment_results/*.bam

}
