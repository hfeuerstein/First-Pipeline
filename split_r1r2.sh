#! usr/bin/env bash

# Splits input bam files into r1 and r2 bams and converts them into fastq files 


function split_r1_r2 {
	for file in $@
	do
		echo "Chromosome of file $file" ;
		read;
		samtools view -b -f 0x40 $file > ${REPLY}_r1.bam
		samtools view -b -f 0x80 $file > ${REPLY}_r2.bam
		samtools fastq ${REPLY}_r1.bam > ${REPLY}_r1.fastq
		samtools fastq ${REPLY}_r2.bam > ${REPLY}_r2.fastq
	done
}
