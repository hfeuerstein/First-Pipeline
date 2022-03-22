#! usr/bin/env bash 

# Input bam is converted into two fastq files (r1 and r2) and saves output into folder bam_files or fastq_files, respectivels 
# requires "samtools" 


function bam2fastq {
	mkdir fastq_files
	mkdir bam_files
	for bam_file in $@
	do
		samtools index $bam_file

		read -p "Shall only a particular chromsome be retranslated into fastq? [y/n]" answer
		if [[ $answer = y ]]
		then
			read -p "Name chromosome [chrxyz]" chromosome
			read -p "Name of the sample:" sample

			# Creates subsetted bam file of chromosome in folder "bam_files" 
			samtools view -b $bam_file $chromosome > ./bam_files/$chromosome.bam

			# Splits chr bam_file into read 1 and read 2 according to flags (0x40 = r1, 0x80 = r2)
			samtools view -b -f 0x40 ./bam_files/${chromosome}.bam > ./bam_files/${sample}_${chromosome}_r1.bam
			samtools view -b -f 0x80 ./bam_files/${chromosome}.bam > ./bam_files/${sample}_${chromosome}_r2.bam

			# Converts r1/2 bam files into fastq files and saves them in folder "fastq_files" 
			samtools fastq ./bam_files/${sample}_${chromosome}_r1.bam > ./fastq_files/${sample}_${chromosome}_r1.fastq
			samtools fastq ./bam_files/${sample}_${chromosome}_r2.bam > ./fastq_files/${sample}_${chromosome}_r2.fastq
		else
		    # If the entire bam file shall be retranslated into fastq files 
			read -p "Name of the sample" sample
			samtools view -b -f 0x40 $bam_file > ./bam_files/${sample}_r1.bam
			samtools view -b -f 0x80 $bam_file > ./bam_files/${sample}_r2.bam
			samtools fastq ./bam_files/${sample}_r1.bam > ./fastq_files/${sample}_r1.fastq
            samtools fastq ./bam_files/${sample}_r2.bam > ./fastq_files/${sample}_r2.fastq
		fi
	done
}
