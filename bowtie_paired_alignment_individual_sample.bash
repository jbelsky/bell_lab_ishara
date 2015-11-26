#!/bin/bash

# Set the dataset types
# dataset=("IA2_SWI-SNF" "IA3_RSC" "IA4_A-_B-_ISW1a" "IA5_ISW1a_Mcm" "IA6_SWI-SNF_Mcm" "IA7_RSC_Mcm")
dataset_id=($(seq 4 1 10))
dataset_name=("7568_ISW1b" "7569_ISW2" "7570_INO80" "7571_SWI-SNF" "7572_RSC" "7573_CHD1" "7574_Nap1")

# for(( d=0; d<${#dataset_id[@]}; d++ ))
for(( d=4; d<6; d++ ))
do

	# Set the output directory
	work_dir="/data/collaborator_data/aligned_experiments/CD${dataset_id[$d]}/"

	# Set the file_directory
	fastq_file_R1="/data/collaborator_data/sequencing_archive/CD${dataset_id[$d]}/CD${dataset_id[$d]}_${dataset_name[$d]}_R1.fastq"
	fastq_file_R2="/data/collaborator_data/sequencing_archive/CD${dataset_id[$d]}/CD${dataset_id[$d]}_${dataset_name[$d]}_R2.fastq"

	# Set the output name
	output="CD${dataset_id[$d]}"

	# Set the alignment genome
	align_genome="yeast/plasmid/bowtie_index/pUC19_ARS1_WT_updated"

	# Create the work dir directory if it does not exist
	if [ ! -d $work_dir ]
	then
		echo "There is no work_dir directory here!"
		mkdir -v $work_dir
	fi

	# Enter into the work directory
	cd $work_dir

	# Create the log file directory if it does not exist
	if [ ! -d "log_files" ]
	then
		echo "There is no log files directory here!"
		mkdir -v "log_files"
	fi

	# Create the log file
	log_file=log_files/${output}_bowtie_log.txt
	temp_log_file=log_files/${output}_bowtie_log_temp.txt

	# Copy the fastq filees to the working directory
	cp -v -t $work_dir ${fastq_file_R1}.bz2 ${fastq_file_R2}.bz2

	# Get the base name
	fastq_file_R1=$(basename $fastq_file_R1)
	fastq_file_R2=$(basename $fastq_file_R2)

	# Unzip the fastq files
	pbzip2 -v -d -p6 ${fastq_file_R1}.bz2 ${fastq_file_R2}.bz2

	# Run bowtie
	echo -e "Aligning \n\t$fastq_file_R1 and \n\t$fastq_file_R2..."

	bowtie -p 6 -n 2 -l 20 --phred33-quals --best --strata -M 1 -k 1 -X 1000 -S -y --chunkmbs 256 \
	/data/illumina_pipeline/genomes/${align_genome} \
	-1 $fastq_file_R1 -2 $fastq_file_R2 ${output}.sam >> $temp_log_file 2>&1

	# Filter only the bowtie output of the log file
	grep -v "^Warning" $temp_log_file >> $log_file
	rm $temp_log_file

	# Convert to BAM only keeping the positive aligned reads
	echo -e "\tConverting from SAM to BAM..."
	samtools view -@ 4 -bhS -f 32 -o ${output}.bam ${output}.sam

	# Sort the bam file
	echo -e "\tSorting..."
	samtools sort -@ 4 -l 9 -o ${output}.bam -T ${output}_temp ${output}.bam

	# Index the final bam file
	samtools index ${output}.bam

	# Remove the fastq and SAM files
	rm -v ${output}.sam $fastq_file_R1 $fastq_file_R2

done
