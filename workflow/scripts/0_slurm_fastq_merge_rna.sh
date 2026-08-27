#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --array 1-32

#SBATCH --account=pfranken_bxd_atac

#SBATCH --mem-per-cpu 10 # Megabytes
#SBATCH --time 02:00:00 

# Get sample name from list
file=$1
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $file)
echo "Merging $SAMPLE on `hostname` (`date`)"

# Get FASTQ files for sample
fastqs=$(ls rawdata/nrf1/* | grep "/"$SAMPLE"_" )
echo "Merging files:"
echo $fastqs
# Merge reads
zcat $fastqs | gzip -nc  > rawdata/nrf1/fastq/$SAMPLE.fastq.gz 

echo "DONE $SAMPLE"


#rename files for publication
#for file in $(ls rawdata/nrf1/fastq/NRF*.fastq.gz);
#do 
#	f=$(basename $file | sed 's/.fastq.gz//g')
#	n=$(cat rawdata/nrf1/nrf_metadata.csv | grep $f | cut -d',' -f2 | sed 's/"//g')
#
#	echo cp $file rawdata/nrf1/fastq/$n.fastq.gz
#	cp $file rawdata/nrf1/fastq/$n.fastq.gz
#
#	echo $(zcat $file | wc -l)
#	echo $(zcat rawdata/nrf1/fastq/$n.fastq.gz | wc -l)
#done