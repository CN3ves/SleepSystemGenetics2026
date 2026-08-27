#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --array 1-226

#SBATCH --account=pfranken_bxd_atac

#SBATCH --mem-per-cpu 10 # Megabytes
#SBATCH --time 02:00:00 

# Get sample name from list
file=$1
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p $file)
echo "Merging $SAMPLE on `hostname` (`date`)"

# Get FASTQ files for sample
fastqs=$(ls rawdata/atac/* | grep "/"$SAMPLE"_" )
echo "Merging files:"
echo $fastqs
# Merge reads
zcat $fastqs | gzip -nc  > rawdata/atac/fastq/$SAMPLE.fastq.gz 

echo "DONE $SAMPLE"
