#!/bin/bash

#The job should run on the testing partition
#SBATCH -p testing

#The name of the job is test_job
#SBATCH -J nextera

#The job requires 4 compute nodes
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a half hour
#SBATCH -t 00:30:00

#Actual commands
module load python
python scripts/extractNexteraBarcode.py --fastq /gpfs/hpchome/evelin95/fastq/babk_A_ATAC.1.fastq.gz --type read1 > /gpfs/hpchome/evelin95/trimmed/babk_A_ATAC.1.barcode.txt && python scripts/extractNexteraBarcode.py --fastq /gpfs/hpchome/evelin95/fastq/babk_A_ATAC.2.fastq.gz --type read2 > /gpfs/hpchome/evelin95/trimmed/babk_A_ATAC.2.barcode.txt
 
