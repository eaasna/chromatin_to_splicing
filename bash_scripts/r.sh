#!/bin/sh

#The job should run on the main partition
#SBATCH -p main

#The name of the job is rsquared
#SBATCH -J rsquared

#The job requires 1 compute node
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a half hour
#SBATCH -t 1-00:00:00

module load R-3.3.0
module load samtools-1.3
Rscript scripts/Rsquared.R
