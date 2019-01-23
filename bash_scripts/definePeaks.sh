#!/bin/sh

#Partition
#SBATCH -p main

#The name of the job
#SBATCH -J definePeaks

#The job requires 1 compute node
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a hour
#SBATCH -t 1-00:00:00

#Actual commands
module load R-3.3.0
Rscript scripts/test.R cqn_permutations_10000
Rscript scripts/test.R cqn_permutations_100000
Rscript scripts/test.R fpkm_permutations_10000
Rscript scripts/test.R fpkm_permutations_100000