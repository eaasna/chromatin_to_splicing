#!/bin/sh

#The job should run on the main partition
#SBATCH -p main

#The name of the job
#SBATCH -J pca

#The job requires 1 compute nodes
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a half hour
#SBATCH -t 1-00:00:00

/gpfs/hpchome/a72094/software/bin/QTLtools pca --bed /gpfs/rocket/home/a72094/projects/chromatin_to_splicing/normalized/counts_matrix_cqn_genotypes.bed.gz --scale --center --out /gpfs/rocket/home/a72094/projects/chromatin_to_splicing/normalized/counts_matrix_cqn_40

