#!/bin/bash

#the job should run on the main partition
#SBATCH -p main

#The name of the job is test_job
#SBATCH -J tabix

#The job requires 4 compute nodes
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a half hour
#SBATCH -t 00:30:00

#Actual commands
module load samtools-1.6
bgzip /gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_contained_enriched.sorted.txt
tabix -s9 -b10 -e11 -f /gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_permutations_contained_enriched.sorted.txt.gz



