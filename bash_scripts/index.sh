#!/bin/sh

#The job should run on the main partition
#SBATCH -p main

#The name of the job
#SBATCH -J sort

#The job requires 20 tasks. Specifying only tasks and not nodes defaults to one task being run per available CPU
#SBATCH --ntasks=20

#The maximum walltime of the job
#SBATCH -t 1-00:00:00

module load samtools-1.6
zcat /gpfs/rocket/home/a72094/projects/chromatin_to_splicing/QTL/cqn_nominal_full.txt.gz | awk -v OFS='\t' '{{$1=$1; print $0}}' | srun sort --parallel=20 -k9,9n -k10,10n -k11,11n | bgzip > /gpfs/rocket/home/a72094/projects/chromatin_to_splicing/tabix/cqn_nominal.sorted.txt.gz