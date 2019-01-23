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
python scripts/trim_adapters.py

