#!/bin/sh

#Partition
#SBATCH -p main

#The name of the job
#SBATCH -J merger

#The job requires 1 compute nodes
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a half hour
#SBATCH -t 1-00:00:00

bash cat_2.txt
