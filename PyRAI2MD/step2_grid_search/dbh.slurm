#!/bin/sh
## script for PyRAI2MD
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=23:59:59
#SBATCH --job-name=dbh-search
#SBATCH --partition=lopez
#SBATCH --mem=11000mb
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm

