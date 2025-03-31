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

export INPUT=input.json
export WORKDIR=/scratch/adaogomes.l/DBH/unsub/ML-NAMD/step2_grid_search/grid-search/NN-dbh-3_300_128_1e-08_1e-08_0.005

cd $WORKDIR
pyrai2md $INPUT
