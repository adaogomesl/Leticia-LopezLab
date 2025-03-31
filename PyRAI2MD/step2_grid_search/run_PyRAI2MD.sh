#!/bin/sh

#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=0-23:59:00
#SBATCH --job-name=dbh-grid
#SBATCH --partition=short
#SBATCH --mem=30Gb
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm

export PATH=/work/lopez/Python-3.7.4/bin/:$PATH
export LD_LIBRARY_PATH=/work/lopez/Python-3.7.4/lib:$LD_LIBRARY_PATH

export INPUT=input
export WORKDIR=$PWD


cd $WORKDIR
pyrai2md $INPUT
