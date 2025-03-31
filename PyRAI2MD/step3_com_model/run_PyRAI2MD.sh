#!/bin/sh
## script for PyRAI2MD
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:59
#SBATCH --job-name=dbh_comm
#SBATCH --partition=short
#SBATCH --mem=11000mb
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm
 
export PATH=/work/lopez/Python-3.7.4/bin/:$PATH
export LD_LIBRARY_PATH=/work/lopez/Python-3.7.4/lib:$LD_LIBRARY_PATH

export INPUT=input
export WORKDIR=$PWD

cd $WORKDIR
pyrai2md $INPUT

