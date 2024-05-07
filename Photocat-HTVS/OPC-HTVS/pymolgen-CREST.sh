#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=23:00:00
#SBATCH --job-name=pymolgen-10conf
#SBATCH --partition=short
#SBATCH --mem=10Gb
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm


python pymolgen-CREST.py A1	'CC(C1=CC=C2C([N+](CC3=CC=C(C(C)(C)C)C=C3)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
python pymolgen-CREST.py A2	'CC(C1=CC=C2C([N+](CC3=CC=C(OC)C=C3)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
python pymolgen-CREST.py A3	'CC(C1=CC=C2C([N+](C3CCCCC3)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
python pymolgen-CREST.py A4	'CC(C1=CC=C2C([N+](CC(C)(C)C)=C(C=C(C(C)(C)C)C=C3)C3=C2C4=C(C)C=C(C)C=C4C)=C1)(C)C'
python pymolgen-CREST.py A5	'CC(C1=CC=C2C([N+](C3CCCCCC3)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
python pymolgen-CREST.py A6	'CC(C1=CC=C2C([N+](C3=CC(OC)=C(OC)C(OC)=C3)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
python pymolgen-CREST.py A7	'CC(C1=CC=C2C([N+](C3=NC(C)=C(C)S3)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
python pymolgen-CREST.py A8	'CC(C1=CC=C2C([N+](C3=NC(C=CC=C4)=C4O3)=C(C=C(C(C)(C)C)C=C5)C5=C2C6=C(C)C=C(C)C=C6C)=C1)(C)C'
python pymolgen-CREST.py A9	'CC(C1=CC=C2C([N+](C3=C(C)C=CC=C3OC)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
python pymolgen-CREST.py A10 'CC(C1=CC=C2C([N+](C3=CC=C(N=CC=N4)C4=C3)=C(C=C(C(C)(C)C)C=C5)C5=C2C6=C(C)C=C(C)C=C6C)=C1)(C)C'
python pymolgen-CREST.py A11 'CC(C1=CC=C2C([N+](C3=CC=C(C(C)(C)C)C=C3)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
python pymolgen-CREST.py A12 'CC(C1=CC=C2C([N+](C3=CC=C(C(F)(F)F)C=C3)=C(C=C(C(C)(C)C)C=C4)C4=C2C5=C(C)C=C(C)C=C5C)=C1)(C)C'
