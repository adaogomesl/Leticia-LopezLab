# 1. Single-point HF - Generate initial orbitals
## References
GATEWAY: https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/gateway.html
SEWARD: https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/seward.html
SCF: https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/scf.html


## Files necessary
```.inp```,```.xyz``` (from dft opt) and ```.sh``` 

#### .inp example
```
>> EXPORT MOLCAS_MAXITER=10000
&GATEWAY
 coord=dbhme-T-sp-HF.xyz
 basis=ano-s-vdzp
 Group=c1
 RICD

&SEWARD
doanalytic

&SCF
```

#### .xyz example
```
18
dbh
 C                  1.49605900   -0.55258700    0.78740200
 C                  1.49540400   -0.56520300   -0.77950200
 C                  0.05621200   -0.05471400   -1.11121600
 C                  0.05686800   -0.03605400    1.11192000
 C                 -0.79550600   -0.71081500    0.00565100
 N                  0.09210500    1.40950500   -0.64606700
 N                  0.09271600    1.41816200    0.62448100
 H                 -0.63772200   -1.79563800    0.01546000
 H                 -0.26299900   -0.09421300    2.14897300
 H                 -0.26448300   -0.12705800   -2.14708400
 H                  2.25600200    0.09028200   -1.21102000
 H                  2.25642900    0.11074400    1.20732500
 H                  1.64624000   -1.55223700    1.20659400
 H                  1.64467900   -1.57196100   -1.18189700
 C                 -2.29233100   -0.39597800    0.00323900
```

#### .sh example
```
#!/bin/sh
## script for OpenMalCas
## $INPUT and $WORKDIR do not belong to OpenMolCas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=23:00:00
#SBATCH --job-name=dbh-me-T-sp-HF
#SBATCH --partition=short
#SBATCH --mem=6Gb
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm

export INPUT=dbhme-T-sp-HF
export WORKDIR=$PWD

export MOLCAS_NPROCS=$SLURM_NTASKS
export MOLCAS_MEM=3200
export MOLCAS_PRINT=3
export MOLCAS_PROJECT=$INPUT
export MOLCAS_WORKDIR=$WORKDIR
export OMP_NUM_THREADS=1
export MOLCAS=/work/lopez/Molcas

module load python/3.7.1

cd $WORKDIR
$MOLCAS/bin/pymolcas -f $INPUT.inp -b 1
rm -r $MOLCAS_PROJECT
```
## Analysis
Check .log and .status should have ```Happy landing!```

--- 
# 2. Single-point CASSCF - Generate CAS

## References
RASSCF: [https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/gateway.html](https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/rassi.html)

## Files necessary
```.inp```,```.xyz```,```.sh``` and ```.StrOrb``` (ScfOrb  from previous step)

#### .inp example
```
&GATEWAY
 coord=$Project.xyz
 basis=ano-s-vdzp
 group=C1
 RICD

&SEWARD
 doanalytic

&RASSCF
 FileOrb=$Project.StrOrb
 Spin=1
 Charge=0
 NActel=8
 Ras2=9
 CIRoot=8 8 1
```

## Analysis
1. Check .log and .status should have ```Happy landing!```
2. Open .rasscf.molden on ChemCraft to vizualize orbitals in the active space (Energy = 0)

--- 

