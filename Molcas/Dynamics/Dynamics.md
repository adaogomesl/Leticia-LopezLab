## Submission

### Necessary files
*From Optimization:*
1. .freq.molden
2. .xyz
3. StrOrb

*Script/File*
1. Gen-FSH.py
2. Control (number of molecules = nodes x cores)

*Input file*
```
&GATEWAY
 coord=$MOLCAS_PROJECT.xyz
 basis=ano-s-vdzp 
 Group=c1
 RICD

>> FOREACH ITER in (1 .. 2000)

&SEWARD
doanalytic

&RASSCF
 Fileorb=$MOLCAS_PROJECT.StrOrb
 Spin=1
 Nactel=10 0 0 
 Charge=0
 Ras1=0
 Ras2=9
 Ras3=0
 ITERATIONS=200,100
 CIRoot=8 8 1
 MDRLXR=2

&Surfacehop
 tully
 decoherence = 0.1

&ALASKA

&Dynamix
 velver
 dt = 20
 velo = 3
 ther= 2
 temp = 300

>> End do
```

### Running dynamics
1. Request resources
```
srun -N 1 --exclusive --partition=short --time=23:59:59 --pty /bin/bash
```
2. Load python
```
module load python/3.7.1
```
3. Edit the name on the control file
```
      input      dbh-unsub.freq.molden
      temp       300
      method     wigner
      partition  lopez
      time       30-00:00:00
      memory     2000
      nodes      20
      cores      25
      jobs       25
      index      1
```
5. Generate input file and wigner samples - It will generate the wigner.xyz
```
python Gen-FSSH.py -x control
```
5. To run dynamics, add in the first line of the runall.sh:  
```
 #!/bin/sh
```
6. Submit runall

## Analysis
1. Request resources
```
srun -N 1 --exclusive --partition=short --time=23:59:59 --pty /bin/bash
```
2. Load python
```
module load python/3.7.1
```









