## Submission
### Model scripts
```
/work/lopez/share_from_Leticia/model-scripts/dynamics
```


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
### Diagnostic - Check final timestep
```
python3 HOP-FSSH.py -x 01-diagnostic
```
1. Create txt file with timesteps printed on the screen
2. Make a histogram to get the distribution of the timestep
   
### Determine Final State
```
python3 HOP-FSSH.py -x 02-
```
1. Classify final products using classify_product.py
```
python3 classify_product.py Prod.S0.xyz
```
2. Use results on ```.csv``` to classify products


### State Population Analysis
```
python3 HOP-FSSH.py -x 03-
```
### Spaghetti Plots
```
python3 HOP-FSSH.py -x 04-
```








