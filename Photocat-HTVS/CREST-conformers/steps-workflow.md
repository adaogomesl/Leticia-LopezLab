1. Create folders
for i in {2..9}; do mkdir 0$i; done
for i in {10..21}; do mkdir $i; done

2. Submit by using sbatch with commands, example
cd /scratch/adaogomes.l/forsteven/workflow-Apr24/015_step1/01 && cp ../workflow-step1.py . && sbatch workflow-step1.py LZMYYVFMLBATPX-UHFFFAOYSA-N 'O=C(C(C(OC)=O)=C1C(OC)=O)C(C(OC)=O)=C(C(OC)=O)C1=O'

3.Copy xyz files
cp ../0*/*.xyz .
cp ../1*/*.xyz .
cp ../2*/*.xyz .
cp ../3*/*.xyz .

4. Convert log to pdb
	a. activate pyflow
        b. Convert .xyz to .pdb: obabel *.xyz -O *.pdb

6. Create folder just with pdb

