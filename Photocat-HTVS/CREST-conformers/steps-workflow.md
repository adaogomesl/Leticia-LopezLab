1. Create folders for each core
for i in {10..21}; do mkdir $i; done

2. Copy necessary files: pymolgen-id-CREST.py and text file with SMILES from step1, step2 and step3

3. Generate SMILES, InchKey, ID and set CREST sbatch
	1. Request resources
 	2. Activate pyflow
	3. Copy ```pymolgen-id-CREST.py```, path: ```/work/lopez/share_from_Leticia/forSteven/scripts/pymolgen-id-CREST.py ```
 	4. Submit script: ```python pymolgen-id-CREST.py molecule_name smiles.txt```
  	5. Submit CREST search: ```sbatch runall.py```

	

4.Check if all conf_search folders generated an XYZ: ```for dir in  */; do [[ $(ls "$dir"*.xyz 2> /dev/null) ]] || echo "$dir has no .xyz files"; done```


5. Copy xyz files: ```cp ../*/*.xyz .```


4. Convert xyz to pdb
	a. activate pyflow
        b. Convert .xyz to .pdb: ```obabel *.xyz -O *.pdb```

6. Create folder just with pdb ```mv ../xyz/*pdb .```

