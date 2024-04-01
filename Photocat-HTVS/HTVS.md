# Installing Pyflow 

#### Reference VERDE Paper: https://pubs.acs.org/doi/10.1021/acs.jpclett.9b02577
This Github page contains instructions for installing Pyflow: 
```
https://github.com/northeastern-rc/pyflow-installation-script
```

### Reminders:
- Make sure to have PyFlow.txt in the same directory as install_pyflow.sh.
- If you have not loaded Gaussian and GAMESS on your .bashrc, you must uncomment lines 58 to 66 in the install_pyflow.sh file
- Run pyflow on your scratch

### Testing installation

#### Test 1 - Creating molecules
```
python pymolgen-bench.py bench-MB 'CN(C1=CC=C(N=C2C=C([U])C(N(C)C)=CC2=[S+]3)C3=C1)C'
```

#### Test 2 - Set up Workflow
```
pyflow setup bench_round001 --config_file config-001.json
```

#### Test 3 - Submit Workflow
```
cp bench-MB/* bench_round001/unopt_pdbs/ && cd bench_round001/ && pyflow begin --do_not_track
```

#### Test 4 - Check progress
```
pyflow progress
```
--- 

# Generating Molecules

#### Draw molecule on Chemdraw and select substituent location by using U for spacers
<img width="149" alt="Example Core" src="https://github.com/Kimpton22/Tutorials-And-Guides/assets/100699955/c88389c5-64fc-41dc-9a27-c6f020c07565">

#### Changing spacers and terminals
Below are the spacers and terminal on the shared pymolgen script, if you need to change them, please edit the script

<img width="387" alt="terminal and spacers" src="https://github.com/Kimpton22/Tutorials-And-Guides/assets/100699955/c6599344-0b81-451d-9aa9-5a2715cfcc70">

#### Generating molecule
1. Remember to source your pyflow environment and request resources
2. To generate pdb files, use the following command, replacing "SMILES" with the actual SMILES string.
   ```
   python pymolgen.py 'SMILES'
   ```
--- 

# Creating and submitting workflows

#### 1. Set up workflow: change XXX for the workflow name
   ```
pyflow setup XXX --config_file verde-config.json
   ```

#### 2. Copy molecules generated to unopt_pdbs directory inside the workflow folder 
_Limit of 1000 pdbs per workflow, conformers must be on the same workflow_

#### 3. Go inside workflow directory and submit the following command
   ```
pyflow begin
   ```
_If workflow with same name was submitted before, add --do_not_track flag_
 ```
pyflow begin --do_not_track flag
   ```

#### 4. Check progress
   ```
pyflow progress
   ```

Scheme of performed jobs. Part 1 will be performed when generating molecules.
<img width="1161" alt="workflow" src="https://github.com/Kimpton22/Tutorials-And-Guides/assets/100699955/0fe723f7-a8d0-492c-a831-ea51a9d07731">

If need to change jobs type, edit the config file. Details on the keywords on PyFlow GitHub: https://github.com/kuriba/PyFlow

---
## Clean up workflow
Remove temporary files for failed jobs, for completed jobs them are automatically deleted
```
rm -r workflow*/*/*/failed/*.chk
```
```
rm -r workflow*/*/*/failed/*.rwf
```

## Extract results from Workflow using gather-results-withsp.py
1. Activate pyflow
```
conda activate pyflow
```
2. Run extraction script
```
python3 gather-results-withsp.py name-workflow
```

## Generate plots
### S0 Potentials + ColorMap showing wavelength
CSV details
  1. Need to be named S0.csv
  2. Collumn must have specific names: Oxi_S0,Red_S0,Wavelength

Usage
```
python3 plot-s0-colormap.py S0
```

Resulted plot
![plot_S0](https://github.com/adaogomesl/Leticia-LopezLab/assets/100699955/2a48e439-f486-4e13-8054-7ffd50c43fbd)

### Excited State Potentials + ColorMap showing wavelength
CSV details
  1. Need to be named S1.csv or T1.csv
  2. Collumn must have specific names: Oxi_S1,Red_S1
     
_in case of T1, replace S1 for T1_

Usage
```
python3 plot-scatter.py S1
```
Resulted plot
![plot_S1](https://github.com/adaogomesl/Leticia-LopezLab/assets/100699955/91ecc5c8-9702-404a-8524-cc4ce6e13bf2)

