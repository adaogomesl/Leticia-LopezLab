# Benchmarking Density Functional Methods for Accurate Prediction of Ground and Excited State Reduction and Oxidation Potentials of Organic Photoredox Catalysts
This page contains all the files necessary to run an automated workflow to predict ground and excited state potentials. Also, contains a script to extract all properties.

All optimized structures output files are available at 
```
10.6084/m9.figshare.25234249.
```


# Installing Pyflow 

This Github page contains instructions for installing Pyflow: 
```
https://github.com/kuriba/PyFlow
```

# Best Performance model chemistries
In the mentioned manuscript, we evaluated 147 model chemistries to predict the ground- and excited-state redox potentials for cyanoarenes, benzophenones, xanthene, and acridinium catalysts. Top-performing model chemistries were identified for predicting excited state oxidation and reduction potentials; N12-SX/6-311+G(d,p) and PBE0-D3BJ/6-311+G(d,p) are best for S1 potentials, and BHandH/6-31+G(d,p) and ωB97X/6-311+G(d,p) are best for T1 potentials.

# Workflow
<img width="600" alt="image" src="https://github.com/adaogomesl/Leticia-LopezLab/assets/100699955/2e146a2a-2a06-4186-bc45-53b71bf93e62">


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

# Extract Potentials
#### 1. Copy extraction script to the same directory of the workflow directory - gather-results.py
```
cp /work/lopez/share_from_Leticia/verde-pyflow/gather-results.py .
```

#### 2. Extract the results - Replace workflow_name by the workflow directory name
_Reminders: Request resources and have PyFlow environment sourced_
```
python gather-results.py workflow_name
```

A CSV file will be generated with computed properties. Example of CSV is below:

<img width="1018" alt="Screenshot 2024-01-17 at 12 51 52 PM" src="https://github.com/Kimpton22/Tutorials-And-Guides/assets/100699955/a4b7a7ef-856f-46a4-b440-de139220a057">
