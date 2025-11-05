# Step 1: Geodesic Interpolation
This step will generate the initial training data set for us to use with PyRAI2MD to conduct machine learning - non-adiabatic molecular dynamics. 

NOTE: You need to do the following steps TWICE, one is for the reactant -> a minimum energy conical intersection (MECI), and then another from MECI -> product. Each of these will have their folders. 

Geometry info:
MECI - Partially Inverted
Product

1. `react_meci` -> contains the reactant -> MECI geometries
2. `meci_product` -> contains the MECI -> product geometries

Lets start with the first_part folder (i.e. reactant -> MECI): 

# Step 1: Interpolation
Files that are needed: 
1. interpolation.sh - shell script that will run the geodesic interpolation.
2. point1-point2.xyz - an xyz file that contains 2 geometries; (1) reactant geometry (i.e. point 1) and (2) the product geometry (i.e. point 2). These are the two points where the geodesic interpolation will be done for.

How to run:
``` 
bash interpolation.sh {xyz file containing both geometries} {number of steps}
```

Example xyz file: 
```
  13
 Reactant
 N     1.63675933    -0.62894955     0.01595707
 N     1.63676838     0.62893842     0.01583925
 B    -0.74947645     0.00000540    -0.02846574
 C    -2.32215207     0.00000076     0.01015706
 H    -2.77251644    -0.89032567    -0.43410182
 H    -2.60923154    -0.00021811     1.07166147
 H    -2.77254750     0.89047523    -0.43376812
 C     0.21778192    -1.23896749    -0.02003786
 H     0.16055291    -1.87661794    -0.90737989
 H     0.12548282    -1.90431280     0.84328345
 C     0.21779631     1.23896583    -0.02015607
 H     0.12550241     1.90442613     0.84307539
 H     0.16056991     1.87649977    -0.90758418
  13
 MECI
 N        1.039489830     -0.373762820      0.524923480
 N        1.635186900      0.646109510     -0.105555570
 B       -0.540856090     -0.037479010     -0.174452830
 C       -1.973082260     -0.111198330      0.496792310
 H       -2.739015160     -0.110176770     -0.289608730
 H       -2.122361370     -1.021791180      1.083564220
 H       -2.179193310      0.747141350      1.142434880
 C        0.374969060     -1.334625360     -0.478862250
 H        0.977044250     -1.504652610     -1.367792250
 H        0.087323390     -2.245091190      0.034563530
 C        0.284862650      1.283476440     -0.609281680 
 H        0.082785900      2.197595280     -0.052745160
 H        0.454405890      1.524236820     -1.661439060
```

Example Command:
```
bash interpolation.sh reactant_meci.xyz 20
```

Outputs: 
1. output.xyz - a xyz file that contains all the interpolated steps.

For the example above you would expect the output.xyz to have 20 geometries of equidistant steps between the reactant and MECI.

# Step 2: Separation
Files that are needed: 
1. separate.py - python script that will separate the xyz files from a combined xyz (i.e. output.xyz, wigner.300.xyz).

How to run: 
```
python3 separate.py {file_name.xyz} {output name}
```
Example: 
```
python3 separate.py output.xyz bd-td
```

Outputs: 
1. {output_name}[1-n].xyz - numbered xyz files containing the individual xyz from the interpolations, where n is the number of interpolated geometries.

For the example above you would get a list of xyz files title bd-td-#.xyz, where # would be 1-20.

# Step 3: Alignment
Files that are needed:
1. overlay.py - python script that will align the geometries with a reference structure.
2. list.txt - text file that contains the filenames of the reference and interpolated geometries; you can generate the list with the following bash line:
   ```
   for ((a=1;a<={# of interpolated geometries};a++)); do echo {filename}-$a.xyz>> list.txt; done
   ```
   Example of list.txt file, the first entry needs to be the reference structure:
   ```
   bd-td.xyz
   bd-td-1.xyz
   bd-td-2.xyz
   bd-td-3.xyz
   bd-td-4.xyz
   bd-td-5.xyz
   bd-td-6.xyz
   bd-td-7.xyz
   bd-td-8.xyz
   ```
4. reference.xyz - a reference xyz file with the original orientation (i.e. S0 optimized structure).

How to run:
```
python3 overlay.py list.txt
```

Outputs:
1. {filename}[1-n]-overlay.xyz - numbered xyz files aligned with the original xyz file; where n is the number of the interpolated geometries.

For the example above you would expect an overlaid version of the 20 interpolated geometries.

NOTE: Make a directory for the overlaid geometries (i.e. `overlaid_geos`) and ove all the overlay geometries into this folder and use the following bash command to delete the '-overlay' from the file name. 
```
for i in *; do mv "$i" "${i%%-overlay.xyz}.xyz"; done
```

# Step 4: Wigner Sampling
1. seperate.py - python scripts that separate the Wigner sampled geometries.
2. wigner-{filename}.300.0.xyz - an xyz file that contains all the geometries generated from the Wigner sampling.

NOTE: Change the following line in the `separate.py` script:
`ngeom = # of wigner sample geometries that you want`

How to run:
```
python3 seperate.py wigner-{filename}-300.0.xyz {output name}
```

Example: 
```
python3 seperate.py wigner-bd-td-300.0.xyz bd-td
```

Outputs: 
1. {output name}-[1-n].xyz - numbered xyz files for the individual Wigner sampled geometries.

For the example above you would expect to have the first n geometries of the Wigner sampling of the reactant.

# Step 5: Organization
This step is to organize the files to be able to set up single calculations for all the interpolated structures. 

Once you have complete steps 1-4 for both ```first_part``` (i.e. reactant -> MECI) and ```second_part``` (i.e. MECI -> product). 

1. Make directory title ```interpolation```
2. `cd interpolation`
3. Make directories for each part (i.e. `first_part` and `second_part`)
4. Copy the overlaid geometries from the `overlaid_geos` folders into their respective folders
5. Move the directory that contains the Wigner sample geometries into the `interpolation` folder

# Step 6: Set-up Calculations
Files needed: 
1. Gen-wigner-interp.py - python script that generates the MOLCAS single-point calculations.
2. {filename}.inp - MOLCAS input file.
3. {filename}.StrOrb - MOLCAS starting orbitals, usually from the S0 optimized calculations.

How to run: 
```
python3 Gen-wigner-interp.py {filename}
```
Outputs: 
1. {filename}[1-n]/bd-td-[1-m] - folders containing the single point calculations, where n is the number of Wigner sample geometries and m is the number of interpolated steps.
2. runall.sh - shell script that will run all the calculations; usually split this due to hitting the limit of the job submission on the Discovery Cluster.

# Step 6: Data Extraction
Files that are needed: 
1. check.sh.
2. input 
3. training_data_tool.py - collects the single point energies into a json file.
4. data_shuffle_tool.py - shuffles the json file, so that it NN can train better.

How to run: 
The following line will let you know of any single-point calculations that failed.
```
bash check.sh
```
The following lines collect the data. 

```
python3 training_data_tool.py -i input -n 20
```

Where the -i flags an input file and -n flags the number of computers that will be used to collect the data. 

Example input: 

```
INSERT INPUT EXAMPLE
```

The following will shuffle it

```
python3 data_shuffle_tool.py data.json
```



