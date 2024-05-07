import os
import sys
import time
import shutil
from datetime import timedelta
import pandas as pd
import openbabel
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import rdBase
from itertools import chain
rdBase.DisableLog('rdApp.warning')


"""
Update - May 6, 2024
This is a script which substitutes a given core molecule with the standard set
of terminals developed by Leticia Gomes. The script will perform substituions on 
the core at position indicated using Y and U. Then a sbatch job will be prepared
to submit CREST search.
Function of substitute a give core molecules with terminals adapted from
by Biruk Abreha and Steven Lopez. 

USAGE: python pymolgen.py molecule_name 'smiles'
"""

# validate inputs
try:
    molecule_name = sys.argv[1]
except:
    print("Please specify a molecule name and a SMILES string.")
    exit()
try:
    parent_smiles = sys.argv[2]
except:
    print("Please specify a SMILES string.")
    exit()

parent_mol = Chem.MolFromSmiles(parent_smiles)
if not parent_mol:
    print("Invalid SMILES string provided.")
    exit()

start_time = time.time()
current_directory = os.getcwd()

# reaction SMILES for linkers
linker_rxns = {'unsubstituted': '[*:1]([U])>>[*:1]([H])'}

# placeholder for linker addition
linker_place_holder = '[*:1]([U])'
linker_place_holder_mol = Chem.MolFromSmarts(linker_place_holder)

# append linkers to parent molecule to generate unsubstituted cores
unsubstituted_cores = []
place_holder_count = len(
    parent_mol.GetSubstructMatches(linker_place_holder_mol))
for linker in linker_rxns:
    rxn = AllChem.ReactionFromSmarts(linker_rxns[linker])
    core = parent_mol
    for i in range(place_holder_count):
        new_mols = list(chain.from_iterable(rxn.RunReactants((core,))))
        core = new_mols[0]
        Chem.SanitizeMol(core)
    unsubstituted_cores.append(core)


# reaction SMILES for terminal groups
terminal_rxns = {'hydrogen': '[*:1]([Y])>>[*:1]([H])',
                 'hydroxy': '[*:1]([Y])>>[*:1]([OH])',
                 'trifluoromethyl': '[*:1]([Y])>>[*:1][C](F)(F)F',
                 'trifluoromethoxy': '[*:1]([Y])>>[*:1][O][C](F)(F)F',
                 'methyl': '[*:1]([Y])>>[*:1][C]',
                 'methoxy': '[*:1]([Y])>>[*:1][O][C]',
                 'nitro': '[*:1]([Y])>>[*:1][N+]([O-])=O',
                 'thiol': '[*:1]([Y])>>[*:1]([SH])',
                 'fluoro': '[*:1]([Y])>>[*:1][F]',
                 'chloro': '[*:1]([Y])>>[*:1][Cl]',
                 'cyano': '[*:1]([Y])>>[*:1]C#N'}

substituent_place_holder = '[*:1]([Y])'
substituent_place_holder_mol = Chem.MolFromSmarts(substituent_place_holder)

# append terminal groups
all_mols = []
for core in unsubstituted_cores:
    place_holder_count = len(
        core.GetSubstructMatches(substituent_place_holder_mol))
    if place_holder_count == 0:
        all_mols.append(core)
        continue
    for terminal in terminal_rxns:
        new_mol = core
        rxn = AllChem.ReactionFromSmarts(terminal_rxns[terminal])
        for i in range(place_holder_count):
            new_mols = list(chain.from_iterable(rxn.RunReactants((new_mol,))))
            new_mol = new_mols[0]
            Chem.Cleanup(new_mol)
        all_mols.append(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)))

# canonicalize smiles to remove duplicates
all_mols = [Chem.MolFromSmiles(smiles) for smiles in [
    Chem.MolToSmiles(mol) for mol in all_mols]]
all_smiles = list(set([Chem.MolToSmiles(mol) for mol in all_mols]))

#Create ID
ID_list=list(range(1, len(all_smiles)+1))   


# create directory to store molecules
if not os.path.exists(molecule_name):
    os.makedirs(molecule_name)
out_folder = os.path.abspath(molecule_name)

# write list of SMILES to a cvs file with ID
smiles_id=pd.DataFrame(list(zip(ID_list,all_smiles)),columns=['ID','SMILES'])
smiles_id.to_csv(molecule_name+'_smiles-ID.csv', index=False)

#Generate InchKey
index=0
inchkey=[]
id_list=[]
path=[]
for smile in all_smiles:
    mol = Chem.AddHs(Chem.MolFromSmiles(smile))
    id_smile=ID_list[index]
    path_subdir=out_folder+'/'+str(id_smile)
    os.makedirs(path_subdir)
    path.append(path_subdir)
    id_list.append(id_smile)
    inch_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
    inchkey.append(inch_key)
    index+=1
    
inchkey_id=pd.DataFrame(list(zip(id_list,inchkey)),columns=['ID','InchKey'])
inchkey_id.to_csv(molecule_name+'_inchkey-ID.csv', index=False)
path_id=pd.DataFrame(list(zip(id_list,path)),columns=['ID','path'])

### Create dataframe with ID, InchKey and SMILES
step=pd.merge(inchkey_id,smiles_id, on= 'ID', how = 'inner')
step_final=pd.merge(path_id,step, on= 'ID', how = 'inner')

### Copy workflow CREST
#step 1
source_file = current_directory+'/workflow-CREST.py'
destination_file = out_folder+'/workflow-CREST.py'
shutil.copy(source_file, destination_file)

###Create sbatch files to submit CREST search
# step1
list_sleep=[20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1040, 1060, 1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320, 1340, 1360, 1380, 1400, 1420, 1440, 1460, 1480, 1500, 1520, 1540, 1560, 1580, 1600, 1620, 1640, 1660, 1680, 1700, 1720, 1740, 1760, 1780, 1800, 1820, 1840, 1860, 1880, 1900, 1920, 1940, 1960, 1980, 2000, 2020, 2040, 2060, 2080, 2100, 2120, 2140, 2160, 2180, 2200, 2220, 2240, 2260, 2280, 2300, 2320, 2340, 2360, 2380, 2400, 2420, 2440, 2460, 2480, 2500, 2520, 2540, 2560, 2580, 2600, 2620, 2640, 2660, 2680, 2700, 2720, 2740, 2760, 2780, 2800, 2820, 2840, 2860, 2880, 2900, 2920, 2940, 2960, 2980, 3000, 3020, 3040, 3060, 3080, 3100, 3120, 3140, 3160, 3180, 3200, 3220, 3240, 3260, 3280, 3300, 3320, 3340, 3360, 3380, 3400, 3420, 3440, 3460, 3480, 3500, 3520, 3540, 3560, 3580, 3600, 3620, 3640, 3660, 3680, 3700, 3720, 3740, 3760, 3780, 3800, 3820, 3840, 3860, 3880, 3900, 3920, 3940, 3960, 3980, 4000, 4020, 4040, 4060, 4080, 4100, 4120, 4140, 4160, 4180, 4200, 4220, 4240, 4260, 4280, 4300, 4320, 4340, 4360, 4380, 4400, 4420, 4440, 4460, 4480, 4500, 4520, 4540, 4560, 4580, 4600, 4620, 4640, 4660, 4680, 4700, 4720, 4740, 4760, 4780, 4800, 4820, 4840, 4860, 4880, 4900, 4920, 4940, 4960, 4980, 5000]
sbatch_title1=molecule_name+'_CREST.sbatch'
with open(sbatch_title1, 'w') as sbatch:
    sbatch.write("#!/bin/bash\n")
    sbatch.write("#SBATCH --job-name=CREST-search\n")
    sbatch.write("#SBATCH --partition=long\n")
    sbatch.write("#SBATCH --time=4-00:00:00\n")
    sbatch.write("#SBATCH --output=%j.o.slurm\n")
    sbatch.write("#SBATCH --error=%j.e.slurm\n")
    sbatch.write("source ~/.bashrc\n")
    sbatch.write("conda activate workflowV2_env\n")
    row_df=len(step_final.axes[0])
    for i in range(row_df):
        path=(step_final.iloc[i,1])
        inchkey=(step_final.iloc[i,2])
        smiles=(step_final.iloc[i,3])
        sbatch.write('cd '+path+' && cp ../workflow-CREST.py . && sbatch workflow-CREST.py '+inchkey+' \''+smiles+'\''+'\n')
        if i in list_sleep:
            sbatch.write("sleep 1500\n")
sbatch.close()     
      
