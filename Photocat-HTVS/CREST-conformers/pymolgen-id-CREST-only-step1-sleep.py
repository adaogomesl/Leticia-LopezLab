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
Update - April 12, 2024
This is a script which substitutes a given core molecule with the standard set
of terminals developed by Leticia Gomes. The script will perform substituions on 
the core at position indicated using Y and U. Then a sbatch job will be prepared
to submit CREST search.
Function of substitute a give core molecules with terminals adapted from
by Biruk Abreha and Steven Lopez. 

USAGE: python pymolgen-id-CREST molecule_name smiles.txt
"""

# Have a text file where each line is the smiles string for each step, in
# the following order part1 > part2 >part3

# validate inputs
try:
    molecule_name = sys.argv[1]
except:
    print("Please specify a molecule name and a SMILES string.")
    exit()
try:
    filename_smiles = sys.argv[2]
except:
    print("Please specify a txt file with SMILES string - Look script for formatting")
    exit()

#Extract smiles string
smiles_string=[] 
with open(filename_smiles, 'r') as file:
        for line in file:
            smiles_string.append(line.split()[0])   

part1_mol=Chem.MolFromSmiles(smiles_string[0])
part2_mol=Chem.MolFromSmiles(smiles_string[1])
part3_mol=Chem.MolFromSmiles(smiles_string[2])

if not part1_mol:
    print("Invalid SMILES string provided for part1")
elif not part2_mol:
    print("Invalid SMILES string provided for part2")
elif not part3_mol:
    print("Invalid SMILES string provided for part3")
    exit()

start_time = time.time()

# reaction SMILES for terminals
terminal1_rxns = {
    "hydrogen": "[*:1]([U])>>[*:1]([H])",
"methyl": "[*:1][U]>>[*:1]C",
"methoxy": "[*:1][U]>>[*:1]OC",
"tert-butyl": "[*:1][U]>>[*:1]C(C)(C)C",
"phenyl": "[*:1][U]>>[*:1]c1ccccc1", 
"fluoro": "[*:1]([U])>>[*:1][F]",
"chloro": "[*:1]([U])>>[*:1][Cl]",
"methoxycarbonyl": "[*:1][U]>>[*:1]C(OC)=O",
"carbamic": "[*:1][U]>>[*:1]C(N)=O",
"dimethylcarbamic": "[*:1][U]>>[*:1]C(N(C)C)=O",
"benzoyl": "[*:1][U]>>[*:1]C(C1=CC=CC=C1)=O",
"cyano": "[*:1]([U])>>[*:1]C#N",
"nitro": "[*:1]([U])>>[*:1][N+]([O-])=O",
"Acetoxy": "[*:1][U]>>[*:1]OC(C)=O",
"acetamido": "[*:1][U]>>[*:1]NC(C)=O",
"methylacetamido": "[*:1][U]>>[*:1]N(C)C(C)=O",
"trifluoromethyl": "[*:1][U]>>[*:1]C(F)(F)F",
"phenoxy": "[*:1][U]>>[*:1]OC1=CC=CC=C1",
"pentafluorophenyl": "[*:1][U]>>[*:1]C1=C(F)C(F)=C(F)C(F)=C1F",
"amino": "[*:1][U]>>[*:1]N",
"methylamino": "[*:1][U]>>[*:1]NC",
"dimethylamino": "[*:1][U]>>[*:1]N(C)C",
"pentafluorophenoxy": "[*:1][U]>>FC(C(F)=C(C(F)=C1F)F)=C1O[*:1]",
"methyl(phenyl)amino": "[*:1][U]>>[*:1]N(C)C1=CC=CC=C1",
"phenylamino": "[*:1][U]>>[*:1]NC1=CC=CC=C1",
"4-phenyl-1,2,3-triazolyl(1)": "[*:1][U]>>[*:1]N1N=NC(C2=CC=CC=C2)=C1",
"benzo[1,2,3]triazolyl": "[*:1][U]>>[*:1]N1N=NC2=C1C=CC=C2",
"tetrazol-5-yl": "[*:1][U]>>[*:1]C1=NNN=N1",
"pyrrol-1-yl": "[*:1][U]>>[*:1]N1C=CC=C1",
"imidazol-1-yl": "[*:1][U]>>[*:1]N1C=NC=C1",
"methylimidazol-1-yl": "[*:1][U]>>[*:1]C1=NC=CN1C",
"imidazol-2-yl": "[*:1][U]>>[*:1]C1=NC=CN1",
"oxazol-2-yl": "[*:1][U]>>[*:1]C1=NC=CO1",
"1-methyltetrazol-5-yl": "[*:1][U]>>[*:1]C1=NN(C)N=N1",
"indol-1-yl": "[*:1][U]>>[*:1]N1C(C=CC=C2)=C2C=C1",
"thiazol-2-yl": "[*:1][U]>>[*:1]C1=NC=CS1",
"fur-2-yl": "[*:1][U]>>[*:1]C1=CC=CO1",
"thien-2-yl": "[*:1][U]>>[*:1]C1=CC=CS1",
"indol-3-yl": "[*:1][U]>>[*:1]C1=CNC2=C1C=CC=C2",
"phenylureido": "[*:1][U]>>[*:1]NC(NC1=CC=CC=C1)=O",
"phenylsulfonyl": "[*:1][U]>>[*:1]S(=O)(C1=CC=CC=C1)=O",
"methanesulfonyl": "[*:1][U]>>[*:1]S(=O)(C)=O",
"sulfonamido": "[*:1][U]>>[*:1]S(=O)(N)=O",
"dimethylsulfonamido": "[*:1][U]>>[*:1]S(=O)(N(C)C)=O",
"1,4-trifluoromethylphenyl": "[*:1][U]>>[*:1]C1=CC=C(C(F)(F)F)C=C1",
"1,4-methoxyphenyl": "[*:1][U]>>[*:1]C1=CC=C(OC)C=C1",
"piperidyl": "[*:1][U]>>[*:1]N1CCCCC1",
"phosphonate": "[*:1][U]>>[*:1]P(OC)(OC)=O",
"morpholinyl": "[*:1][U]>>[*:1]N1CCOCC1",
    }
# placeholder for terminal1 addition
terminal1_place_holder = '[*:1]([U])'
terminal1_place_holder_mol = Chem.MolFromSmarts(terminal1_place_holder)

# append terminal1s to parent molecule to generate unsubstituted cores
unsubstituted_cores_part1 = []
unsubstituted_cores_part2 = []
unsubstituted_cores_part3 = []


#Check how many [U] there is on the core SMILES string 
place_holder_count = len(part1_mol.GetSubstructMatches(terminal1_place_holder_mol))

for terminal1 in terminal1_rxns:
    rxn = AllChem.ReactionFromSmarts(terminal1_rxns[terminal1])
    core1 = part1_mol
    core2 = part2_mol
    core3 = part3_mol
    for i in range(place_holder_count):
        new_mols1 = list(chain.from_iterable(rxn.RunReactants((core1,))))
        new_mols2 = list(chain.from_iterable(rxn.RunReactants((core2,))))
        new_mols3 = list(chain.from_iterable(rxn.RunReactants((core3,))))

        core1 = new_mols1[0]
        core2 = new_mols2[0]
        core3 = new_mols3[0]

        Chem.SanitizeMol(core1)
        Chem.SanitizeMol(core2)
        Chem.SanitizeMol(core3)

    unsubstituted_cores_part1.append(core1)
    unsubstituted_cores_part2.append(core2)
    unsubstituted_cores_part3.append(core3) 


# reaction SMILES for terminal groups
terminal2_rxns = {
    "hydrogen": "[*:1]([Y])>>[*:1]([H])",
"methyl": "[*:1][Y]>>[*:1]C",
"methoxy": "[*:1][Y]>>[*:1]OC",
"tert-butyl": "[*:1][Y]>>[*:1]C(C)(C)C",
"phenyl": "[*:1][Y]>>[*:1]c1ccccc1", 
"fluoro": "[*:1]([Y])>>[*:1][F]",
"chloro": "[*:1]([Y])>>[*:1][Cl]",
"methoxycarbonyl": "[*:1][Y]>>[*:1]C(OC)=O",
"carbamic": "[*:1][Y]>>[*:1]C(N)=O",
"dimethylcarbamic": "[*:1][Y]>>[*:1]C(N(C)C)=O",
"benzoyl": "[*:1][Y]>>[*:1]C(C1=CC=CC=C1)=O",
"cyano": "[*:1]([Y])>>[*:1]C#N",
"nitro": "[*:1]([Y])>>[*:1][N+]([O-])=O",
"Acetoxy": "[*:1][Y]>>[*:1]OC(C)=O",
"acetamido": "[*:1][Y]>>[*:1]NC(C)=O",
"methylacetamido": "[*:1][Y]>>[*:1]N(C)C(C)=O",
"trifluoromethyl": "[*:1][Y]>>[*:1]C(F)(F)F",
"phenoxy": "[*:1][Y]>>[*:1]OC1=CC=CC=C1",
"pentafluorophenyl": "[*:1][Y]>>[*:1]C1=C(F)C(F)=C(F)C(F)=C1F",
"amino": "[*:1][Y]>>[*:1]N",
"methylamino": "[*:1][Y]>>[*:1]NC",
"dimethylamino": "[*:1][Y]>>[*:1]N(C)C",
"pentafluorophenoxy": "[*:1][Y]>>FC(C(F)=C(C(F)=C1F)F)=C1O[*:1]",
"methyl(phenyl)amino": "[*:1][Y]>>[*:1]N(C)C1=CC=CC=C1",
"phenylamino": "[*:1][Y]>>[*:1]NC1=CC=CC=C1",
"4-phenyl-1,2,3-triazolyl(1)": "[*:1][Y]>>[*:1]N1N=NC(C2=CC=CC=C2)=C1",
"benzo[1,2,3]triazolyl": "[*:1][Y]>>[*:1]N1N=NC2=C1C=CC=C2",
"tetrazol-5-yl": "[*:1][Y]>>[*:1]C1=NNN=N1",
"pyrrol-1-yl": "[*:1][Y]>>[*:1]N1C=CC=C1",
"imidazol-1-yl": "[*:1][Y]>>[*:1]N1C=NC=C1",
"methylimidazol-1-yl": "[*:1][Y]>>[*:1]C1=NC=CN1C",
"imidazol-2-yl": "[*:1][Y]>>[*:1]C1=NC=CN1",
"oxazol-2-yl": "[*:1][Y]>>[*:1]C1=NC=CO1",
"1-methyltetrazol-5-yl": "[*:1][Y]>>[*:1]C1=NN(C)N=N1",
"indol-1-yl": "[*:1][Y]>>[*:1]N1C(C=CC=C2)=C2C=C1",
"thiazol-2-yl": "[*:1][Y]>>[*:1]C1=NC=CS1",
"fur-2-yl": "[*:1][Y]>>[*:1]C1=CC=CO1",
"thien-2-yl": "[*:1][Y]>>[*:1]C1=CC=CS1",
"indol-3-yl": "[*:1][Y]>>[*:1]C1=CNC2=C1C=CC=C2",
"phenylureido": "[*:1][Y]>>[*:1]NC(NC1=CC=CC=C1)=O",
"phenylsulfonyl": "[*:1][Y]>>[*:1]S(=O)(C1=CC=CC=C1)=O",
"methanesulfonyl": "[*:1][Y]>>[*:1]S(=O)(C)=O",
"sulfonamido": "[*:1][Y]>>[*:1]S(=O)(N)=O",
"dimethylsulfonamido": "[*:1][Y]>>[*:1]S(=O)(N(C)C)=O",
"phosphonate": "[*:1][Y]>>[*:1]P(OC)(OC)=O",
"1,4-trifluoromethylphenyl": "[*:1][Y]>>[*:1]C1=CC=C(C(F)(F)F)C=C1",
"1,4-methoxyphenyl": "[*:1][Y]>>[*:1]C1=CC=C(OC)C=C1",
"piperidyl": "[*:1][Y]>>[*:1]N1CCCCC1",
"morpholinyl": "[*:1][Y]>>[*:1]N1CCOCC1",
    }
substituent_place_holder = '[*:1]([Y])'
substituent_place_holder_mol = Chem.MolFromSmarts(substituent_place_holder)

# append terminal groups
all_mols1 = []
all_mols2 = []
all_mols3 = []

for core in unsubstituted_cores_part1:
    place_holder_count = len(
        core.GetSubstructMatches(substituent_place_holder_mol))
    if place_holder_count == 0:
        all_mols1.append(core)
        continue
    for terminal in terminal2_rxns:
        new_mol = core
        rxn = AllChem.ReactionFromSmarts(terminal2_rxns[terminal])
        for i in range(place_holder_count):
            new_mols = list(chain.from_iterable(rxn.RunReactants((new_mol,))))
            new_mol = new_mols[0]
            Chem.Cleanup(new_mol)
        all_mols1.append(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)))
        
for core in unsubstituted_cores_part2:
    place_holder_count = len(
        core.GetSubstructMatches(substituent_place_holder_mol))
    if place_holder_count == 0:
        all_mols2.append(core)
        continue
    for terminal in terminal2_rxns:
        new_mol = core
        rxn = AllChem.ReactionFromSmarts(terminal2_rxns[terminal])
        for i in range(place_holder_count):
            new_mols = list(chain.from_iterable(rxn.RunReactants((new_mol,))))
            new_mol = new_mols[0]
            Chem.Cleanup(new_mol)
        all_mols2.append(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)))
        
        
for core in unsubstituted_cores_part3:
    place_holder_count = len(
        core.GetSubstructMatches(substituent_place_holder_mol))
    if place_holder_count == 0:
        all_mols3.append(core)
        continue
    for terminal in terminal2_rxns:
        new_mol = core
        rxn = AllChem.ReactionFromSmarts(terminal2_rxns[terminal])
        for i in range(place_holder_count):
            new_mols = list(chain.from_iterable(rxn.RunReactants((new_mol,))))
            new_mol = new_mols[0]
            Chem.Cleanup(new_mol)
        all_mols3.append(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)))      
        
        
# canonicalize smiles and remove duplicates
all_mols1 = [Chem.MolFromSmiles(smiles) for smiles in [Chem.MolToSmiles(mol) for mol in all_mols1]]
all_mols2 = [Chem.MolFromSmiles(smiles) for smiles in [Chem.MolToSmiles(mol) for mol in all_mols2]]
all_mols3 = [Chem.MolFromSmiles(smiles) for smiles in [Chem.MolToSmiles(mol) for mol in all_mols3]]


all_smiles1 = ([Chem.MolToSmiles(mol) for mol in all_mols1])
all_smiles2 = ([Chem.MolToSmiles(mol) for mol in all_mols2])
all_smiles3 = ([Chem.MolToSmiles(mol) for mol in all_mols3])



checked_smiles=[]
index_to_remove=[]
for i in range(len(all_smiles1)):
    mol_to_find=all_smiles1[i]
    if mol_to_find not in checked_smiles:
       checked_smiles.append(mol_to_find)
       indices = [idx for idx in range(len(all_smiles1)) if all_smiles1[idx] == mol_to_find]
       first_idx=indices[0]
       for i in indices:
           if i!=first_idx and all_smiles2[first_idx]==all_smiles2[i] and all_smiles3[first_idx]==all_smiles3[i]:
              index_to_remove.append(i)

all_smiles1 = [elem for idx, elem in enumerate(all_smiles1) if idx not in index_to_remove]
all_smiles2 = [elem for idx, elem in enumerate(all_smiles2) if idx not in index_to_remove]
all_smiles3 = [elem for idx, elem in enumerate(all_smiles3) if idx not in index_to_remove]
ID_list=list(range(1, len(all_smiles1)+1))   
     

# create directory to store molecules
molecule_name1=molecule_name+'-part1'
if not os.path.exists(molecule_name1):
    os.makedirs(molecule_name1)
out_folder1 = os.path.abspath(molecule_name1)

##molecule_name2=molecule_name+'-part2'
##if not os.path.exists(molecule_name2):
##    os.makedirs(molecule_name2)
##out_folder2 = os.path.abspath(molecule_name2)
##
##molecule_name3=molecule_name+'-part3'
##if not os.path.exists(molecule_name3):
##    os.makedirs(molecule_name3)
##out_folder3 = os.path.abspath(molecule_name3)

# write list of SMILES to a cvs file with ID
smiles_id1=pd.DataFrame(list(zip(ID_list,all_smiles1)),columns=['ID','SMILES-part1'])
smiles_id1.to_csv(molecule_name1+'_smiles-ID.csv', index=False)

##smiles_id2=pd.DataFrame(list(zip(ID_list,all_smiles2)),columns=['ID','SMILES-part2'])
##smiles_id2.to_csv(molecule_name2+'_smiles-ID.csv', index=False)
##
##smiles_id3=pd.DataFrame(list(zip(ID_list,all_smiles3)),columns=['ID','SMILES-part3'])
##smiles_id3.to_csv(molecule_name3+'_smiles-ID.csv', index=False)

#Generate InchKey
#step1
index=0
inchkey=[]
id_list=[]
path1=[]
for smile in all_smiles1:
    mol = Chem.AddHs(Chem.MolFromSmiles(smile))
    id_smile=ID_list[index]
    path_subdir=out_folder1+'/'+str(id_smile)
    os.makedirs(path_subdir)
    path1.append(path_subdir)
    id_list.append(id_smile)
    inch_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
    inchkey.append(inch_key)
    index+=1
    
inchkey_id1=pd.DataFrame(list(zip(id_list,inchkey)),columns=['ID','InchKey'])
inchkey_id1.to_csv(molecule_name1+'_inchkey-ID.csv', index=False)
path_id1=pd.DataFrame(list(zip(id_list,path1)),columns=['ID','path'])

###step2 
##index=0
##inchkey=[]
##id_list=[]
##path2=[]
##for smile in all_smiles2:
##    mol = Chem.AddHs(Chem.MolFromSmiles(smile))
##    id_smile=ID_list[index]
##    path_subdir=out_folder2+'/'+str(id_smile)
##    os.makedirs(path_subdir)
##    path2.append(path_subdir)
##    id_list.append(id_smile)
##    inch_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
##    inchkey.append(inch_key)
##    index+=1
##    
##inchkey_id2=pd.DataFrame(list(zip(id_list,inchkey)),columns=['ID','InchKey'])
##inchkey_id2.to_csv(molecule_name2+'_inchkey-ID.csv', index=False)
##path_id2=pd.DataFrame(list(zip(id_list,path2)),columns=['ID','path'])
##
###step3
##index=0
##inchkey=[]
##id_list=[]
##path3=[]
##for smile in all_smiles3:
##    mol = Chem.AddHs(Chem.MolFromSmiles(smile))
##    id_smile=ID_list[index]
##    path_subdir=out_folder3+'/'+str(id_smile)
##    os.makedirs(path_subdir)
##    path3.append(path_subdir)
##    id_list.append(id_smile)
##    inch_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
##    inchkey.append(inch_key)
##    index+=1
##    
##inchkey_id3=pd.DataFrame(list(zip(id_list,inchkey)),columns=['ID','InchKey'])
##inchkey_id3.to_csv(molecule_name3+'_inchkey-ID.csv', index=False)
##path_id3=pd.DataFrame(list(zip(id_list,path3)),columns=['ID','path'])

### Create dataframe with ID, InchKey and SMILES
step1=pd.merge(inchkey_id1,smiles_id1, on= 'ID', how = 'inner')
step1_final=pd.merge(path_id1,step1, on= 'ID', how = 'inner')

##step2=pd.merge(inchkey_id2,smiles_id2, on= 'ID', how = 'inner')
##step2_final=pd.merge(path_id2,step2, on= 'ID', how = 'inner')
##
##step3=pd.merge(inchkey_id3,smiles_id3, on= 'ID', how = 'inner')
##step3_final=pd.merge(path_id3,step3, on= 'ID', how = 'inner')



### Copy workflow CREST
#step 1
source_file = '/work/lopez/share_from_Leticia/forSteven/scripts/workflow-step1.py'
destination_file = out_folder1+'/workflow-step1.py'
shutil.copy(source_file, destination_file)

###step 2
##source_file = '/work/lopez/share_from_Leticia/forSteven/scripts/workflow-step2.py'
##destination_file = out_folder2+'/workflow-step2.py'
##shutil.copy(source_file, destination_file)
##
###step 3
##source_file = '/work/lopez/share_from_Leticia/forSteven/scripts/workflow-step3.py'
##destination_file = out_folder3+'/workflow-step3.py'
##shutil.copy(source_file, destination_file)


###Create sbatch files to submit CREST search
# step1
list_sleep=[50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500, 2550, 2600, 2650, 2700, 2750, 2800, 2850, 2900, 2950, 3000, 3050, 3100, 3150, 3200, 3250, 3300, 3350, 3400, 3450, 3500, 3550, 3600, 3650, 3700, 3750, 3800, 3850, 3900, 3950, 4000]
sbatch_title1=molecule_name+'-part1.sbatch'
with open(sbatch_title1, 'w') as sbatch:
    sbatch.write("#!/bin/bash\n")
    sbatch.write("#SBATCH --job-name=CREST-search\n")
    sbatch.write("#SBATCH --partition=long\n")
    sbatch.write("#SBATCH --time=2-00:00:00\n")
    sbatch.write("#SBATCH --output=%j.o.slurm\n")
    sbatch.write("#SBATCH --error=%j.e.slurm\n")
    sbatch.write("source ~/.bashrc\n")
    sbatch.write("conda activate workflowV2_env\n")
    row_df=len(step1_final.axes[0])
    for i in range(row_df):
        path=(step1_final.iloc[i,1])
        inchkey=(step1_final.iloc[i,2])
        smiles=(step1_final.iloc[i,3])
        sbatch.write('cd '+path+' && cp ../workflow-step1.py . && sbatch workflow-step1.py '+inchkey+' \''+smiles+'\''+'\n')
        if i in list_sleep:
            sbatch.write("sleep 1800\n")
      
sbatch.close()

### step2
##sbatch_title2=molecule_name+'-part2.sbatch'
##with open(sbatch_title2, 'w') as sbatch:
##    sbatch.write("#!/bin/bash\n")
##    row_df=len(step2_final.axes[0])
##    for i in range(row_df):
##        path=(step2_final.iloc[i,1])
##        inchkey=(step2_final.iloc[i,2])
##        smiles=(step2_final.iloc[i,3])
##        sbatch.write('cd '+path+' && cp ../workflow-step2.py . && sbatch workflow-step2.py '+inchkey+' \''+smiles+'\''+'\n')
##sbatch.close()
##
### step3
##sbatch_title3=molecule_name+'-part3.sbatch'
##with open(sbatch_title3, 'w') as sbatch:
##    sbatch.write("#!/bin/bash\n")
##    row_df=len(step3_final.axes[0])
##    for i in range(row_df):
##        path=(step3_final.iloc[i,1])
##        inchkey=(step3_final.iloc[i,2])
##        smiles=(step3_final.iloc[i,3])
##        sbatch.write('cd '+path+' && cp ../workflow-step3.py . && sbatch workflow-step3.py '+inchkey+' \''+smiles+'\''+'\n')
##sbatch.close()

##### Write runall
##sbatch_title='runall.sh'
##with open(sbatch_title, 'w') as sbatch:
##    sbatch.write("#!/bin/bash\n")
##    sbatch.write("#SBATCH --job-name=CREST-search\n")
##    sbatch.write("#SBATCH --partition=short,lopez\n")
##    sbatch.write("#SBATCH --time=1-00:00:00\n")
##    sbatch.write("#SBATCH --nodes=1\n")
##    sbatch.write("#SBATCH --ntasks=16\n")
##    sbatch.write("#SBATCH --mem=300G\n")
##    sbatch.write("#SBATCH --output=%j.o.slurm\n")
##    sbatch.write("#SBATCH --error=%j.e.slurm\n")
##    sbatch.write("source ~/.bashrc\n")
##    sbatch.write("conda activate workflowV2_env\n")
##    sbatch.write("sbatch "+sbatch_title1+'\n')
##    sbatch.write("sbatch "+sbatch_title2+'\n')
##    sbatch.write("sbatch "+sbatch_title3+'\n')
##sbatch.close()
    

