#!/usr/bin/env python
#SBATCH --job-name=workflow_CREST
#SBATCH --output=out
#SBATCH --error=error
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=20:00:00


solvent = 'Acetonitrile'
method = 'PBE1PBE/6-311+G(d,p)'
nproc = 16
mem = 10
max_confs=10

"""
USAGE: sbatch workflow.py molecule_name smiles
"""

###################################################################

from workflowV2 import molecule
from workflowV2.calculator import Run
from workflowV2.software.GAUSSIAN import GAUSSIAN
from workflowV2.software.CREST import CREST
from workflowV2 import message
import pandas as pd
import sys
import os


# validate inputs
try:
    input_name = sys.argv[1]
except:
    print("Please specify a molecule name and a SMILES string.")
    exit()
try:
    parent_smiles = sys.argv[2]
except:
    print("Please specify a SMILES string.")
    exit()


#log the output to a file
message.logtofile('{0}.log'.format(input_name))

##################################################################

#set up the molecule object
mol = molecule.SmilesToMol (parent_smiles,charge=0,mult=1)

#conformational search
os.makedirs('conf_search',exist_ok=True)
os.chdir('conf_search')

conformer_calculator = CREST(mol,jobname='{0}-confsearch'.format(input_name),runtype='confsearch',nproc=nproc,mem=mem,gbsa=solvent)

mol = Run(conformer_calculator)

os.chdir('../')

for i in range(max_confs):
   mol.conformers[i].ToXYZ(input_name+'_'+str(i)+'.xyz'.format(input_name))

print('Done!')


