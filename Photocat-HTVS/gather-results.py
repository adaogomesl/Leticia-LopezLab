import sys,os
import pandas as pd
from openbabel import openbabel
from rdkit import Chem

### Create by: Leticia A. Gomes
### Goal: Extract results from workflow and compute properties with sp correction
### Properties: VEE, Reduction and Oxidation Potential S0, S1 and T1
### Usage: python3 gather-results.py NAME-WORKFLOW


title='%s' % sys.argv[1].split('.')[0]
logpath=os.getcwd()
inchkey_vee=[]


#Lists to save VEE Results
vee=[]
wavelength=[]
osc=[]

#Lists to save S0 G Results
key_S0=[]
G_S0=[]
ZPVE_S0=[]

#Lists to save Anion G Results
key_AN=[]
G_AN=[]

#Lists to save Cation G Results
key_CAT=[]
G_CAT=[]


#Lists to save S1 ZPVE Results
key_S1=[]
ZPVE_S1=[]

#Lists to save T1 ZPVE Results
key_T1=[]
ZPVE_T1=[]

#Define types calculations to extract results
tddft='s0-sp-tddft-solv'
S0='s0-opt-freq-solv'
S1='s1-opt-freq-solv'
T1='t1-opt-freq-solv'
cation='cat-opt-freq-solv'
anion='1an-opt-freq-solv'


def get_path_pdb(title):
  name_workflow='%s' % (title)
  path_workflow='%s/%s' % (logpath,name_workflow)
  pdb_folder='unopt_pdbs'
  path='%s/%s' % (path_workflow,pdb_folder) 
  path_pdb=[]
  list_dir=os.listdir(path)
  for file in list_dir:
      if file.endswith('.pdb'):
         pathfile='%s/%s' % (path,file)
         path_pdb.append(pathfile)
  return path_pdb

def convert_pdb_files_to_smiles(path_pdb):
    inchkey_list=[]
    smiles_list=[]
    for file in path_pdb:
        # Create an Open Babel molecule
        ob_mol=openbabel.OBMol()

        # Read the PDB file
        conv=openbabel.OBConversion()
        conv.SetInFormat("pdb")
        conv.ReadFile(ob_mol,file)
        
        #Convert the molecule to SMILES
        conv.SetOutFormat("can")
        smiles=conv.WriteString(ob_mol).split('\t')[0]
        smiles_list.append(smiles)
        inchkey=conv.WriteString(ob_mol).split('unopt_pdbs/')[1].split('_')[0]
        inchkey_list.append(inchkey)
    smile_df=pd.DataFrame(list(zip(inchkey_list,smiles_list)),columns =['InchKey','SMILES'])
    smile_noduplicate=smile_df.drop_duplicates(subset=['InchKey'],keep='first',inplace=False)
    #os.chdir(logpath)
    #smile_noduplicate.to_csv('smiles.csv',index=False)
    return smile_noduplicate

#GET PATH
def get_path(calc):
  name_workflow='%s' % (title)
  path_workflow='%s/%s' % (logpath,name_workflow)
  path='%s/%s' % (path_workflow,calc) 
  wave=[]
  path_file=[]
  list_dir=os.listdir(path)
  for item in list_dir:
     check=item.__contains__('sbatch')
     if check==False:
        wave.append(item)
  for item in wave:
     path_wave='%s/%s/completed' % (path,item)
     os.chdir(path_wave)
     for file in os.listdir():
      if file.endswith('.log'):
         pathfile='%s/%s' % (path_wave,file)
         path_file.append(pathfile)
  return path_file


#VEE
def get_vee(path_vee):
   with open(path_vee, "r") as file:
      for line in file:
         if "Singlet" in line:
            if (float(line.split()[8][2:]))!=0:
               inchkey_vee.append((file.name.split("/")[-1]).split("_")[0])
               vee.append(float(line.split()[4]))
               wavelength.append(float(line.split()[6]))
               osc.append(line.split()[8][2:])
               break
            else:
               continue
   file.close()

def results_vee():
   path_file=get_path(tddft)
   for item in path_file:
    get_vee(item)
   vee_df = pd.DataFrame(list(zip(inchkey_vee,vee,wavelength,osc)),columns =['InchKey','VEE (eV)','Wavelength (nm)','Osc. Strength'])
   vee_df_noduplicate=vee_df.drop_duplicates(subset=['InchKey'],keep='first',inplace=False)
   #os.chdir(logpath)
   #vee_df_noduplicate.to_csv('vee.csv',index=False)
   return vee_df_noduplicate


#Get ZPVE
def get_ZPVE(path_file,key,energy):
    for item in path_file:
       with open(item, "r") as file:
        for line in file:
          if "Sum of electronic and zero-point Energies="in line:
              energy.append(float(line.split()[6]))
              key.append((file.name.split("/")[-1]).split("_")[0])
       file.close()
       
#Get Free Energy
def get_G(path_file,key,energy):
    for item in path_file:
       with open(item, "r") as file:
        for line in file:
          if "Sum of electronic and thermal Free Energies="in line:
              energy.append(float(line.split()[7]))
              key.append((file.name.split("/")[-1]).split("_")[0])
       file.close()
       

def red_S0():
    key_red_S0=[]
    red_S0=[]

    #Create dataframe with the extracted results
    S0_G=pd.DataFrame(list(zip(key_S0,G_S0)),columns =['InchKey','G_S0 (Hartree)'])
    AN_G=pd.DataFrame(list(zip(key_AN,G_AN)),columns =['InchKey','G_AN (Hartree)'])

    #Combine dataframes by the InchKey
    S0_red_raw = pd.merge(S0_G,AN_G, on= 'InchKey', how = 'inner')

    #Compute Reduction potential on ground state
    row_pd=len(S0_red_raw.axes[0])
    for i in range(row_pd):
        key_red_S0.append(S0_red_raw.iloc[i,0])
        #Compute Free Energy
        GS0=(S0_red_raw.iloc[i,1])
        GAN=(S0_red_raw.iloc[i,2])
        red_S0.append(((-(GAN-GS0))*27.2114)-4.429)
        
    S0red=pd.DataFrame(list(zip(key_red_S0,red_S0)),columns =['InchKey','Red_S0 (eV)'])
    S0_red=pd.merge(S0_red_raw,S0red, on= 'InchKey', how = 'inner')
    S0_red_noduplicate=S0_red.drop_duplicates(subset=['InchKey'],keep='first',inplace=False)
    S0red_noduplicate=S0red.drop_duplicates(subset=['InchKey'],keep='first',inplace=False)
    #os.chdir(logpath)
    #S0_red_noduplicate.to_csv('S0_red.csv',index=False)
    return S0red_noduplicate


def oxi_S0():
    key_oxi_S0=[]
    oxi_S0=[]

    #Create dataframe with the extracted results
    S0_G=pd.DataFrame(list(zip(key_S0,G_S0)),columns =['InchKey','G_S0 (Hartree)'])
    CAT_G=pd.DataFrame(list(zip(key_CAT,G_CAT)),columns =['InchKey','G_CAT (Hartree)'])

    #Combine dataframes by the InchKey
    S0_oxi_raw = pd.merge(S0_G,CAT_G, on= 'InchKey', how = 'inner')

    #Compute Oxidation potential on ground state
    row_pd=len(S0_oxi_raw.axes[0])
    for i in range(row_pd):
        key_oxi_S0.append(S0_oxi_raw.iloc[i,0])
        #Compute Free Energy
        GS0=(S0_oxi_raw.iloc[i,1])
        GCAT=(S0_oxi_raw.iloc[i,2])
        oxi_S0.append((((GCAT-GS0))*27.2114)-4.429)
        
    S0oxi=pd.DataFrame(list(zip(key_oxi_S0,oxi_S0)),columns =['InchKey','Oxi_S0 (eV)'])
    S0_oxi=pd.merge(S0_oxi_raw,S0oxi, on= 'InchKey', how = 'inner')
    S0_oxi_noduplicate=S0_oxi.drop_duplicates(subset=['InchKey'],keep='first',inplace=False)
    S0oxi_noduplicate=S0oxi.drop_duplicates(subset=['InchKey'],keep='first',inplace=False)
    #os.chdir(logpath)
    #S0_oxi_noduplicate.to_csv('S0_oxi.csv',index=False)
    return S0oxi_noduplicate

def Excited_potentials(key,energy,red_S0_final,oxi_S0_final,title):
    key_exc=[]
    zpve=[]
    key_final=[]
    oxi_exc=[]
    red_exc=[]

    #Create dataframe with the extracted results
    S0_ZPVE=pd.DataFrame(list(zip(key_S0,ZPVE_S0)),columns =['InchKey','ZPVE_S0 (Hartree)'])
    EXC_ZPVE=pd.DataFrame(list(zip(key,energy)),columns =['InchKey','ZPVE_Excited (Hartree)'])

    #Combine dataframes by the InchKey
    EXC_raw = pd.merge(S0_ZPVE,EXC_ZPVE, on= 'InchKey', how = 'inner')
    
    
    #Compute E0-0 Transition
    row_pd=len(EXC_raw.axes[0])
    for i in range(row_pd):
        key_exc.append(EXC_raw.iloc[i,0])
        ZPVES0=(EXC_raw.iloc[i,1])
        ZPVE_EXC=(EXC_raw.iloc[i,2])
        zpve.append((ZPVE_EXC-ZPVES0)*27.2114)
    zpvedf=pd.DataFrame(list(zip(key_exc,zpve)),columns =['InchKey','E00 (eV)'])
    zpve_df=pd.merge(EXC_raw,zpvedf, on= 'InchKey', how = 'inner')
                  
    #Compute Excited Potentials
    Ground_pot=pd.merge(red_S0_final,oxi_S0_final, on= 'InchKey', how = 'inner')
    EXC_pot=pd.merge(Ground_pot,zpve_df, on= 'InchKey', how = 'inner')
    row_pd=len(EXC_pot.axes[0])
    for i in range(row_pd):
        key_final.append(EXC_pot.iloc[i,0])
        #Compute Potentials
        EXCOxi=((EXC_pot.iloc[i,2])-(EXC_pot.iloc[i,5]))
        oxi_exc.append(EXCOxi)
        EXCRed=((EXC_pot.iloc[i,1])+(EXC_pot.iloc[i,5]))
        red_exc.append(EXCRed)
    excpot_final=pd.DataFrame(list(zip(key_final,oxi_exc,red_exc)),columns =['InchKey',title+'Oxi Potential (eV)',title+'Red Potential (eV)'])
    Exc_Potential=pd.merge(EXC_pot,excpot_final, on= 'InchKey', how = 'inner')
    Exc_Potential_noduplicate=Exc_Potential.drop_duplicates(subset=['InchKey'],keep='first',inplace=False)
    excpot_final_noduplicate=excpot_final.drop_duplicates(subset=['InchKey'],keep='first',inplace=False)
    #os.chdir(logpath)
    #titlecsv=title+'_Pot.csv'
    #Exc_Potential_noduplicate.to_csv(titlecsv,index=False)
    return excpot_final_noduplicate  

def main():
    smiles_results=convert_pdb_files_to_smiles(get_path_pdb(title))
    vee_df=results_vee()

    #S0 extraction
    get_ZPVE(get_path(S0),key_S0,ZPVE_S0)
    get_G(get_path(S0),key_S0,G_S0)

    #AN extraction
    get_G(get_path(anion),key_AN,G_AN)
    red_S0_final=red_S0()
    
    #CAT extraction
    get_G(get_path(cation),key_CAT,G_CAT)
    oxi_S0_final=oxi_S0()
 
    #S1 extraction
    get_ZPVE(get_path(S1),key_S1,ZPVE_S1)
    
    #T1 extraction
    get_ZPVE(get_path(T1),key_T1,ZPVE_T1)
    
    S1_pot=Excited_potentials(key_S1,ZPVE_S1,red_S0_final,oxi_S0_final,'S1')
    T1_pot=Excited_potentials(key_T1,ZPVE_T1,red_S0_final,oxi_S0_final,'T1')
    
    ground=pd.merge(oxi_S0_final,red_S0_final, on= 'InchKey', how = 'inner')
    excited=pd.merge(S1_pot,T1_pot, on= 'InchKey', how = 'inner')
    final_result=pd.merge(ground,excited, on= 'InchKey', how = 'inner')
    all_result=pd.merge(vee_df,final_result, on= 'InchKey', how = 'inner')
    all_result_withsmiles=pd.merge(smiles_results,all_result, on= 'InchKey', how = 'inner')
    os.chdir(logpath)
    all_result_withsmiles.to_csv('results-'+title+'.csv',index=False)
  
    
if __name__ == "__main__":
    main()
    
