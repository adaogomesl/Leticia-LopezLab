import shutil
import sys,os
import pandas as pd

### Create by: Leticia A. Gomes
### Goal: Automatize preparation of files to resubmission
### Command to run: python3 failed-organize-resub.py NAME-WORKFLOW
# Change NAME-WORKFLOW to your workflow flow names

#***********************
#idim=your workflow number 
#idim=10
title='%s' % sys.argv[1].split('.')[0]
logpath=os.getcwd()
filenames=[]

#Define types of calculations that failed will be checked
jobs=['s0-sp-tddft-solv','s0-opt-freq-solv','s1-opt-freq-solv','t1-opt-freq-solv','cat-opt-freq-solv','1an-opt-freq-solv','cat-opt-freq-vac', 's0-opt-freq-vac']


#******************
#Do for loop
#for i in range(idim):
#name_workflow='%s-%d' % (title,i+1)
 

#GET PATH
def get_path(jobs):
  name_workflow='%s' % (title)
  path_workflow='%s/%s' % (logpath,name_workflow)
  path_file=[]
  wave=[]
  files=[]
  for item in jobs:
    path_folders='%s/%s' % (path_workflow,item)
    list_dir=os.listdir(path_folders)
    for item in list_dir:
      check=item.__contains__('sbatch')
      if check==False:
        wave.append(item)
    for item in wave:
      try#:
        path_wave='%s/%s/failed' % (path_folders,item)
        os.chdir(path_wave)
        for file in os.listdir():
          if file.endswith('.com'):
            pathfile='%s/%s' % (path_wave,file)
            path_file.append(pathfile)
            files.append(file.split(".")[0])
      except(FileNotFoundError):
        pass
    #remove duplicates
    file_path=[]
    for item in path_file:
      if item not in file_path:
        file_path.append(item)
    for item in files:
      if item not in filenames:
        filenames.append(item)
  return file_path

def copy_files(file_path):
  destination_dir = "/scratch/adaogomes.l/photocatalyst/workflow-photocat/resub/"
  destination_path='%s%s' % (destination_dir,title)
  os.makedirs(destination_path)

  for item in file_path:
    path_file=item
    shutil.copy(path_file, destination_path)
  return destination_path

def create_sbatch(filenames,destination_path):
  os.chdir(destination_path)
  sbatch_titles=[]
  for item in filenames:
    sbatch='%s.sbatch' % (item)
    sbatch_titles.append(sbatch)
  with open("runall.py", 'w') as runall:
    runall.write("#!/bin/sh\n")
    for item in sbatch_titles:
      runall.write('sbatch '+item+'\n')
  runall.close()

  for item in filenames:
    name=item+'.sbatch'
    with open(name, 'w') as sbatch:
     sbatch.write("#!/bin/bash\n")
     sbatch.write("#SBATCH --job-name="+item+'\n')
     sbatch.write("#SBATCH --input="+item+'.com\n')
     sbatch.write("#SBATCH --partition=loNG\n")
     sbatch.write("#SBATCH --time=10-00:00:00\n")
     sbatch.write("#SBATCH --nodes=1\n")
     sbatch.write("#SBATCH --ntasks=16\n")
     sbatch.write("#SBATCH --mem=300G\n")
     sbatch.write("#SBATCH --output=%j.o.slurm\n")
     sbatch.write("#SBATCH --error=%j.e.slurm\n")
     sbatch.write('\n')
     sbatch.write("hostname\n")
     sbatch.write("\n")
     sbatch.write("export g16root=/work/lopez/\n")
     sbatch.write(". $g16root/g16/bsd/g16.profile\n")
     sbatch.write("\n")
     sbatch.write("work=`pwd`\n")
     sbatch.write("\n")
     sbatch.write("export GAUSS_SCRDIR=$WORKING_DIR\n")
     sbatch.write("\n")
     sbatch.write("cd $work\n")
     sbatch.write("time $g16root/g16/g16 "+item+".com\n")
  sbatch.close()

            


def main():
    path=get_path(jobs)
    destination_path=copy_files(path)
    create_sbatch(filenames,destination_path)


if __name__ == "__main__":
    main()

