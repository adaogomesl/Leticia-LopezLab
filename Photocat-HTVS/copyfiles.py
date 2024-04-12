import shutil 

source_dir = "/scratch/adaogomes.l/photocatalyst/workflow-photocat/"
workflow="bench_round"
#for i in range (9,10):
lista=['131','132','096']
for i in lista:
   #name_workflow='%s%d' % (workflow,i)
    name_workflow='%s%s' % (workflow,i)
    path_workflow='%s%s' % (source_dir,name_workflow)
    destination_dir = "/work/lopez/share_from_Leticia/photocatalyst/bench_potential_excited/"
    destination_path='%s%s' % (destination_dir,name_workflow)
    print(destination_path)
    shutil.copytree(path_workflow, destination_path)

