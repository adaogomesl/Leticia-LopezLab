export INPUT=input.json
export WORKDIR=/scratch/adaogomes.l/DBH/unsub/ML-NAMD/step2_grid_search/grid-search/NN-dbh-3_200_64_1e-08_1e-08_0.005

cd $WORKDIR
pyrai2md $INPUT
