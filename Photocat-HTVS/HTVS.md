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

