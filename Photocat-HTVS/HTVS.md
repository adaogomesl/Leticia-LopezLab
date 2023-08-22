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

