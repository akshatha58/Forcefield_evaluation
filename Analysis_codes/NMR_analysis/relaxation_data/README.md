```get_correlations.py``` computes the N-H bond correlations from the trajectory file using code from 
https://github.com/achicks15/CorrFunction_NMRRelaxation.git (with a few modifications to superpose trajectory frames)
```
python3 get_correlations.py <pdbpath>
```

```analyse_data.py``` computes the time averaged mean and standard deviations of N-H bond correlations for each residue in the protein, and plots them.
```
python3 analyse_data.py <pdbid> <flag> <forcefield>
```
