# Analysis_codes

Main script that calls all analysis codes
```
./main.sh <directory> <offset> <force field>
```

## Important subdirectories
1. process_trajectories: Centers the protein in trajectory file, extracts 0.5 ns frames, removes tags, standardises residue names, and renumbers residues from 1.
2. basic_properties: Gets RMSD, RMSF, Rg, Ramachandran plots, and thermodynamic properties for each protein. Does clustering analysis on structures and finds RMSDs decomposed according to secondary structure elements.
3. NMR_analysis: Calculates chemical shifts and J3 constants from the trajectory, compares them with experimental data and calculates statistics and maps relative errors to a .pkl file for visualising in pymol.
4. plotting_codes: Plots data obtained from basic_properties and NMR_analysis.
