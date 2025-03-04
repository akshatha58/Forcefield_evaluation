# Forcefield evaluation
Codes for implementing workflows involved in running simulations, processing trajectories, and analysing data for various force field/water model combinations.

Main directories: 
1. process_input_files: Has code for processing topology files to tag residues and atom type
2. input_files: Has .mdp files and a script for running gromacs commands
3. Analysis_codes: Has codes to analyse structural, thermodynamic, and NMR data for all systems

Data directories required:
1. NMR_BMRB: Store BMRB Files in STAR format for each protein
2. offset_lists and ss_lists: Store secondary structure information obtained from the DSSP (Helix and sheet residues for various proteins)
3. MD2NMR_amber and MD2NMR_charmm: Get from the repository: https://github.com/dkoes/MD2NMR 

Directory format for simulation files used in all codes:
<forcefield_name> / <pdbid> / *all simulation input files 

Directory format created during analysis:
<forcefield_name> / <pdbid> / Analysis / basic_properties
<forcefield_name> / <pdbid> / Analysis / clusters
<forcefield_name> / <pdbid> / Analysis / processed_trajs
<forcefield_name> / <pdbid> / Analysis/ NMR_analysis
<forcefield_name> / <pdbid> / Analysis/ NMR_analysis / chem_shifts
<forcefield_name> / <pdbid> / Analysis/ NMR_analysis / J3_constants
                                       

   
