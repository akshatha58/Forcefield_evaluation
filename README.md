# Forcefield evaluation
Codes for implementing workflows involved in running simulations, processing trajectories, and analysing data for various force field/water model combinations.

## Main directories: 
1. process_input_files: Has code for processing topology files to tag residues and atom type
2. input_files: Has .mdp files and a script for running gromacs commands
3. Analysis_codes: Has codes to analyse structural, thermodynamic, and NMR data for all systems

## Data directories required:
1. NMR_BMRB: Store BMRB Files in STAR format for each protein
2. offset_lists and ss_lists: Store secondary structure information obtained from the DSSP (Helix and sheet residues for various proteins)
3. MD2NMR: Get from the repository: https://github.com/dkoes/MD2NMR
   Two separate repositories are created, one each for Amber and CHARMM (MD2NMR_amber and MD2NMR_charmm). dump.py and shifts.py have been modified to suit the use cases in this study. Get the other files and dictionaries from the above repository into each of these folders.

## Directory format for simulation files used in all codes:

<forcefield_name> / <pdb_id> / all simulation input files 

## Directory format created during analysis:

1. <forcefield_name> / <pdb_id> / Analysis / basic_properties
2. <forcefield_name> / <pdb_id> / Analysis / clusters
3. <forcefield_name> / <pdb_id> / Analysis / processed_trajs
4. <forcefield_name> / <pdb_id> / Analysis / NMR_analysis
5. <forcefield_name> / <pdb_id> / Analysis / NMR_analysis / chem_shifts
6. <forcefield_name> / <pdb_id> / Analysis / NMR_analysis / J3_constants
                                       

   
