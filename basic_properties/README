## basic_properties directory

1. `main_properties.sh` : main call (./main_properties.sh <directory> )
    a. Computes potential, total energy, temperature, pressure and density across simulation time. 
       Statistics are stored in `metrics.txt` in the results directory
    b. Computes Ramachandran plot, backbone RMSD and RMSFs, and radii of gyration from the 
       trajectory
    c. Computes backbone RMSDs decomposed according to secondary structure (calls `rmsd_ss.sh`)
    d. Computes clusters of structures with and RMSD less than or equal to 0.25 (calls `get_clusters.sh`)

2. `get_ss_groups.py` and `ndx_groups.py` aid in obtaining and formatting residue ranges for secondary structures (according to the NMR structure) to update the index file. Output is stored in `groups.txt` in the working directory. 
