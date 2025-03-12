"""
compare_shifts.py
Akshatha
7 October 2024

Imports chemical shifts data from MD and experiment and compares them (residue-wise), with std devs
ID  Res  N_exp  H_exp  C_exp  N_sd_exp  H_sd_exp  C_sd_exp  N_md  H_md  C_md  N_sd_md  H_sd_md  C_sd_md  ss_info

python3 compare_shifts.py <path>
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import global_plots as pltfn
import os

path = sys.argv[1]

def plot_shifts(path):
    """
    Plot chemical shift comparisons between NMR and MD 
    """
    pdb = os.path.splitext(os.path.basename(path))[0]
    datafile = path+"/Analysis/NMR_analysis/chem_shifts/"+pdb+"_expshifts_combined.txt"

    # Plot relative errors and exp/md scatter
    pltfn.consolidate_plots(datafile, "scatter")
    pltfn.consolidate_plots(datafile, "all")

def get_statistics(path):
    """
    Get chi squared errors, RMSD and Rsquared values from data
    """
    pdb = os.path.splitext(os.path.basename(path))[0]
    datafile = path+"/Analysis/NMR_analysis/chem_shifts/"+pdb+"_expshifts_combined.txt"
    directory_path = os.path.dirname(datafile)

    params = ['H', 'C', 'N']
    whichss = ['all', 'H', 'S', 'F']
    with open(directory_path+"/"+pdb+"_statistics.txt", "w") as f: 
        # f.write("res\tatom\tchi-sq\trmsd\tR-sq\trel-err\tmean_nmr\tmean_md\n")
        f.write("res\tatom\tchi-sq\trmsd\trel-err\tmean_nmr\tmean_md\n")
        for param in params:
            for ss in whichss:
                # chi_sq, rmsd, r_sq, avg_rel_err, mean_nmr, mean_md = pltfn.stats(datafile, param, ss, "shifts")
                chi_sq, rmsd, avg_rel_err, mean_nmr, mean_md = pltfn.stats(datafile, param, ss, "shifts")
                if chi_sq > 0:
                    # f.write(f"{ss}\t{param}\t{chi_sq:.2f}\t{rmsd:.2f}\t{r_sq:.2f}\t{avg_rel_err:.2f}\t{mean_nmr:.2f}\t{mean_md:.2f}\n")
                    f.write(f"{ss}\t{param}\t{chi_sq:.2f}\t{rmsd:.2f}\t{avg_rel_err:.2f}\t{mean_nmr:.2f}\t{mean_md:.2f}\n")
    f.close()
    
plot_shifts(path)
get_statistics(path)


