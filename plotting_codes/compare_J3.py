"""
compare_J3.py
Akshatha
6 November 2024

Imports J3 data from MD and experiment and compares them (residue-wise), with std devs

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
    datafile = path+"/Analysis/NMR_analysis/J3_constants/"+pdb+"_Jcombined.txt"
    # Plot relative errors and exp/md scatter
    pltfn.plot_J3_results(datafile)
    pltfn.plot_J3_scatter(datafile)


# MENTION RESIDUES THAT HAVE THE MAXIMUM RELATIVE ERROR BEYOND CERTAIN STANDARD DEVIATIONS FROM THE OVERALL MEAN VALUE

def get_statistics(path):
    """
    Get chi squared errors, RMSD and Rsquared values from data
    """
    pdb = os.path.splitext(os.path.basename(path))[0]
    datafile = path+"/Analysis/NMR_analysis/J3_constants/"+pdb+"_Jcombined.txt"
    directory_path = os.path.dirname(datafile)

    param="N" #dummy param
    whichss = ['all', 'H', 'S', 'F']
    with open(directory_path+"/"+pdb+"_statistics.txt", "w") as f: 
        f.write("res\tchi-sq\trmsd\tR-sq\trel-err\tmean_nmr\tmean_md\n")
        for ss in range(len(whichss)):
            chi_sq, rmsd, r_sq, avg_rel_err, mean_nmr, mean_md = pltfn.stats(datafile, param, whichss[ss], "J3")
            if chi_sq > 0:
                f.write(f"{whichss[ss]}\t{chi_sq:.2f}\t{rmsd:.2f}\t{r_sq:.2f}\t{avg_rel_err:.2f}\t{mean_nmr:.2f}\t{mean_md:.2f}\n")

    f.close()


    
plot_shifts(path)
get_statistics(path)


