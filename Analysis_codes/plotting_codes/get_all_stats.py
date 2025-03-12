"""
get_all_stats.py
6 March 2025

Tabulates RMSDs and Relative errors across proteins and force fields, for various NMR observables, 
decomposed according to DSSP secondary structure elements for further analysis

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.patches as mpatches
import sys
import os
import global_plots as gplots

def compute_stat(stat, nmrdata, mddata, mdsd):
    """
    computes various statistics for the given data
    """
   
    if len(nmrdata != 0):
        if stat == "RMSD":
            return np.sqrt(np.mean((nmrdata - mddata) ** 2))
        elif stat == "chisq":
            chi_sq = 0
            for idx in range(len(nmrdata)):
                chi_sq += ((nmrdata[idx] - mddata[idx])**2)/mdsd[idx]**2
            return chi_sq/len(nmrdata)
    else:
        return -800

def stats(datafile, params, whichss, prop, stat):
    """
    Get statistics for the data
    """
    # Get parameter-specific data 
    if prop == "shifts":
        parameters = gplots.return_paramdata(datafile, params, 'shifts')
        data = pd.DataFrame(parameters, columns=["ID", "res", params+"_exp", params+"_sd_exp", params+"_md", params+"_sd_md", "ss_info", "rel_err"+params])

        # Rename columns to make names easier to read and handle
        if whichss == 0:
            nmr_data = np.array(data[params+"_exp"])
            md_data = np.array(data[params+"_md"])
            md_sd = np.array(data[params+"_sd_md"])
            nmr_sd =np.array(data[params+"_sd_exp"])
            rel_err = np.array(data["rel_err"+params])

        else:
            # select only those rows whose last column equals whichss
            # filtered_df = df[df.iloc[:, -1] == 'H']
            filtered_data = data[data.iloc[:,-2] == whichss]
            nmr_data = np.array(filtered_data[params+"_exp"])
            md_data = np.array(filtered_data[params+"_md"])
            md_sd = np.array(filtered_data[params+"_sd_md"])
            nmr_sd = np.array(filtered_data[params+"_sd_exp"])
            rel_err = np.array(filtered_data["rel_err"+params])
            
    elif prop == "J3":
        params='NA' # dummy variable
        parameters = gplots.return_paramdata(datafile, params, "J3")
        data = pd.DataFrame(parameters, columns=["ID", "Res", "J_exp", "J_sd_exp", "J_md", "J_sd_md", "ss_info", "rel_err"])

         # Rename columns to make names easier to read and handle
        if whichss == 0:
            nmr_data = np.array(data["J_exp"])
            md_data = np.array(data["J_md"])
            md_sd = np.array(data["J_sd_md"])
            nmr_sd = data["J_sd_exp"]
            rel_err = np.array(data["rel_err"])

        else:
            # select only those rows whose last column equals whichss
            # filtered_df = df[df.iloc[:, -1] == 'H']
            filtered_data = data[data.iloc[:,-2] == whichss]
            nmr_data = np.array(filtered_data["J_exp"])
            md_data = np.array(filtered_data["J_md"])
            nmr_sd = np.array(filtered_data["J_sd_exp"])
            md_sd = np.array(filtered_data["J_sd_md"])
            rel_err = np.array(filtered_data["rel_err"])

    if stat == "relerr":
        if len(nmr_data != 0):
            return np.mean(rel_err)
        else:
            return -800
    else:
        return compute_stat(stat, nmr_data, md_data, md_sd)
    
def init_dict():
    """
    Initialises a dictionary having average statistics for all secondary structure elements
    """
    metric = {}
    for i in range(0,13):
        metric[i] = -800
    return metric

def get_all_statistics(stat, prop, params):
    """
    Write stats (RMSD or Relative errors) for proton chemical shifts for each protein and force field into a file
    """
    pdblist = ['1ozi', '1r36','1z9b', '2gmo', '2juo', '2jxy', '2k0q', '2l4x', '2l6b', '2lht', '2lro', '2m68', '4bwh']
    fflist = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    outputfile = "comparisons/heatmaps/"+stat+"_"+params+".txt"
    whichss = range(0, 13)

    with open(outputfile, "w") as f: 
        f.write("pdb\tff\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\n")
        
        for pdb in pdblist:
            for ff in fflist:
                
                metric = init_dict()

                if prop == "shifts":
                    datafile = "../../"+ff+"/"+pdb+"/Analysis/NMR_analysis/chem_shifts/"+pdb+"_expshifts_combined.txt"
                elif prop == "J3":
                    datafile = "../../"+ff+"/"+pdb+"/Analysis/NMR_analysis/J3_constants/"+pdb+"_Jcombined.txt"

                if os.path.exists(datafile):
                   
                    for ss in whichss:
                        metric[ss] = stats(datafile, params, ss, prop, stat)

                    values = list(metric.values())
                    f.write(pdb+"\t"+ff+"\t")
                    f.write(' '.join(f'{value:.3f}' for value in values))
                    f.write("\n")
        f.close()

get_all_statistics("chisq", "J3", "J3")
get_all_statistics("chisq", "shifts", "H")
get_all_statistics("chisq", "shifts", "C")
get_all_statistics("chisq", "shifts", "N")

get_all_statistics("relerr", "J3", "J3")
get_all_statistics("relerr", "shifts", "H")
get_all_statistics("relerr", "shifts", "C")
get_all_statistics("relerr", "shifts", "N")

get_all_statistics("RMSD", "J3", "J3")
get_all_statistics("RMSD", "shifts", "H")
get_all_statistics("RMSD", "shifts", "C")
get_all_statistics("RMSD", "shifts", "N")
