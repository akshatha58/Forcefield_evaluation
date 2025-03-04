"""
map_devs.py
15 November 2024

Maps residues that have an absolute value deviation in J3 coupling constant
To be used in pymol functionality

python3 map_devs.py <path>
"""
import pandas as pd
import numpy as np
import sys, os
import pickle

basepath=sys.argv[1]

sys.path.insert(0, basepath+'/Analysis_optimised/plotting_codes')

import global_plots as gplt

proteins = ['1ozi', '2gmo', '1z9b', '2juo']
# classes = { 'all_beta: ['1ozi'], 'alpha_beta':['2gmo', '1z9b','2juo']}
# classes = {'all_alpha': ['2l6b', '1r36], 'all_beta: ['2jxy', '2k0q','2m68','1ozi','2lro'], 'alpha_beta':['2gmo', '1z9b', '2l4x','4bwh','2juo'], 'spl':['2lht'], 'few_ss':['2ec7']}

prots_all = ['1ozi', '1r36', '2ec7', '2gmo', '2juo', '2jxy', '2k0q', '2l4x', '2l6b', '2lro', '4bwh']


def return_idx(protein):
    """
    Returns offset indices so that the right residue numbers are mapped to the pdb structure file
    """
    if protein=="1ozi":
        offset = 1
    elif protein=="2gmo":
        offset = -7
    elif protein=="1z9b":
        offset = 1
    elif protein=="2juo":
        offset = 0
    return offset
    
resdict={}

for protein in proteins:
    datafile = basepath+protein+"/Analysis/NMR_analysis/J3_constants/"+protein+"_Jcombined.txt"
    directory_path = os.path.dirname(datafile)
    params='NA' # dummy variable
    parameters = gplt.return_paramdata(datafile, params, "J3")
    data = pd.DataFrame(parameters, columns=["ID", "Res", "J_exp", "J_sd_exp", "J_md", "J_sd_md", "ss_info", "rel_err"])

    offset = return_idx(protein)
    # Rename columns to make names easier to read and handle
    resids = data["ID"]
    resname = data["Res"]
    nmr_data = data["J_exp"]
    nmr_sd = data["J_sd_exp"]
    md_data = data["J_md"]
    md_sd = data["J_sd_md"]
    ss_info = data["ss_info"]
    rel_err = data["rel_err"]

    # If abs(rel_err) >= 0.25 add resids and resname to arrays and add them to a file to visualise deviations later on. 
    dev_resids = []
    dev_residnames = []
    dev_ss = []
    for idx in range(len(resids)):
            dev_resids.append((resids[idx]+offset, rel_err[idx]))
            dev_residnames.append(resname[idx])
            dev_ss.append(ss_info[idx])

    resdict[protein] = dev_resids

with open('residue_errors.pkl', 'wb') as f:
     pickle.dump(resdict, f)
