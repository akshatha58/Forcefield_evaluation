"""
map_devs.py
15 November 2024

Maps residues that have an absolute value deviation in chemical shifts
python3 map_devs.py <basepath>
"""
import pandas as pd
import numpy as np
import sys, os
import pickle

basepath=sys.argv[1]
filename="errors_amber.pkl"

sys.path.insert(0, basepath+'/../Analysis_codes/plotting_codes')

import global_plots as gplt

# proteins = ['1ozi', '2gmo', '1z9b', '2juo']
# classes = { 'all_beta: ['1ozi'], 'alpha_beta':['2gmo', '1z9b','2juo']}
# classes = {'all_alpha': ['2l6b', '1r36], 'all_beta: ['2jxy', '2k0q','2m68','1ozi','2lro'], 'alpha_beta':['2gmo', '1z9b', '2l4x','4bwh','2juo'], 'spl':['2lht'], 'few_ss':['2ec7']}
# proteins=['2gmo']

proteins = ['1ozi', '1r36', '2ec7', '2gmo', '2juo', '2jxy', '2k0q', '2l4x', '2l6b', '2lro', '4bwh']


# def return_idx(protein):
#     """
#     Returns offset indices so that the right residue numbers are mapped to the pdb structure file
#     """
#     offset = 0
#     if protein=="1ozi":
#         offset = 1
#     elif protein=="2gmo":
#         offset = -7
#     elif protein=="1z9b":
#         offset = 1
#     elif protein=="2juo":
#         offset = 0
#     return offset
    
resdict={}

for protein in proteins:
    datafile = basepath+protein+"/Analysis/NMR_analysis/chem_shifts/"+protein+"_expshifts_combined.txt"
    directory_path = os.path.dirname(datafile)
    params='H' 
    parameters = gplt.return_paramdata(datafile, params, "shifts")
    data = pd.DataFrame(parameters, columns=["ID", "Res", params+"_exp", params+"_sd_exp", params+"_md", params+"_sd_md", "ss_info", "rel_err_"+params])

    # offset = return_idx(protein)
    # Rename columns to make names easier to read and handle

    resids = data["ID"]
    resname = data["Res"]
    nmr_data = data[params+"_exp"]
    nmr_sd = data[params+"_sd_exp"]
    md_data = data[params+"_md"]
    md_sd = data[params+"_sd_md"]
    ss_info = data["ss_info"]
    rel_err = data["rel_err_"+params]

    # If abs(rel_err) >= 0.25 add resids and resname to arrays and add them to a file to visualise deviations later on. 
    dev_resids = []
    dev_residnames = []
    dev_ss = []
    for idx in range(len(resids)):
            dev_resids.append((resids[idx], rel_err[idx]))
            dev_residnames.append(resname[idx])
            dev_ss.append(ss_info[idx])

    resdict[protein] = dev_resids
    print(resdict)

with open(filename, 'wb') as f:
     pickle.dump(resdict, f)
