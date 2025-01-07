"""
modularise_data.py 
Akshatha 
24 August 2024

Gets visual representations for secondary structure-wise information for RMSF data

Usage: python3 rmsf_modularise.py 1r36
"""

import pandas as pd
import matplotlib.pyplot as plt
import warnings
import numpy as np
warnings.simplefilter('ignore', DeprecationWarning)
import get_ss_groups as get_ss
import sys
import os

############################################################### Main #####################################################################

pdb = sys.argv[1]
pdbfiles = ['1ozi', '1r36','1z9b', '2ec7', '2gmo', '2juo', '2jxy', '2k0q', '2l4x', '2l6b', '2lht', '2lro', '2m68', '4bwh']



# for pdb in pdbfiles:
filenamecharmm = "/home/dynamics/akshatha/Simfiles/"+pdb+"/Analysis/basic_properties/rmsf_backbone.xvg"
filenameamber = "/home/dynamics/akshatha/Final_simfiles/"+pdb+"/Analysis/basic_properties/rmsf_backbone.xvg"
legends=["CHARMM36m", "Amberff99SB-ILDN"]
colors=['black', 'red']

dummycount=0

# Read data from xmgrace file
for filename in [filenamecharmm, filenameamber]:

    data = pd.read_csv(filename, skiprows=18, sep='\s+', names=['Residue', 'Backbone'])
    offset = data['Residue'].iloc[0]
    data['Residue'] = data['Residue'] - offset + 1

    residues = data['Residue']
    rmsf = data['Backbone']

    # Get secondary structure residues from pdb_ss.txt

    helix_ranges, sheet_ranges = get_ss.get_ss_data(pdb, 0)

    # Plot the RMSF data, with helix and sheet regions highlighted
    plt.rcParams['font.size'] = 13
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']


    plt.plot(residues, rmsf,color=colors[dummycount], linewidth=0.75, label=legends[dummycount])

    for i in range(0, len(helix_ranges), 2):
        if len(helix_ranges) != 0:
            plt.fill_between(residues, rmsf, where=(residues >= helix_ranges[i]) & (residues <= helix_ranges[i+1]), color='red', alpha=0.5)
    for i in range(0, len(sheet_ranges), 2):
        if len(sheet_ranges) != 0:
            plt.fill_between(residues, rmsf, where=(residues >= sheet_ranges[i]) & (residues <= sheet_ranges[i+1]), color='blue', alpha=0.5)

    dummycount += 1
    plt.fill_between(residues, rmsf, color='grey', alpha=0.3)
    plt.xlabel('Residue')
    plt.ylabel('RMSF (nm)')
    plt.legend()
    if pdb not in ["2jxy", "2lro", "2m68", "1r36"]:
        plt.xticks(np.arange(0, residues.iloc[-1], 10))
    plt.xlim(min(residues), max(residues))
    plt.title('Backbone RMSF: '+pdb)

plt.savefig("RMSF_corrected/"+pdb+"_RMSFcomp.png")
plt.close()

