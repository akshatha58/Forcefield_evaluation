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

# pdb = sys.argv[1]
# prop = sys.argv[2]
pdbfiles = ['1ozi', '1r36','1z9b', '2gmo', '2juo', '2jxy', '2k0q', '2l4x', '2l6b', '2lht', '2lro', '2m68', '4bwh']

def get_files(foldernames, pdb, prop):
    """
    Get the files to be plotted in plot_prop
    """
    files = []
    for ff in foldernames:
        print(ff)
        if prop == "rmsd" or prop == "rg":
            files.append("../../"+ff+"/"+pdb+"/Analysis/basic_properties/"+prop+"_all.xvg")
        elif prop == "rmsf":
            files.append("/../../"+ff+"/"+pdb+"/Analysis/basic_properties/"+prop+"_backbone.xvg")

    return files

def plot_rmsd(pdb, foldernames):
    # Read data from xmgrace file
    dummycount=0
    files = get_files(foldernames, pdb, "rmsd")

    colors=['black', 'red', 'green', 'darkviolet']
    legends=["CHARMM/3P", "Amber/4P", "CHARMM/4P", "Amber/3P"]
    # legends = ["Amber/4P", "Amber/3P"]
    # colors=['red', 'darkviolet']

    for filename in files:

        data = pd.read_csv(filename, sep='\s+', names=['Residue', 'Backbone'])
        offset = data['Residue'].iloc[0]
        data['Residue'] = data['Residue'] - offset + 1

        residues = data['Residue']
        rmsd = data['Backbone']

        # Get secondary structure residues from pdb_ss.txt

        helix_ranges, sheet_ranges = get_ss.get_ss_data(pdb, 0)

        # Plot the RMSF data, with helix and sheet regions highlighted
        plt.rcParams['font.size'] = 13
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']


        plt.plot(residues, rmsd,color=colors[dummycount], linewidth=0.75, label=legends[dummycount])
        dummycount += 1

        plt.xlabel('Time (ns)')
        plt.ylabel('RMSD (nm)')
        plt.legend()
        plt.xlim([0,500])
        # if pdb not in ["2jxy", "2lro", "2m68", "1r36"]:
        #     plt.xticks(np.arange(0, residues.iloc[-1], 10))
        # plt.xlim(min(residues), max(residues))
        plt.title('Backbone RMSD: '+pdb)

    plt.savefig("comparisons/"+pdb+"_RMSDcomp.png")
    plt.close()

def plot_rmsf(pdb, foldernames):

    # Read data from xmgrace file
    dummycount=0
    files = get_files(foldernames, pdb, "rmsf")

    colors=['black', 'red', 'green', 'darkviolet']
    legends=["CHARMM/3P", "Amber/4P", "CHARMM/4P", "Amber/3P"]
    # legends = ["Amber/4P", "Amber/3P"]
    # colors=['red', 'darkviolet']


    for filename in files:

        data = pd.read_csv(filename, sep='\s+', names=['Residue', 'Backbone'])
        offset = data['Residue'].iloc[0]
        data['Residue'] = data['Residue'] - offset + 1

        residues = data['Residue']
        rmsf = data['Backbone']

        # Get secondary structure residues from pdb_ss.txt
        print(pdb)
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
        plt.xlabel('Residues')
        plt.ylabel('RMSF (nm)')
        plt.xlim([0,500])
        plt.legend()
        if pdb not in ["2jxy", "2lro", "2m68", "1r36"]:
            plt.xticks(np.arange(0, residues.iloc[-1], 10))
        plt.xlim(min(residues), max(residues))
        plt.title('Backbone RMSF: '+pdb)

    plt.savefig("comparisons/"+pdb+"_RMSFcomp.png")
    plt.close()

        # print(pdb)
        # filenameamber = "/home/dynamics/akshatha/Final_simfiles/"+pdb+"/Analysis_equi/basic_properties/"+prop+"_backbone.xvg"
        # colors=['black']
        # legends=["Amberff99SB-ILDN"]

        # for filename in [filenameamber]:

        #     data = pd.read_csv(filename, sep='\s+', names=['Residue', 'Backbone'], skiprows=17)
        #     offset = data['Residue'].iloc[0]
        #     data['Residue'] = data['Residue'] - offset + 1

        #     residues = data['Residue']
        #     rmsf = data['Backbone']

        #     # Get secondary structure residues from pdb_ss.txt

        #     helix_ranges, sheet_ranges = get_ss.get_ss_data(pdb, 0)

        #     # Plot the RMSF data, with helix and sheet regions highlighted
        #     plt.rcParams['font.size'] = 13
        #     plt.rcParams['font.family'] = 'serif'
        #     plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']


        #     plt.plot(residues, rmsf,color=colors[dummycount], linewidth=0.75, label=legends[dummycount])

        #     for i in range(0, len(helix_ranges), 2):
        #         if len(helix_ranges) != 0:
        #             plt.fill_between(residues, rmsf, where=(residues >= helix_ranges[i]) & (residues <= helix_ranges[i+1]), color='red', alpha=0.5)
        #     for i in range(0, len(sheet_ranges), 2):
        #         if len(sheet_ranges) != 0:
        #             plt.fill_between(residues, rmsf, where=(residues >= sheet_ranges[i]) & (residues <= sheet_ranges[i+1]), color='blue', alpha=0.5)

        #     dummycount += 1
        #     plt.fill_between(residues, rmsf, color='blue', alpha=0.3)
        #     plt.xlabel('Residues')
        #     plt.ylabel('RMSF (nm)')
        #     plt.xlim([0,500])
        #     plt.legend()
        #     if pdb not in ["2jxy", "2lro", "2m68", "1r36"]:
        #         plt.xticks(np.arange(0, residues.iloc[-1], 10))
        #     plt.xlim(min(residues), max(residues))
        #     plt.title('Backbone RMSF: '+pdb)

        # plt.savefig("newplots/comparisons/"+pdb+"_RMSFcomp.png")
        # plt.close()

def plot_rg(pdb, foldernames):
    dummycount=0
    files = get_files(foldernames, pdb, "rg")
    
    colors=['black', 'red', 'green', 'darkviolet']
    legends=["CHARMM/3P", "Amber/4P", "CHARMM/4P", "Amber/3P"]
    # legends = ["Amber/4P", "Amber/3P"]
    # colors=['red', 'darkviolet']

    for filename in files:

        data = pd.read_csv(filename, sep='\s+', names=['Time', 'Backbone', 'x', 'y', 'z'])

        time = data['Time']
        rmsd = data['Backbone']

        # Get secondary structure residues from pdb_ss.txt
        # helix_ranges, sheet_ranges = get_ss.get_ss_data(pdb, 0)

        # Plot the RMSF data, with helix and sheet regions highlighted
        plt.rcParams['font.size'] = 13
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

        plt.plot(time, rmsd,color=colors[dummycount], linewidth=0.75, label=legends[dummycount])

        dummycount += 1
        # plt.fill_between(time, rmsd, color='grey', alpha=0.3)
        plt.xlabel('Time')
        plt.ylabel('Rg (nm)')
        plt.xlim([0,500])
        plt.legend()
        # if pdb not in ["2jxy", "2lro", "2m68", "1r36"]:
        #     plt.xticks(np.arange(0, residues.iloc[-1], 10))
        # plt.xlim(min(residues), max(residues))
        # plt.ylim([1.3,1.7])
        plt.title('Rg: '+pdb)

    plt.savefig("comparisons/"+pdb+"_Rgcomp.png")
    plt.close()


pdb=["2k0q"]
# "CHARMM/3P", "Amber/4P", "CHARMM/4P", "Amber/3P"
foldernames = ["Simfiles","Final_simfiles", "charmm_tip4p", "amber_tip3p"]

for prot in pdb:
    plot_rmsd(prot, foldernames)
    plot_rmsf(prot, foldernames)
    plot_rg(prot, foldernames)
