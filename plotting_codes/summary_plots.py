"""
summary_plots.py
15 December 2024
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys


pdb_listmid = ['4bwh', '2l6b', '2k0q']
pdb_listhigh = ['1z9b', '2ec7', '2lht', '2l4x']
pdb_listlow = ['1ozi', '1r36', '2gmo', '2juo', '2jxy', '2lro', '2m68']
dir = sys.argv[1]

plt.rcParams['font.size'] = 25
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

def colorfold(pdb):
    """
    Assigns plot colour based on the class of the protein
    """
    classdict = {'a':['1r36', '2l6b'], 'b':['1ozi', '2lro', '2jxy', '2k0q', '2m68'], 'c':['2gmo', '1z9b', '2l4x', '4bwh', '2juo'], 'd':['2lht'], 'e':['2ec7']}
    if pdb in classdict['a']:
        return 'red'
    elif pdb in classdict['b']:
        return 'green'
    elif pdb in classdict['c']:
        return 'black'
    elif pdb in classdict['d']:
        return 'blue'
    else:
        return 'orange'
    
    
def plot(pdblist, prop, plotfile):
    """
    Plots structural properties for proteins
    """
    iter = 0
    fig, ax = plt.subplots(len(pdblist), 1, figsize=(15, 15))
    if prop == "rmsd":
        # plt.suptitle("RMSD (Backbone)", fontsize=30)
        cols = ['time', 'property']
    elif prop == "rg":
        # plt.suptitle("Radius of Gyration", fontsize=30)
        cols= ['time', 'property', 'propx', 'propy', 'propz']
                 
    for pdb in pdblist:
        colour = colorfold(pdb)
        iter += 1
        filename = dir+pdb+"/Analysis/basic_properties/"+prop+".xvg"
        data = pd.read_csv(filename, sep="\s+", names=cols)
        time = data['time']
        yaxis = data['property']
        # xorient = data['propx']
        # yorient = data['propy']
        # zorient = data['propz']

        plt.subplot(len(pdblist), 1, iter)
        plt.plot(time, yaxis, linewidth=1.5, label=pdb, color=colour)
        # plt.plot(time, xorient, label='x')
        # plt.plot(time, yorient, label='y')
        # plt.plot(time, zorient, label='z')
        # plt.ylabel("RMSD (nm)", fontsize=16)
        plt.grid('both')
        plt.legend(fontsize=25, loc="lower right")
        
    # plt.xlabel("Time (ns)", fontsize=20)
    fig.text(0.5, 0.04, 'Time (ns)', ha='center', fontsize=30)
    if prop == "rmsd":
        fig.text(0.04, 0.5, 'RMSD (nm)', va='center', rotation='vertical', fontsize=30)
    elif prop == "rg":
        fig.text(0.04, 0.5, 'Rg (nm)', va='center', rotation='vertical', fontsize=30)
    plt.savefig("plots/"+plotfile)

    


# plot(pdb_listhigh, "rmsd", "rmsd_high.png")
# plot(pdb_listmid, "rmsd", "rmsd_mid.png")
# plot(pdb_listlow, "rmsd", "rmsd_low.png")
# plot(pdb_listhigh, "rg", "rg_high.png")
# plot(pdb_listmid, "rg", "rg_mid.png")
# plot(pdb_listlow, "rg", "rg_low.png")

