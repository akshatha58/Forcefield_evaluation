"""
summary_plots.py
15 December 2024
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os

num_resids = {'1r36':140, '2l6b':115,'1ozi':99,'2lro':172,'2jxy':194,'2k0q':74,'2m68':142,'2gmo':83,'1z9b':135,'2l4x':128,'4bwh':88,'2juo':89,'2lht':122}

pdb_listmid = ['4bwh', '2l6b', '2k0q', '2gmo']
pdb_listhigh = ['1z9b', '2ec7', '2lht', '2l4x']
pdb_listlow = ['1ozi', '1r36', '2juo', '2jxy', '2lro', '2m68']

allproteins = pdb_listmid + pdb_listhigh + pdb_listlow
pdb_both_ff = ['1ozi', '1r36', '2gmo', '2k0q']

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
       
def plot_ffprop(ff_file, pdb, prop, cols):
    if os.path.isfile(ff_file):
        data_amber = pd.read_csv(ff_file, sep="\s+", names=cols, skiprows=18)
        time = data_amber['time']
        yaxis_amber = data_amber['property']
        print(pdb, " Amber")
        print("Average"+prop+": ",np.mean(yaxis_amber))
        print("Std Dev"+prop+": ",np.std(yaxis_amber))
        if prop == 'rmsd':
            rmsd_100 = np.mean(yaxis_amber)/(1+np.log(num_resids[pdb]/100))
            print("Normalised RMSD_100: ",rmsd_100)
        plt.plot(time, yaxis_amber, linewidth=1.5, label=pdb)

def plot(pdblist, prop, plotfile):
    """
    Plots structural properties for proteins
    """
    iter = 0
    fig, ax = plt.subplots(len(pdblist), 1, figsize=(17, 15))
    if prop == "rmsd":
        # plt.suptitle("RMSD (Backbone)", fontsize=30)
        cols = ['time', 'property']
    elif prop == "rg":
        # plt.suptitle("Radius of Gyration", fontsize=30)
        cols= ['time', 'property', 'propx', 'propy', 'propz']
                 
    for pdb in pdblist:
        colour = colorfold(pdb)
        iter += 1

        # Extracting property information from Amber force field
        # for ss in ['helix', 'sheet', 'flexible']:
        filename_amber = dir+"/Final_simfiles/"+pdb+"/Analysis_equi/basic_properties/"+prop+".xvg"
        plt.subplot(len(pdblist), 1, iter)
        plot_ffprop(filename_amber, pdb, prop, cols)
        

        if pdb in pdb_both_ff:
            # For proteins where CHARMM simulations are complete
            # for ss in ['helix', 'sheet', 'flexible']:
            filename_charmm = dir+"/Simfiles/"+pdb+"/Analysis_equi/basic_properties/"+prop+".xvg"
            plot_ffprop(filename_charmm, pdb, prop, cols)
            filename_charmm_tip4p = dir+"/charmm_tip4p"+pdb+"/Analysis_equi/basic_properties/"+prop+".xvg"
            plot_ffprop(filename_charmm_tip4p, pdb, prop, cols)
            
        plt.grid('y')
        plt.legend(fontsize=25)

        if iter != len(pdblist):
            plt.xticks([])
        
    # plt.xlabel("Time (ns)", fontsize=20)
    fig.text(0.5, 0.04, 'Time (ns)', ha='center', fontsize=30)
    if prop == "rmsd":
        fig.text(0.04, 0.5, 'RMSD (nm)', va='center', rotation='vertical', fontsize=30)
    elif prop == "rg":
        fig.text(0.04, 0.5, 'Rg (nm)', va='center', rotation='vertical', fontsize=30)
    plt.savefig("plots/"+plotfile)

def annotation(prop, data_amber, data_charmm, colour):
    # for i, label in enumerate(data_charmm['Protein']):
        # plt.annotate(label, (data_amber[prop+'_mean'].iloc[i], data_charmm[prop+'_mean'].iloc[i]), fontsize=5)
    plt.errorbar(data_amber[prop+'_mean'], data_charmm[prop+'_mean'], yerr=data_charmm[prop+'_sd'], xerr=data_amber[prop+'_sd'], alpha=0.5, color='violet', linestyle='')
    plt.scatter(data_amber[prop+'_mean'], data_charmm[prop+'_mean'], color=[colour, colour, colour, colour],s=20)
   
def plot_averages(plotfile):
    """
    Plots structural properties for proteins
    """
    cols = ['Protein',	'helix_mean',	'sheet_mean',	'flexible_mean',	'rmsd_mean', 'helix_sd',	'sheet_sd',	'flexible_sd',	'rmsd_sd',	'rg_mean',	'rg_sd']

    filename_amber = dir+"/Analysis_codes/plotting_codes/plots/rgrmsd_amber.txt"
    amber = pd.read_csv(filename_amber, sep="\s+")
    filename_charmm = dir+"/Analysis_codes/plotting_codes/plots/rgrmsd_charmm.txt"
    data_charmm = pd.read_csv(filename_charmm, sep="\s+")
    data_amber = amber[amber.iloc[:, 0].isin(data_charmm['Protein'])]

    # annotation('helix', data_amber, data_charmm, 'red')
    # annotation('sheet', data_amber, data_charmm, 'blue')
    # annotation('flexible', data_amber, data_charmm, 'black')
    # annotation('rmsd', data_amber, data_charmm, 'green')
    annotation('rg', data_amber, data_charmm, 'black')

    refline = np.arange(1,1.6,0.001)
    plt.plot(refline, refline, "--")
    # plt.xlim([0,1.4])
    # plt.ylim([0,1.4])
    plt.ylabel("Rg (nm) CHARMM")
    plt.xlabel("Rg (nm) Amber")
    plt.grid('both')
    plt.savefig("plots/"+plotfile)

plot(['1ozi'], "rmsd", "rmsd_1ozi.png")
# plot(pdb_listmid, "rmsd", "rmsd_mid.png")
# plot(pdb_listlow, "rmsd", "rmsd_low.png")
# plot(pdb_listhigh, "rg", "rg_high.png")
# plot(pdb_listmid, "rg", "rg_mid.png")
# plot(pdb_listlow, "rg", "rg_low.png")


# plot_averages('rg_all.png')
