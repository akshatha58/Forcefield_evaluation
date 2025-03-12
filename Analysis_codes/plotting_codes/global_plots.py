"""
global_plots.py
Akshatha
15 October 2024
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.patches as mpatches
import sys
import os

# Specify complete file path as input
# shift_data = sys.argv[1]

# Trialfiles
# datafile = "1ozi_expshifts_combined.txt"

# Specify general plot styles

plt.rcParams['font.size'] = 18
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# Helper functions related to data importing and extraction

def load_data(datafile):
    """
    Loads the datafile and returns its header and data separately
    Input: path of datafile (string)
    Output: header (list), data (2D array)
    """
    # print("Loading data...")
    d = pd.read_csv(datafile, sep="\s+")
    header = d.columns.tolist()
    data = d.values
    return header, data

def return_paramdata(datafile, params, prop):
    """
    Returns selected_params, a 2D array of the averaged NMR and MD values and secondary structures of the available residues of a particular isotope
    """
    if prop=="shifts":
        cols_to_select=["ID", "Res", params+"_exp", params+"_sd_exp", params+"_md", params+"_sd_md", "ss_info", "rel_err_"+params]
    if prop=="J3":
        cols_to_select=["ID", "Res", "J_exp", "J_sd_exp","J_md", "J_sd_md", "ss_info", "rel_err"]

    header, data = load_data(datafile)

    header_idx = []
    for cols in cols_to_select:
        header_idx.append(header.index(cols))
        
    selected_params = data[:,header_idx]

    if prop=="shifts":
        boolean_mask = selected_params[:,2] >= 0
        selected_params = selected_params[boolean_mask]

    return selected_params

# Helper functions related to plotting data

def return_absolute_shiftrange(params):
    """
    Returns fixed y axis ranges according to whether H, C or N shifts are calculated
    """
    if params=="N":
        return [90, 150]
    elif params=="C":
        return [40, 70]
    elif params=="H":
        return [2, 12] 

def return_colours(labels):
    """
    Define a colour scheme for residues depending on their secondary structure
    """
    colours=[]
    label_dict = {1:"red", 2:"maroon", 3:"brown", 5: "blue", 4: "grey", 6:"grey", 7:"grey", 8:"pink",  9:"green", 10:"cyan", 11:"violet", 12:"black"}

    for label in labels:
        colours.append(label_dict[label])

    return colours

# Main plotting function

def plot_shift_results(datafile, params):
    """
    # Subplot 1 (top): bar plot of relative errors wrt residues
    # Subplot 2 (bottom): Line/scatter plot of MD and NMR shifts (NMR in dotted lines)
    """
    # Extract directory name for saving files
    directory_path = os.path.dirname(datafile)
    savefile=directory_path+"/"+params+"_chemshifts.png"

    # Get parameter-specific data
    parameters = return_paramdata(datafile, params, "shifts")
    data = pd.DataFrame(parameters, columns=["ID", "res", params+"_exp", params+"_sd_exp", params+"_md", params+"_sd_md", "ss_info", "rel_err_"+params])

    # Rename columns to make names easier to read and handle
    resids = data["ID"]
    nmr_data = data[params+"_exp"]
    nmr_sd = data[params+"_sd_exp"]
    md_data = data[params+"_md"]
    md_sd = data[params+"_sd_md"]
    ss_info = data["ss_info"]
    rel_err = data["rel_err_"+params]
    # print(rel_err)

    # Make subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8))
    fig.suptitle('NMR versus MD: '+params+" chemical shift values", fontsize=20)
    colours = return_colours(ss_info)

    # Bar plot of relative errors (residue-wise)
    plt.subplot(2, 1, 1)
    plt.bar(resids, rel_err, color=colours)
    plt.grid(axis='y', alpha=0.5)
    plt.ylabel("Relative error", fontsize=18)
    plt.ylim([-0.5, 0.5])

    # Line plot of chemshift values (NMR and MD)
    plt.subplot(2, 1, 2)
    plt.scatter(resids, nmr_data, marker="o", linestyle="--", c=colours, label="NMR shifts")
    plt.errorbar(resids, nmr_data, yerr=nmr_sd, alpha=0.5, color='green')
    plt.scatter(resids, md_data,  marker="^", c=colours, label="MD shifts")
    plt.errorbar(resids, md_data, yerr=md_sd, alpha=0.5, color='grey')
    plt.grid(axis='y', alpha=0.5)
    plt.ylabel(params+" shifts (ppm)", fontsize=18)
    plt.xlabel("Residues", fontsize=18)
    plt.ylim(return_absolute_shiftrange(params))
    plt.legend(ncol=2, loc="upper right")
    

    # Save the plot after adding secondary structure labels in the legend

    legend_handles = [
        mpatches.Patch(color='red', label='Alpha Helix'),
        mpatches.Patch(color="maroon", label="3-10 Helix"),
        mpatches.Patch(color="brown", label="Pi Helix"),
        mpatches.Patch(color="grey", label="Turns/bends"),
        mpatches.Patch(color='blue', label='Sheets'),
        mpatches.Patch(color="green", label="Alpha Termini"),
        mpatches.Patch(color="cyan", label="3-10 Termini"),
        mpatches.Patch(color="violet", label="Pi Termini"),
        mpatches.Patch(color='pink', label='Kappa Helix'),
        mpatches.Patch(color='black', label='Flexible'),
    ]
    
    # fig.legend(handles=legend_handles, ncol=4, loc='upper center', bbox_to_anchor=(0.5, 0.95), fontsize=14)
    plt.savefig(savefile)
    plt.close()

def plot_J3_results(datafile):
    """
    # Subplot 1 (top): bar plot of relative errors wrt residues
    # Subplot 2 (bottom): Line/scatter plot of MD and NMR shifts (NMR in dotted lines)
    """
    # Extract directory name for saving files
    directory_path = os.path.dirname(datafile)
    savefile=directory_path+"/J3.png"

    # ID   Res   J_exp  J_sd_exp  ID  J_md  J_sd_md  ss_info
    params='NA' # dummy variable
    parameters = return_paramdata(datafile, params, "J3")
    data = pd.DataFrame(parameters, columns=["ID", "Res", "J_exp", "J_sd_exp", "J_md", "J_sd_md", "ss_info", "rel_err"])

    # Rename columns to make names easier to read and handle
    resids = data["ID"]
    nmr_data = data["J_exp"]
    nmr_sd = data["J_sd_exp"]
    md_data = data["J_md"]
    md_sd = data["J_sd_md"]
    ss_info = data["ss_info"]
    rel_err = data["rel_err"]        

    # Make subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8))
    fig.suptitle('NMR versus MD: J3 (HN-HA coupling constant values)', fontsize=20)
    colours = return_colours(ss_info)

    # Bar plot of relative errors (residue-wise)
    plt.subplot(2, 1, 1)
    plt.bar(resids, rel_err, color=colours)
    plt.grid(axis='y', alpha=0.5)
    plt.ylabel("Relative error", fontsize=18)
    plt.ylim([-1, 1])

    # Line plot of chemshift values (NMR and MD)
    plt.subplot(2, 1, 2)
    plt.scatter(resids, nmr_data, marker="o", linestyle="--", c=colours, label="NMR results")
    plt.errorbar(resids, nmr_data, yerr=nmr_sd, alpha=0.5, color='green')
    plt.scatter(resids, md_data,  marker="^", c=colours, label="MD results")
    plt.errorbar(resids, md_data, yerr=md_sd, alpha=0.5, color='grey')
    plt.grid(axis='y', alpha=0.5)
    plt.ylabel("J3 coupling constants (Hz)", fontsize=18)
    plt.xlabel("Residues", fontsize=18)
    plt.legend(ncol=2, loc="upper right")
    # plt.ylim([1,14])
    

    # Save the plot after adding secondary structure labels in the legend

    legend_handles = [
        mpatches.Patch(color='red', label='Alpha Helix'),
        mpatches.Patch(color="maroon", label="3-10 Helix"),
        mpatches.Patch(color="brown", label="Pi Helix"),
        mpatches.Patch(color="grey", label="Turns/bends"),
        mpatches.Patch(color='blue', label='Sheets'),
        mpatches.Patch(color="green", label="Alpha Termini"),
        mpatches.Patch(color="cyan", label="3-10 Termini"),
        mpatches.Patch(color="violet", label="Pi Termini"),
        mpatches.Patch(color='pink', label='Kappa Helix'),
        mpatches.Patch(color='black', label='Flexible'),
    ]
    
    # fig.legend(handles=legend_handles, ncol=4, loc='upper center', bbox_to_anchor=(0.5, 0.95), fontsize=14)
    plt.savefig(savefile)
    plt.close()

def plot_shifts_scatter(datafile, params):

    """
    Plots relative errors between MD and NMR chemical shifts.
    param_dict = {params+"_exp": expdata, params+"_sd_exp": exp_sd, params+"_md": mddata, params+"_sd_md": md_sd}
    """
    parameters = return_paramdata(datafile, params, "shifts")
    data = pd.DataFrame(parameters, columns=["ID", "res", params+"_exp", params+"_sd_exp", params+"_md", params+"_sd_md", "ss_info", "rel_err_"+params])

    # Rename columns to make names easier to read and handle
    nmr_data = data[params+"_exp"]
    md_data = data[params+"_md"]
    ss_info = data["ss_info"]

    colours = return_colours(ss_info)
    plt.scatter(nmr_data, md_data, color=colours)
    refline = np.arange(min(min(md_data), min(nmr_data)), max(max(md_data), max(nmr_data)))
    plt.plot(refline, refline, "--")

    plt.grid(axis='y', alpha=0.5)
    plt.ylabel(params+" shifts-MD (ppm)", fontsize=18)
    plt.xlabel(params+" shifts-NMR (ppm)", fontsize=18)
    # plt.xticks(np.arange(0,110,10))
    # plt.yticks(np.arange(20,110,10))
    return nmr_data, md_data

def plot_J3_scatter(datafile):
    """
    Plots a scatter plot of J3 coupling constant data
    """
    legend_handles = [
        mpatches.Patch(color='red', label='Alpha Helix'),
        mpatches.Patch(color="maroon", label="3-10 Helix"),
        mpatches.Patch(color="brown", label="Pi Helix"),
        mpatches.Patch(color="grey", label="Turns/bends"),
        mpatches.Patch(color='blue', label='Sheets'),
        mpatches.Patch(color="green", label="Alpha Termini"),
        mpatches.Patch(color="cyan", label="3-10 Termini"),
        mpatches.Patch(color="violet", label="Pi Termini"),
        mpatches.Patch(color='pink', label='Kappa Helix'),
        mpatches.Patch(color='black', label='Flexible'),
    ]

    params='NA' # dummy variable
    parameters = return_paramdata(datafile, params, "J3")
    data = pd.DataFrame(parameters, columns=["ID", "Res", "J_exp", "J_sd_exp", "J_md", "J_sd_md", "ss_info", "rel_err"])

    # Rename columns to make names easier to read and handle
    resids = data["ID"]
    nmr_data = data["J_exp"]
    nmr_sd = data["J_sd_exp"]
    md_data = data["J_md"]
    md_sd = data["J_sd_md"]
    ss_info = data["ss_info"]

    figure = plt.figure(figsize=(10,8))
    colours = return_colours(ss_info)
    plt.scatter(nmr_data, md_data, color=colours)
    refline = np.arange(min(min(md_data), min(nmr_data)), max(max(md_data), max(nmr_data)))
    plt.plot(refline, refline, "--")

    plt.grid(axis='y', alpha=0.5)
    plt.ylabel("J3-MD (Hz)", fontsize=18)
    plt.xlabel("J3-NMR (Hz)", fontsize=18)
    # plt.xticks(np.arange(0,110,10))
    # plt.yticks(np.arange(20,110,10))

    figure.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=4, fontsize=14)
    figure.suptitle("J3 Coupling Constant values (HN - HA): MD versus NMR", fontsize=20)

    directory_path = os.path.dirname(datafile)
    savefile=directory_path+"/J3_comp.png"
    figure.savefig(savefile)
    plt.close(figure)

def consolidate_plots(datafile, what2plot):
    """
    Consolidates data from 1H, 13C and 15N chemical shifts into a single three-part plot
    what2plot toggles between "separate", "errors" and "shifts"

    separate: plot shifts versus residue information (NMR and MD separately) with standard deviations
    errors: plot relative errors with residue numbers (bar plot)
    shifts: plot NMR versus MD shifts (scatter plot)

    Consolidates data from H, C, and N shifts in one plot
    """
    directory_path = os.path.dirname(datafile)

    legend_handles = [
        mpatches.Patch(color='red', label='Alpha Helix'),
        mpatches.Patch(color="maroon", label="3-10 Helix"),
        mpatches.Patch(color="brown", label="Pi Helix"),
        mpatches.Patch(color="grey", label="Turns/bends"),
        mpatches.Patch(color='blue', label='Sheets'),
        mpatches.Patch(color="green", label="Alpha Termini"),
        mpatches.Patch(color="cyan", label="3-10 Termini"),
        mpatches.Patch(color="violet", label="Pi Termini"),
        mpatches.Patch(color='pink', label='Kappa Helix'),
        mpatches.Patch(color='black', label='Flexible'),
    ]

    params = ['H', 'C', 'N']
    subplot_num = 0
    
    if what2plot=="scatter":
        figure = plt.figure(figsize=(18,8))

        for param in params:
            subplot_num+=1
            plt.subplot(1, 3, subplot_num)
            plot_shifts_scatter(datafile, param)
            
        figure.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=4, fontsize=14)
        figure.suptitle("Chemical shift values: 1H, 13C, 15N: MD versus NMR", fontsize=20)

        directory_path = os.path.dirname(datafile)
        savefile=directory_path+"/chemshifts_comp.png"
        figure.savefig(savefile)
        plt.close(figure)

    elif what2plot == "all":
        for param in params:
            plot_shift_results(datafile, param)
        
# Get statistics from data

def get_r2(md_data, nmr_data):
    """
    Function that calculates r-squared value (agreement of md_data with a reference nmr data)
    """
    md_mean = np.mean(md_data)
    sst = np.sum((md_data - md_mean)**2)
    exp_pred = nmr_data
    ssr = np.sum((md_data - exp_pred)**2)
    r_squared = 1 - (ssr/sst)
    return r_squared

def stats(datafile, params, whichss, prop):
    """
    Get the chi squared values for the data
    """
    # Get parameter-specific data 
    if prop == "shifts":
        parameters = return_paramdata(datafile, params, 'shifts')
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
        parameters = return_paramdata(datafile, params, "J3")
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


    chi_sq = 0
    avg_rel_err = 0 
    # r2 = 0
    mean_nmr = 0
    mean_md = 0

    if len(nmr_data) != 0:
        for idx in range(len(nmr_data)):
            chi_sq += ((nmr_data[idx] - md_data[idx])**2)/md_sd[idx]**2
            mean_nmr += nmr_data[idx]
            mean_md += md_data[idx]

        chi_sq = chi_sq/len(nmr_data)
        avg_rel_err = np.mean(rel_err)
        mean_nmr = mean_nmr/len(nmr_data)
        mean_md = mean_md/len(md_data)

        rmsd = np.sqrt(np.mean((nmr_data - md_data) ** 2))
        # r2 = get_r2(md_data, nmr_data)
    
    else:
        chi_sq = -800
        rmsd = -800
        # r2 = -800
  
    # return chi_sq, rmsd, avg_rel_err, r2, mean_nmr, mean_md
    return chi_sq, rmsd, avg_rel_err, mean_nmr, mean_md





