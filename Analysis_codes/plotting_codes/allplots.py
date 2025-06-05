"""
newplots.py
11 April 2025

Plots properties across all force field/water model combinations for all proteins 
Plots a consolidated pairwise comparison plot for each NMR observable for all proteins

"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.ndimage import gaussian_filter1d
import pandas as pd
import numpy as np
import sys
import os
import global_plots as gplot

plot_dir = os.path.abspath("../initial_files")
if plot_dir not in sys.path:
    sys.path.append(plot_dir)

import get_all_ss as get_ss

def get_pdb_list(pdblist):
    """
    Gets the list of available pdb paths and corresponding legends according to the simulation time
    Change according to simulations conducted, file paths, and preferred sigma values for Gaussian smoothening
    """

    legends = [os.path.splitext(os.path.basename(pdblist[i]))[0] for i in range(len(pdblist))]

    return pdblist, legends

# Plotting parameters and y-axis labels 
plt.rcParams['font.size'] = 17
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plot_yaxes={'rmsd':'RMSD (nm)', 'rmsf': 'RMSF (nm)', 'rg':'Radius of gyration (nm)', 'potential':'Potential energy (kJ/mol)', 'total_energy': 'Total energy (kJ/mol)', 'kinetic':'Kinetic energy (kJ/mol)','density': 'Density (kg/$m^3$)', 'temperature': 'Temperature (K)', 'pressure': 'Pressure (bar)', 'distance':"Distance (nm)", 'surften':'$\\frac{\Delta G_{solv}}{SASA}$'}

# Helper functions

def get_cols(prop):
    if prop in ['rmsd', 'potential', 'total_energy', 'pressure', 'density', 'temperature', 'kinetic', 'rmsf', 'distance','surften','sasa','solv_free_en']:
        cols = ['xaxis', 'property']
    elif prop == "rg":
        cols= ['xaxis', 'property', 'propx', 'propy', 'propz']
    return cols

def skiprows(prop):
    skip = 0
    if prop in ["potential", "kinetic", "total_energy", "pressure", "temperature", "density", "distance", "surften", "sasa", "solv_free_en"]:
        skip = 24
    return skip

def plot_energies(pdblist, prop):
    """
    Plots energies averaged over the entire system (simulation quality analysis)
    """
    cols = get_cols(prop)
    for pdbnum in range(len(pdblist)):
        if prop not in ['rmsd', 'rg']:
            filename = pdblist[pdbnum]+"/Analysis_equi/basic_properties/"+prop+".xvg"
        else:
            filename = pdblist[pdbnum]+"/Analysis_equi/basic_properties/"+prop+"_all.xvg"

        data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(prop))

        xaxis = data['xaxis']
        yaxis = data['property']

        plt.plot(xaxis, yaxis, linewidth=2, alpha=0.5)
        plt.grid('both')
            
def plot_rmsfs(pdb, prop, pdblist):

    import rmsf_modularise as rfm
    dummycount = 0
    cols = get_cols(prop)
    colours = ['blue', 'orange', 'green', 'red']
    for pdbnum in range(len(pdblist)):
        filename = pdblist[pdbnum]+"/Analysis_equi/basic_properties/"+prop+"_backbone.xvg"
        data = pd.read_csv(filename, sep="\s+", names=cols)
        offset = data['xaxis'].iloc[0]
        data['xaxis'] = data['xaxis'] - offset + 1

        xaxis = data['xaxis']
        yaxis = data['property']

        plt.plot(xaxis, yaxis, linewidth=2, color=colours[dummycount])
        plt.grid('both')

        ss_ranges = get_ss.get_ss_data(pdb)

        plt.plot(xaxis, yaxis)
        for key in ss_ranges:
            for i in range(0, len(ss_ranges[key]), 2):
                if len(ss_ranges[key]) >= 2:
                    plt.fill_between(xaxis, yaxis, where=(xaxis >= ss_ranges[key][i]) & (xaxis <= ss_ranges[key][i+1]), color=rfm.return_colours(key), alpha=0.1)
        
        dummycount += 1
        plt.fill_between(xaxis, yaxis, color='grey', alpha=0.1)
        plt.xlabel('Residues')
        plt.ylim([0,1.2])
        plt.xlim(min(xaxis), max(xaxis))
        plt.title('Backbone RMSF: '+pdb)

def merge_relative_errors(file_list, output_file, param, pdb, proton):
    """
    Merges relative error columns from multiple files based on residue IDs.

    Parameters:
    file_list (list): List of file paths to be merged.
    output_file (str): Path to save the merged output file.
    """
    merged_df = None

    for i, file in enumerate(file_list):
        # Read the file into a DataFrame
        if param == "chem_shifts":
            df = pd.read_csv(file, sep="\s+", usecols=["ID", "Res","rel_err_"+proton, "ss_info"], dtype={"ID": int, "Res":str, "rel_err_"+proton: float, "ss_info": str})
            df.rename(columns={"rel_err_"+proton: f"rel_err_file{i+1}"}, inplace=True)
            df.rename(columns={"Res": f"Res{i+1}"}, inplace=True)
            df.rename(columns={"ss_info": f"ss{i+1}"}, inplace=True)

        elif param == "J3_constants":
            df = pd.read_csv(file, sep="\s+", usecols=["ID","Res", "rel_err", "ss_info"], dtype={"ID": int, "Res": str, "rel_err": float, "ss_info": str})
            df.rename(columns={"rel_err": f"rel_err_file{i+1}"}, inplace=True)
            df.rename(columns={"Res": f"Res{i+1}"}, inplace=True)
            df.rename(columns={"ss_info": f"ss{i+1}"}, inplace=True)
       
        # For each line, if the relative error is greater than max, add the protein, ID, residue name, and prop to another file
        outfile = 'pairwise_diffs/pairwise_deviations_full_'+str(i+1)+'_'+proton+'.txt'
        with open(outfile, 'w') as f:
            for index, row in df.iterrows():
                if abs(row[f'rel_err_file{i+1}']) > maxerr(param,proton):
                    f.write(f"{pdb} \t{row['ID']} \t{row[f'Res{i+1}']} \t{row[f'ss{i+1}']}\t{param}\n")
        f.close()

        # Merge with the existing DataFrame
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on="ID", how="outer")
            # Fill NaN values with 0
            merged_df.fillna(-100, inplace=True)

    # Sort by ID for better readability
    merged_df.sort_values(by="ID", inplace=True)

    

    # Save the merged DataFrame to a file
    merged_df.to_csv(output_file, index=False, sep="\t")
    print(f"Merged file saved to {output_file}")

def return_ss(filename):
    df = pd.read_csv(filename, sep="\s+", usecols=["ss_info"], dtype={"ss_info": str})
    # convert ss_info to a list of integers
    ss_info = list(map(int, df['ss_info'].tolist()))
    return ss_info

def maxerr(prop, proton):
    if prop == "chem_shifts":
        if proton == 'H':
            return 0.06
        elif proton == 'C':
            return 0.03
        elif proton == 'N':
            return 0.04
    elif prop == "J3_constants":
        return 0.25

def plot_consolidated_nmr(prop, pdb):

    # Get relative error data across all models 
    ffs = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']

    thing = {'chem_shifts':'expshifts_', 'J3_constants':'J'}
    thing2 = {'chem_shifts':'shifts', 'J3_constants':'J3'}

    params = {'chem_shifts':['H','C','N'], 'J3_constants': 'NA'}

    files = []
        

    for ff in ffs:
        filename = '/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/NMR_analysis/'+prop+'/'+pdb+'_'+thing[prop]+'combined.txt'
        ss_info = return_ss(filename)
        # print(ss_info)
        files.append(filename)

    if prop == 'J3_constants':
        output = 'pairwise_diffs/relerr_combined_'+pdb+'_J.txt'
        param='NA'
        axislabel = "J3 constants"
        merge_relative_errors(files, output, prop, pdb, param)

        # Get pairwise differences between absolute values of relative errors
        merged_data = pd.read_csv(output, sep='\s+')

        # Same force field, different water models
        # remove rows with -800 values in the pair of columns for each comparison
        indices = [['rel_err_file1', 'rel_err_file3'], ['rel_err_file2', 'rel_err_file4'], ['rel_err_file1', 'rel_err_file4'], ['rel_err_file2', 'rel_err_file3'], ['rel_err_file1', 'rel_err_file2'], ['rel_err_file3', 'rel_err_file4']]
        err = []

        plt.subplots(3, 2, figsize=(15, 8))
        plt.suptitle('Differences in error values across models for '+pdb+': '+axislabel, fontsize=20)
        colours = gplot.return_colours(ss_info)

        # Plot them as follows
        # 1-3     2-4
        # 1-4     2-3
        # 1-2     3-4

        titles = ['Amber (TIP4P-Ew versus TIP3P)', 'CHARMM (TIP3Pm versus TIP4P-Ew)', 'TIP4P-Ew (Amber versus CHARMM)', 'TIP3P(m) (CHARMM versus Amber)', 'Amber/TIP4P-Ew versus CHARMM/TIP3Pm', 'Amber/TIP3P versus CHARMM/TIP4P-Ew']

        for index, index_pair in enumerate(indices):
            merged_data = merged_data[(merged_data[index_pair[0]] > -100) & (merged_data[index_pair[1]] > -100)]
            # print(merged_data)
            resids = merged_data["ID"].tolist()
            data1 = merged_data[index_pair[0]].tolist()
            data2 = merged_data[index_pair[1]].tolist()

            err = np.ones(len(data1))

            for j in range(len(data1)):
                if data1[j]*data2[j] > 0:
                    # same sign
                    err[j] = abs(data1[j]) - abs(data2[j])
                else:
                    err[j] = abs(data1[j]) + abs(data2[j])
                    if abs(data1[j]) < abs(data2[j]):
                        err[j]*=-1

            plt.subplot(3, 2, index+1)
            plt.bar(resids, err, color=colours)
            plt.title(titles[index],fontsize=18)

            if index+1 in [1, 3, 5]:
                plt.ylabel('$\Delta$', fontsize=18)
            plt.grid(axis='both', alpha=0.5)
            if index+1 == 5 or index+1 == 6:
                plt.xlabel('Residues',fontsize=18)

            plt.ylim([-1,1])


        plt.subplots_adjust(hspace=0.55)
        plt.savefig('pairwise_diffs/'+prop+'_'+pdb+'.png')

    if prop == 'chem_shifts':
       
        for param in params[prop]:
            output = 'pairwise_diffs/relerr_combined_'+pdb+'_'+param+'.txt'
            axislabel = param+" Chemical shifts"
            merge_relative_errors(files, output, prop, pdb, param)

            # Get pairwise differences between absolute values of relative errors
            merged_data = pd.read_csv(output, sep='\s+')

            # Same force field, different water models
            # remove rows with -800 values in the pair of columns for each comparison
            indices = [['rel_err_file1', 'rel_err_file3'], ['rel_err_file2', 'rel_err_file4'], ['rel_err_file1', 'rel_err_file4'], ['rel_err_file2', 'rel_err_file3'], ['rel_err_file1', 'rel_err_file2'], ['rel_err_file3', 'rel_err_file4']]
            err = []

            plt.subplots(3, 2, figsize=(15, 8))
            plt.suptitle('Differences in error values across models for '+pdb+': '+axislabel, fontsize=20)
            colours = gplot.return_colours(ss_info)

            # Plot them as follows
            # 1-3     2-4
            # 1-4     2-3
            # 1-2     3-4

            titles = ['Amber (TIP4P-Ew versus TIP3P)', 'CHARMM (TIP3Pm versus TIP4P-Ew)', 'TIP4P-Ew (Amber versus CHARMM)', 'TIP3P(m) (CHARMM versus Amber)', 'Amber/TIP4P-Ew versus CHARMM/TIP3Pm', 'Amber/TIP3P versus CHARMM/TIP4P-Ew']

            for index, index_pair in enumerate(indices):
                merged_data = merged_data[(merged_data[index_pair[0]] > -100) & (merged_data[index_pair[1]] > -100)]
                # print(merged_data)
                resids = merged_data["ID"].tolist()
                data1 = merged_data[index_pair[0]].tolist()
                data2 = merged_data[index_pair[1]].tolist()

                err = np.ones(len(data1))

                for j in range(len(data1)):
                    if data1[j]*data2[j] > 0:
                        # same sign
                        err[j] = abs(data1[j]) - abs(data2[j])
                    else:
                        err[j] = abs(data1[j]) + abs(data2[j])
                        if abs(data1[j]) < abs(data2[j]):
                            err[j]*=-1

                plt.subplot(3, 2, index+1)
                plt.bar(resids, err, color=colours)
                plt.title(titles[index],fontsize=18)

                if index+1 in [1, 3, 5]:
                    plt.ylabel('$\Delta$', fontsize=18)
                plt.grid(axis='both', alpha=0.5)
                if index+1 == 5 or index+1 == 6:
                    plt.xlabel('Residues',fontsize=18)

                if param == 'C' or param == 'N':
                    plt.ylim([-0.25,0.25])
                    plt.yticks(np.arange(-0.25, 0.26, 0.1))
                elif param == 'H':
                    plt.ylim([-0.5,0.5])
                    plt.yticks(np.arange(-0.5, 0.56, 0.2))


            plt.subplots_adjust(hspace=0.55)
            plt.savefig('pairwise_diffs/'+prop+'_'+pdb+'_'+param+'.png')

def get_res_avgs(param):
    atomtype = {'H':'H', 'C':'CA', 'N':'N'}
    
    filename = "/home/dynamics/akshatha/mmcif_files/mdshifts_4.txt"
    data = pd.read_csv(filename, sep='\s+', names=['num', 'Res', 'N', 'H', 'C', 'N_sd', 'H_sd', 'C_sd', 'coverage'])
    amino_acid_dict = {"ALA": [0, 0, 0],"ARG": [0, 0,0],"ASN": [0, 0,0],"ASP": [0, 0,0],"CYS": [0, 0,0],"GLN": [0, 0,0],"GLU": [0, 0,0],"GLY": [0, 0,0],"HIS": [0, 0, 0],"ILE": [0, 0, 0],"LEU": [0, 0, 0],"LYS": [0, 0, 0],"MET": [0, 0, 0],"PHE": [0, 0, 0],"PRO": [0, 0, 0],"SER": [0, 0, 0],"THR": [0, 0, 0],"TRP": [0, 0, 0],"TYR": [0, 0, 0],"VAL": [0, 0, 0]}

    grouped_data = data[(data[param] != -1) & (data['coverage'] > 0.9)].groupby('Res')
    for prot, group in grouped_data:
        amino_acid_dict[prot][0] = group[param].mean()
        amino_acid_dict[prot][1] = group[param+'_sd'].mean()
        amino_acid_dict[prot][2] += len(group[param])

    bmrb_file = "/home/dynamics/akshatha/NMR_BMRB/shifts/chemshift_stats.csv"
    # comp_id,atom_id,count,min,max,avg,std,num_outliers
    bmrb_data = pd.read_csv(bmrb_file, sep=",", usecols=["comp_id", "atom_id", "count","avg", "std"], dtype={"comp_id": str, "atom_id": str, "count":int, "avg": float, "std": float})

    with open("/home/dynamics/akshatha/mmcif_files/mdshifts_avg_" + param + "_4.txt", 'w') as f:
        f.write("Res\tMD_Avg\tMD_Std_Dev\tMD_Count\tBMRB_Avg\tBMRB_Std_Dev\tBMRB_Count\n")
        for key in amino_acid_dict:
            bmrb_avg = bmrb_data[(bmrb_data["comp_id"] == key) & (bmrb_data["atom_id"] == atomtype[param])]["avg"].mean()
            bmrb_std_dev = bmrb_data[(bmrb_data["comp_id"] == key) & (bmrb_data["atom_id"] == atomtype[param])]["std"].mean()
            bmrb_count = bmrb_data[(bmrb_data["comp_id"] == key) & (bmrb_data["atom_id"] == atomtype[param])]["count"].mean()
            f.write(f"{key}\t{amino_acid_dict[key][0]:.4f}\t{amino_acid_dict[key][1]:.4f}\t{amino_acid_dict[key][2]}\t{bmrb_avg:.4f}\t{bmrb_std_dev:.4f}\t{int(bmrb_count)}\n")
            # f.write(f"{key}\t{bmrb_avg:.4f}\t{bmrb_std_dev:.4f}\t{int(bmrb_count)}\n")
    f.close()

def plot_err_bandwidths(prop, pdb):

    # Get relative error data across all models 
    ffs = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    params = ['H','C','N']
    indices = {'Final_simfiles':'rel_err_file1', 'Simfiles':'rel_err_file2', 'amber_tip3p':'rel_err_file3', 'charmm_tip4p':'rel_err_file4'}
    
    for param in params:
        mergedfile = '/home/dynamics/akshatha/Analysis_codes/plotting_codes/pairwise_diffs/relerr_combined_'+pdb+'_'+param+'.txt'
        # Get pairwise differences between absolute values of relative errors
        merged_data = pd.read_csv(mergedfile, sep='\s+')
        axislabel = param+" Chemical shifts"
        index = 0
        titles = ['Amber/TIP4P-Ew', 'CHARMM/TIP3Pm', 'Amber/TIP3P', 'CHARMM/TIP4P-Ew']


        plt.subplots(4, 1, figsize=(12, 8))
        plt.suptitle('Relative errors for '+pdb+': '+axislabel, fontsize=20)
        

        for ff in ffs:
            index += 1
            filename = '/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/NMR_analysis/'+prop+'/'+pdb+'_expshifts_combined.txt'
            ss_info = return_ss(filename)
            colours = gplot.return_colours(ss_info)
            merged_data = merged_data[(merged_data[indices[ff]] > -100)]
            # print(merged_data)
            resids = merged_data["ID"].tolist()
            data = merged_data[indices[ff]].tolist()

            plt.subplot(4, 1, index)
            plt.bar(resids, data, color=colours)
            # Bandwidths
            plt.axhspan(-maxerr(prop, param), maxerr(prop, param), color='grey', alpha=0.3) 
            plt.title(titles[index-1],fontsize=15)

            plt.grid(axis='both', alpha=0.5)

            if index == 4:
                plt.xlabel('Residues',fontsize=20)

            if param == 'H':
                plt.ylim([-0.5,0.5])
            else:
                plt.ylim([-0.25, 0.25])
            # plt.yticks(np.arange(-0.5, 0.52, 0.25))

            plt.subplots_adjust(hspace=0.55)
            plt.gcf().text(0.04, 0.5, 'Relative Errors', va='center', rotation='vertical', fontsize=20)
            plt.savefig('pairwise_diffs/'+prop+'_'+pdb+'_'+param+'_bands.png')

def calc_surften(pdb, ff):
    fig, ax = plt.subplots(3, 1, figsize=(17, 15))
    plt.legend(ff, fontsize=18)
    for f in ff:
        basefile = '/home/dynamics/akshatha/'+f+'/'+pdb+'/'+'Analysis_equi/basic_properties/'
        sasa = basefile + 'sasa.xvg'
        surf_en = basefile + 'solv_free_en.xvg'

        cols = get_cols('sasa')
        sasa_data = pd.read_csv(sasa, sep="\s+", names=cols, skiprows=skiprows('sasa'))
        cols = get_cols('solv_free_en')
        surf_en_data = pd.read_csv(surf_en, sep="\s+", names=cols, skiprows=skiprows('solv_free_en'))
        surf_ten = abs(surf_en_data['property']/sasa_data['property'])
        print(f)
        print(np.mean(surf_ten), np.std(surf_ten))

        plt.subplot(3,1, 1)
        # calculate the moving average of sasa_data
        # sasa_data['property'] = sasa_data['property'].rolling(window=200).mean()
        plt.plot(sasa_data['xaxis'], sasa_data['property'], linewidth=2)
        plt.ylabel('SASA ($nm^2$)')
        plt.grid('both')
        plt.xticks(np.arange(100,1100,100))

        plt.subplot(3, 1, 2)
        # surf_en_data['property'] = surf_en_data['property'].rolling(window=200).mean()
        plt.plot(surf_en_data['xaxis'], surf_en_data['property'], linewidth=2)
        plt.ylabel('$\Delta G_{solvation}$')
        plt.grid('both')
        plt.xticks(np.arange(100,1100,100))

        plt.subplot(3, 1, 3)
        # surf_ten = surf_ten.rolling(window=200).mean()
        plt.plot(sasa_data['xaxis'], surf_ten, linewidth=2)
        plt.ylabel('$\gamma = $' + plot_yaxes['surften'])
        plt.xlabel('Time (ns)')
        plt.grid('both')
        fig.suptitle('Surface calculations: '+pdb, fontsize=30)
        plt.xticks(np.arange(100,1100,100))

        # plt.xlabel('Time (ns)', fontsize=15)
        # plt.ylabel(plot_yaxes['surften'], fontsize=15)
    
    plt.savefig('surften_'+pdb+'.png')

def plot_nonrmsf(prop):
    basepath='/home/dynamics/akshatha/'
    ff=['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    plotfile=prop+".png"
    pdb = ['1ozi', '2gmo', '2juo', '2k0q', '2l4x', '1r36','2m68', '4bwh']

    i=0
    fig, ax = plt.subplots(len(pdb), 1, figsize=(17, 15))
    for p in pdb:
        i+=1
        plt.subplot(len(pdb), 1, i)
        plotlist=[]
        fflist = []
        for f in ff:
            base=basepath+f+"/"
            if os.path.exists(base+p+"/Analysis_equi/"):
                plotlist.append(base+p)
                fflist.append(f)
        plot_energies(plotlist, prop)
        # calc_surften(p, fflist)
        plt.legend([p], fontsize=18, loc='right')

    plt.subplots_adjust(hspace=0.55)
    # plt.tight_layout()
    plt.xlabel('Time (ns)', fontsize=28)
    fig.text(0.02, 0.5, plot_yaxes[prop], va='center', rotation='vertical', fontsize=28)
    plt.savefig(plotfile)

def plot_rmsf():
    basepath='/home/dynamics/akshatha/'
    ff=['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    prop = 'rmsf'
    plotfile=prop+".png"
    pdb = ['1ozi', '2gmo', '2juo', '2k0q', '2l4x', '2m68', '4bwh', '1r36']
    i=0

    fig, ax = plt.subplots(2, 4, figsize=(16, 10))
    for p in pdb:
        i+=1
        plt.subplot(2, 4, i)
        plotlist=[]
        for f in ff:
            base=basepath+f+"/"
            if os.path.exists(base+p+"/Analysis_equi/"):
                plotlist.append(base+p)
                print(p, f)
        plot_rmsfs(p, 'rmsf', plotlist)
        if i in [1, 5]:
            plt.ylabel('RMSF (nm)')
        # plt.legend([p], fontsize=18, loc='right')

    plt.subplots_adjust(hspace=0.55)
    plt.tight_layout()
    # plt.xlabel('Time (ps)', fontsize=28)
    # fig.text(0.05, 0.5, plot_yaxes[prop], va='center', rotation='vertical', fontsize=28)
    plt.savefig(plotfile)

def plot_most_deviating():
    plt.figure(figsize=(10, 6))
    params = ['H', 'C', 'N']

    for param in params:
        colour = {1:'blue', 2:'orange', 3:'green', 4:'red'}
        for idx in range(1,3,2):
            data = pd.read_csv('pairwise_diffs/pairwise_deviations_full_'+str(idx)+'_'+param+'.txt', sep='\s+', names=['protein', 'resnum', 'resname', 'ss','prop'])
            # data = pd.read_csv('pairwise_diffs/pairwise_deviations.txt', sep='\s+', names=['protein', 'resnum', 'resname', 'ss','prop'])
            plt.subplot(1, 2, 1)
            merged_data = data[(data['prop'] == 'chem_shifts')]
            bars = plt.bar(
                merged_data['ss'].value_counts().index,  # Residue names
                merged_data['ss'].value_counts().values,  # Counts
                color=colour[idx],  # Bar color
                edgecolor='black',  # Bar edge color
                width=0.8,  # Reduce bar width for separation
                alpha = 0.3

            )

            # Add gridlines
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.xticks(np.arange(0,13,1))

            # Add labels and title
            plt.xlabel('Secondary structures', fontsize=14)
            plt.ylabel('Counts', fontsize=14)
            # plt.title('Highest Deviations )', fontsize=16)
            # Customize ticks
            plt.xticks(fontsize=8, rotation=45)
            plt.yticks(fontsize=12)

            plt.subplot(1, 2, 2)
            merged_data_s = data[(data['prop'] == 'chem_shifts')]
            bars = plt.bar(
                merged_data_s['resname'].value_counts().index,  # Residue names
                merged_data_s['resname'].value_counts().values,  # Counts
                color=colour[idx],  # Bar color
                edgecolor='black',  # Bar edge color
                width=0.8,  # Reduce bar width for separation
                alpha=0.3
            )

            # Add gridlines
            plt.grid(axis='y', linestyle='--', alpha=0.7)

            # Add labels and title
            plt.xlabel('Residues', fontsize=14)
            plt.ylabel('Counts', fontsize=14)
            # plt.title('Highest Pairwise Deviations (H-Shifts)', fontsize=16)

            # Customize ticks
            plt.xticks(fontsize=10, rotation=45)
            plt.yticks(fontsize=12)
            

        # Save the plot
        plt.tight_layout()
        plt.savefig('pairwise_diffs/residues_hs'+param+'.png')

def plot_res_avgs():
    # plot bar plots of residue averages and error bars as standard deviations
    j = 1
    plt.figure(figsize=(10, 8))

    ylim = {'H':[4,10], 'C':[40,70], 'N':[90,145]}
    for param in ['H', 'C', 'N']:
        for i in range(1, 5):
            file = '/home/dynamics/akshatha/mmcif_files/mdshifts_avg_'+param+'_'+str(i)+'.txt'
            data = pd.read_csv(file, sep='\s+')
            plt.subplot(3, 1, j)
            if i == 1:
                plt.bar(data['Res'], data['BMRB_Avg'], yerr=data['BMRB_Std_Dev'], color='green', capsize=7, edgecolor='black', width=0.8, alpha = 0.3, label='BMRB')
            # plt.bar(data['Res'], data['MD_Avg'], yerr=data['MD_Std_Dev'], capsize=3, edgecolor='black', width=0.8, alpha = 0.3, label='MD')
            plt.scatter(data['Res'], data['MD_Avg'], s=20, label=i, marker="o")
            plt.errorbar(data['Res'], data['MD_Avg'], yerr=data['MD_Std_Dev'], alpha=0.5, capsize=5)
        plt.gcf().text(0.02, 0.5, 'Chemical shifts (ppm)', va='center', rotation='vertical', fontsize=18)
        plt.xticks(fontsize=13, rotation=45)
        plt.ylim(ylim[param])
        plt.grid('x')
        plt.ylabel(param)
        j+=1
        # plt.legend()

    plt.xlabel('Residue Name', fontsize=18)
    plt.subplots_adjust(hspace=0.55)
    plt.savefig('pairwise_diffs/BMRB_avgs_all.png')

def plot_res_sasa(pdb):
    """
    Plots the average residue-wise SASA against the reference values in the NMR structure
    """
    import rmsf_modularise as rfm
    plt.figure(figsize=(10, 8))
    # ffs = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    ffs=['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    label={'Final_simfiles':'Amber/4P-Ew', 'Simfiles':'CHARMM/3Pm', 'amber_tip3p':'Amber/3P', 'charmm_tip4p':'CHARMM/4P-Ew'}
    for ff in ffs:
        
        basepath='/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/basic_properties/'
        file=basepath+'res_sasa.xvg'
        data = pd.read_csv(file, sep='\s+', skiprows=25, names=['res', 'avg', 'std'])
        offset = data['res'].iloc[0]
        data['res'] = data['res'] - offset + 1

        plt.subplot(2,1,1)
        plt.plot(data['res'], data['avg'], label=label[ff])
        plt.xticks(fontsize=13, rotation=45)
        plt.grid('both')
        plt.ylabel("SASA Average")

        ss_ranges = get_ss.get_ss_data(pdb)
        plt.fill_between(data['res'], data['avg'], color='grey', alpha=0.1)
        for key in ss_ranges:
            for i in range(0, len(ss_ranges[key]), 2):
                if len(ss_ranges[key]) >= 2:
                    plt.fill_between(data['res'], data['avg'], where=(data['res'] >= ss_ranges[key][i]) & (data['res'] <= ss_ranges[key][i+1]), color=rfm.return_colours(key), alpha=0.1)
        plt.legend()
        plt.subplot(2, 1, 2)
        # plotting the difference between reference and simulation SASAs
        ref_file=basepath+'res_sasa_ref.xvg'
        ref_data = pd.read_csv(ref_file, sep='\s+', skiprows=25, names=['res', 'avg', 'std'])
        offset = ref_data['res'].iloc[0]
        ref_data['res'] = ref_data['res'] - offset + 1
        diff = ref_data['avg'] - data['avg']
        plt.plot(ref_data['res'], diff)
        plt.grid('both')
        plt.xticks(fontsize=13, rotation=45)
        plt.ylabel("$\delta$ SASA (Ref - MD)")
        plt.xlabel("Residue")


    plt.savefig('sasa/'+pdb+'_avg.png')

def plot_2gmodists():
    ffs = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    for ff in ffs:
        file='/home/dynamics/akshatha/'+ff+'/2gmo/Analysis_equi/basic_properties/2gmo_distances.xvg'
        data = pd.read_csv(file, sep='\s+', skiprows=24, names=['time', 'dist'])
        plt.plot(data['time'], data['dist'], label=ff)

    plt.xlabel('Time (ps)')
    plt.ylabel('Distance (nm)')
    plt.title('2GMO distances')
    plt.legend()
    plt.savefig('2gmo_distances.png')


def rg_distributions():
    """
    Plots Rg distributions for all proteins
    """
    pdb= ['1ozi', '2gmo', '2juo', '2k0q', '2l4x', '2m68', '4bwh', '1r36', '2lro','2lht','2jxy','2l6b']
    ffs = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    plt.figure(figsize=(10, 8))

    for p in pdb:
        plt.subplot(3, 4, pdb.index(p)+1)
        for ff in ffs:
            file = '/home/dynamics/akshatha/'+ff+'/'+p+'/Analysis_equi/basic_properties/rg_all.xvg'
            if os.path.exists(file):
                data = pd.read_csv(file, sep='\s+', skiprows=24, names=['time', 'xaxis', 'property', 'propx', 'propy', 'propz'])
            
                # plot a histogram of the Rg values
                plt.hist(data['property'], bins=30, alpha=0.4, density=True)

                if pdb.index(p) + 1 in [1, 5, 9]:
                    plt.ylabel('Density')
                if pdb.index(p) +1 in [9, 10, 11, 12]:
                    plt.xlabel('Rg (nm)')
                
                plt.title(p)
                plt.xlim([0.5,2])
                plt.ylim([0,7.5])
                plt.yticks(fontsize=14)
                plt.xticks(np.arange(0.5, 2.1, 0.5), fontsize=14)
                plt.grid('both')
                # plt.legend()
    
        plt.tight_layout()
    plt.savefig('rg_distributions.png')


rg_distributions()










# plot_consolidated_nmr('chem_shifts', '2l4x')
# plot_consolidated_nmr('J3_constants', '2gmo')
# plot_consolidated_nmr('chem_shifts', '2k0q')

# for prop in ['density','rg', 'rmsd', 'total_energy']:
    # plot_nonrmsf(prop)

# plot_rmsf()
# plot_res_sasa('4bwh')

# for param in ['H', 'C', 'N']:
    # get_res_avgs(param)

# plot_res_avgs()



