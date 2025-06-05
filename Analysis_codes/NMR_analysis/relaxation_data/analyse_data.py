"""
analyse_data.py
27 April 2025

Analyses N-H Correlation data to classify residues that de-correlate faster than the others
Also analyses differences in the four force fields with respect to residue-wise correlations

NOTE: 
Change the basepath, ff list according to the naming convention you have used.
flag = 1 plots/saves results for all force fields together. flag = 0 only saves/plots results for a single force field.

python3 analyse_data.py <pdbid> <flag> <force field>

"""

import pandas as pd
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import get_all_ss as get_ss

pdb=sys.argv[1]
flag=sys.argv[2]
f = sys.argv[3]
basepath='/home/dynamics/akshatha/'
ff = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']


plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

def return_colours(labels):
    """
    Define a colour scheme for residues depending on their secondary structure
    """
    label_dict = {1:"red", 2:"maroon", 3:"brown", 5: "blue", 4: "grey", 6:"grey", 7:"purple", 8:"pink",  9:"green", 10:"cyan", 11:"violet", 12:"grey"}

    return label_dict[labels]

def plot_corrs(pdb, flag, f):
    # Initialise a dataframe to store the xaxis, yaxis and yerror values
    df = pd.DataFrame()
    basepath='/home/dynamics/akshatha/'
    if int(flag) == 1:
        ff = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    else:
        ff = [f]

    for i,f in enumerate(ff):
        results_path=basepath+f+'/'+pdb+'/Analysis_equi/NMR_analysis/relaxation_data/'

        # Load correlation data 
        data = pd.read_csv(results_path+"NH_Ct_new.csv")

        # compute mean and standard deviation for each column 
        mean_values = data.mean()
        std_dev = data.std()

        xaxis = np.arange(0,len(std_dev)-1)
        yaxis=list(mean_values[1:])
        yerror=list(std_dev[1:])

        # print(xaxis, yaxis, yerror)
        
        
        # Append the xaxis and yaxis values to a dataframe for all force fields
        if int(flag) == 1: 
            if f == 'Final_simfiles':
                df = pd.concat([df, pd.DataFrame({'xaxis': xaxis, 'yaxis_1': yaxis, 'yerror_1': yerror})], axis=1)
            else:
                df = pd.concat([df, pd.DataFrame({'yaxis_'+str(i): yaxis, 'yerror_'+str(i): yerror})], axis=1)

            plt.subplot(2, 1, 1)
            plt.title("N-H Autocorrelation vectors (time-averaged): "+pdb)
            plt.plot(xaxis,yaxis, marker='.', linestyle='-', linewidth=1)     
            plt.grid('both')   
            plt.ylabel("$\mu$")
            
            ss_ranges = get_ss.get_ss_data(pdb)

            plt.fill_between(xaxis, yaxis, color='grey', alpha=0.1)
            for key in ss_ranges:
                for i in range(0, len(ss_ranges[key]), 2):
                    if len(ss_ranges[key]) >= 2:
                        plt.fill_between(xaxis, yaxis, where=(xaxis >= ss_ranges[key][i]) & (xaxis <= ss_ranges[key][i+1]), color=return_colours(key), alpha=0.1)
            
            plt.subplot(2, 1, 2)
            plt.plot(xaxis,yerror, linewidth=1)
            plt.grid('both')
            plt.xlabel("Residue number")
            plt.ylabel("$\sigma$")
            
            df = df.round(4)
            df.to_csv("results/NH_Ct_"+pdb+"_new.csv", index=False)
        else:
            df = pd.concat([df, pd.DataFrame({'ID': xaxis, 'mean':yaxis, 'std':yerror})], axis=1)
            df.to_csv(results_path+"NH_Ct_"+pdb+".csv", index=False)
            print("Results saved to: ", results_path+"NH_Ct_"+pdb+".csv")


        # plt.savefig(pdb+'.png')

plot_corrs(pdb, flag, f)

def get_diff_res(pdb):
    """
    Get the residue numbers that have differences in correlation behaviour across force field models
    """
    data = pd.read_csv("results/NH_Ct_"+pdb+"_new.csv")
    # Get the columns for yaxis and yerror
    yaxis_cols = [col for col in data.columns if 'yaxis' in col]
    yerror_cols = [col for col in data.columns if 'yerror' in col]

    # get the standard deviation of the yaxis values across yaxis_cols and print those that are greater than 0.1
    std_dev = data[yaxis_cols].std(axis=1)
    # highdev = std_dev[std_dev > 0.1].index.tolist()
    # print("Residues with high deviations across models: ", highdev)
    plt.bar(range(len(std_dev)), std_dev, alpha=0.5)
    plt.title("Residues with high deviations across models: "+pdb)
    plt.xlabel("Residue number")
    # plt.ylabel("Standard deviation")
    # give space between two subplots
    plt.tight_layout()

    plt.savefig("results/high_dev_res_"+pdb+"_new.png")
    plt.close()

# GET RESIDUES WHICH SHOW DEVIATIONS ACROSS MODELS, AND THEIR SECONDARY STRUCTURES

if int(flag) == 1:
    get_diff_res(pdb)
