"""
hierarchical_analysis.py
16 May 2025

Analyses hierarchical structures in data, from residue level to secondary structure level to protein level.

Features (residue level):
- RMSF
- SASA_avg and SASA_std
- Molecular weight
- Secondary structure type
- Residue type
- NOE, T1, T2 (Experimental data)
- J3 constant relative errors (if available)
- 1H, 15N, 13C Chemical shift relative errors
- N-H Bond correlations (avg and std)
- Configurational entropy

Workflows:
- Residue-level analysis: heatmap of normalised residue-level data
- Pairwise correlation analysis: heatmap of pairwise correlations between residue-level data (cosine similarity)
    - For each protein
    - Cosine similarity averaged over all proteins
"""

import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
import mdtraj as md

import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

plot_dir = os.path.abspath("../initial_files")
if plot_dir not in sys.path:
    sys.path.append(plot_dir)

import get_all_ss as ssgroups

plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

MWT = {'ALA':0.089, 'ARG':0.174, 'ASN':0.132, 'ASP':0.133, 'CYS':0.121, 'GLU':0.147, 'GLN':0.146,
              'GLY':0.075, 'HIS':0.155, 'ILE':0.131, 'LEU':0.131, 'LYS':0.146, 'MET':0.149,
              'PHE':0.165, 'PRO':0.115, 'SER':0.105, 'THR':0.119, 'TRP':0.204, 'TYR':0.181,
              'VAL':0.117}

def return_filename(prop, pdb):
    """
    Returns the filename where each feature is stored. 
    """
    prop_dict = {'RMSF': '/Analysis_equi/basic_properties/rmsf_backbone.xvg', 
                 'SASA': '/Analysis_equi/basic_properties/res_sasa.xvg',
                 'NH_Corr':'/Analysis_equi/NMR_analysis/relaxation_data/NH_Ct_'+pdb+'.csv', 
                 'Shifts': '/Analysis_equi/NMR_analysis/chem_shifts/'+pdb+'_expshifts_combined.txt',
                 'J3': '/Analysis_equi/NMR_analysis/J3_constants/'+pdb+'_Jcombined.txt',
                 'NOE': '/Analysis_equi/NMR_analysis/relaxation_data/'+pdb+'_NOEexp.txt', 
                 'T1': '/Analysis_equi/NMR_analysis/relaxation_data/'+pdb+'_T1exp.txt',
                 'T2': '/Analysis_equi/NMR_analysis/relaxation_data/'+pdb+'_T2exp.txt'}
    skip = 0
    # skiprows = {'RMSF': 0, 'SASA': 25, 'NH_Corr': 0, 'Shifts': 1, 'J3': 1, 'Config_entropy':1}
    if prop == "SASA":
        skip = 25
    elif prop in ['Shifts', 'J3']:
        skip = 1
    elif prop == "NH_Corr":
        skip = 2
    
    return prop_dict[prop], skip

def get_residue_order(pdb, directory):
    """
    Gets the residue number and residue names from the .gro file.
    """
    # Get the residue number/residue name from the .gro file 
    gro_file = os.path.join(directory, 'Analysis_equi', 'processed_trajs', f'{pdb}_nowater.gro')
    resids = {}

    # if gro_file does not exist, run the following shell command:
    # printf "1 1" | /usr/local/gromacs/bin/gmx trjconv -f $mdmod_gro -s $mdmod_gro -o $struct_nowater -pbc whole -center
    if not os.path.exists(gro_file):
        import subprocess
        # change directory
        subprocess.run(['./get_nowater.sh', directory], check=True, cwd=os.path.dirname(__file__))
    top = md.load(gro_file)
    # Extract the residue names and their corresponding indices
    resnames = [top.topology.residue(i).name for i in range(top.topology.n_residues)]
    resids = {i+1: resnames[i] for i in range(len(resnames))}

    return resids

def get_data(directory, prop, pdb, resids):
    
    """
    If prop is Shifts, it extracts rel_err_H, rel_err_N, and rel_err_C, ss_info, and res from the file. 
    If prop is J3, it extracts rel_err and res from the file.
    If prop is NH_Corr, it extracts res, mean, sd from the file.
    If prop is RMSF, it extracts res (starting from 1) and rmsf from the file.
    If prop is SASA, it extracts res, sasa, sasa_sd from the file.
    """
    # Map residue names to residue numbers in the protein 
    filename, skip = return_filename(prop, pdb)
    
    # If filename does not exist, calculate the property
    if os.path.exists(directory+filename) == False:
        import subprocess
        subprocess.run(['./compute_property.sh', prop, directory], check=True, cwd=os.path.dirname(__file__))

    if prop == 'Shifts':
        data = pd.read_csv(directory+filename, sep='\s+')
        data = data[['ID','Res', 'rel_err_H', 'rel_err_N', 'rel_err_C', 'C_exp', 'N_exp', 'H_exp', 'C_md', 'N_md', 'H_md']]

    elif prop == 'J3':
       offset = {'1ozi':1, '2gmo':-7, '2juo':0}
       data = pd.read_csv(directory+filename, sep='\s+')
       data = data[['ID', 'Res', 'rel_err']]
       # add an offset to the ID column 
       data['ID'] = data['ID'] + offset[pdb]

    elif prop == 'NH_Corr':
       data = pd.read_csv(directory+filename, sep=',', skiprows=skip, names=['ID', 'mean_NH', 'std_NH'])
       data['ID'] = data['ID'] - data['ID'].min() + 1
       # Map residue numbers to residue names
       data['Res'] = data['ID'].map(resids)

    elif prop == 'RMSF':
       data = pd.read_csv(directory+filename, sep='\s+', names=['ID', 'rmsf'], skiprows=skip)
       # renumber 'res' to start from 1: subtract the first value of 'res' from all values
       data['ID'] = data['ID'] - data['ID'].min() + 1
       # Map residue numbers to residue names
       data['Res'] = data['ID'].map(resids)

    elif prop == 'SASA':
       data = pd.read_csv(directory+filename, sep='\s+', names=['ID', 'sasa', 'sasa_sd'], skiprows=skip)
         # renumber 'res' to start from 1: subtract the first value of 'res' from all values
       data['ID'] = data['ID'] - data['ID'].min() + 1
            # Map residue numbers to residue names
       data['Res'] = data['ID'].map(resids)

    elif prop == 'NOE':
        data = pd.read_csv(directory+filename, sep="\s+", skiprows=skip, names=['ID', 'Res', 'NOE', 'std_NOE'])
    
    elif prop == 'T1':
        data = pd.read_csv(directory+filename, sep="\s+", skiprows=skip, names=['ID', 'Res', 'T1', 'std_T1'])

    elif prop == 'T2':
        data = pd.read_csv(directory+filename, sep="\s+", skiprows=skip, names=['ID', 'Res', 'T2', 'std_T2'])
    
    if pdb == '2gmo' and prop in ['NOE', 'T1', 'T2']:
        # subtract 8 from the ID column
        data['ID'] = data['ID'] - 8
        
    elif prop == 'Config_entropy':
        data = pd.read_csv(directory+filename, sep=',', skiprows=skip, names=['Res', 'ID', 'entropy', 'std_Ent'])

    return data

def get_raw_data(ff, pdb, features):
    # initialise empty dataframe
    compiled_df = pd.DataFrame()

    # if pdb in ['1ozi', '2gmo', '2juo'] and 'J3' not in features:
    #     features.append('J3')
    # elif pdb in ['2k0q', '2l4x', '2m68', '4bwh'] and 'J3' in features:
    #     features.remove('J3')

    directory = os.path.join(os.path.dirname(__file__), '..', '..', ff, pdb)
    if not os.path.exists(directory+'/Analysis_equi'):
        print(f"Directory {directory} does not exist.")
    
    print("Getting data for", pdb, "and force field", ff)
    resids = get_residue_order(pdb, directory)
    

    for i, prop in enumerate(features):
        print(f"Getting {prop} data")
        data = get_data(directory, prop, pdb, resids)
        if i == 0:
            data['MWT'] = data['Res'].map(MWT)
            ss_data = os.path.join(directory, '..','..','mmcif_files', f'{pdb}_sslists.txt')
            ss_dict = ssgroups.get_ss_residues(ss_data)

            # map ss_info to the corresponding residue number
            data['ss_info'] = 12
            for i in range(len(data)):
                for key in ss_dict.keys():
                    if data.loc[i,'ID'] - 1 in ss_dict[key]:
                        data.loc[i, 'ss_info'] = key
                        break
        
        print(data)
        # Concatenate dataframes with same ID and Res
        if compiled_df.empty:
            compiled_df = data
        else:
            # Merge the dataframes on 'ID' and 'Res'
            # Use outer join to keep all data
            data['ID'] = data['ID'].astype(int)
            data['Res'] = data['Res'].astype(str)
            compiled_df = pd.merge(compiled_df, data, on=['ID', 'Res'], how='outer')
        
    print(compiled_df)

    # # fill empty values with -100
    # compiled_df.fillna(-100, inplace=True)
    # # Convert ss_info to int
    # if 'ss_info' in compiled_df.columns:
        # compiled_df['ss_info'] = compiled_df['ss_info'].astype(int)

    compiled_df.drop(columns=['ss_info'], inplace=True)
    compiled_df.drop(columns=['Res'], inplace=True)
    compiled_filename=  os.path.join('results', ff, f'{pdb}_compiled_raw.txt')
    compiled_df.to_csv(compiled_filename, sep='\t', index=True)
    
def normalise_columns(file, outfile):
    """
    Normalises the columns of the dataframe.
    """
    df = pd.read_csv(file, sep='\t')
    # Drop the first column (ID)
    df.drop(columns=['Unnamed: 0'], inplace=True)
    # Normalise all columns except the first one
    df.iloc[:, 1:] = (df.iloc[:, 1:] - df.iloc[:, 1:].min()) / (df.iloc[:, 1:].max() - df.iloc[:, 1:].min())
    df.to_csv(outfile, sep='\t', index=False)

def save_heatmap(file, pdb, ff):
    data = pd.read_csv(file, sep='\t')

    # Plot heatmap with xaxis as ID and yaxis as features
    data.drop(columns=['ID'], inplace=True)
    # data.drop(columns=['Unnamed: 0'], inplace=True)

    # Replace columns with -100 values with NaN
    data.replace(-100, np.nan, inplace=True)
    sns.heatmap(data.corr(), cmap='coolwarm', square=True, linecolor='black', linewidths=0.5, cbar_kws={"shrink": 0.8})
    # data.drop(columns=['ID'], inplace=True)

    corr = data.corr()
    if pdb == '2l4x':
        corr.replace(np.nan, 0, inplace=True)
    corr.to_csv('results/heatmaps/new/corr_matrix_'+pdb+'_'+ff+'.csv', sep='\t')

    plt.title(pdb + ' ' + ff)

def plot_heatmap(data):

    # Plot heatmap with xaxis as ID and yaxis as features
    data.replace(-100, np.nan, inplace=True)

    sns.heatmap(data, cmap='coolwarm', square=True, linecolor='black', linewidths=0.5, cbar_kws={"shrink": 0.8})
       
    plt.xticks(ticks=np.arange(len(data.transpose().columns))+0.5, labels=data.transpose().columns, rotation=90)
  
def main():
    # DATASETS
    ffs = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    pdbs = ['1ozi','2gmo','2juo','2k0q', '2l4x', '2m68', '4bwh']
    features = ['RMSF', 'SASA', 'NOE', 'T1', 'T2', 'Shifts', 'NH_Corr']

    for pdb in pdbs:
        plt.figure(figsize=(14,14))
        # for i,ff in enumerate(ffs):
        ff = 'amber_tip3p'
        print("Analysing correlations for "+pdb+" and force field "+ff)
        get_raw_data(ff, pdb, features)
        normalise_columns('results/'+ff+'/'+pdb+'_compiled_raw.txt', 'results/'+ff+'/'+pdb+'_norm.txt')
        print("======================================================================================")

        # plt.subplot(2,2,i+1)
        save_heatmap('results/'+ff+'/'+pdb+'_norm.txt', pdb, ff)
        # plt.title(pdb + ' ' + ff)
    
        # leave space between subplots
        # plt.subplots_adjust(hspace=0.3, wspace=0.3)
        # plt.suptitle('Heatmap of Normalised Residue-Level Data Correlations for '+pdb, fontsize=18)
        # plt.subplots_adjust(top=0.9)
        plt.savefig('results/heatmaps/new/'+pdb+'_heatmap.png')
        plt.close()

def get_pairwise_differences(ff):
    """
    Get the average and std dev of the pairwise differences in the correlation matrix 
    """
    pdb = ['1ozi','2gmo','2juo','2k0q', '2m68', '4bwh', '2l4x']
    

    # Create all pairwise combinations of pdbs
    pairs = []
    for i in range(len(pdb)):
        for j in range(i+1, len(pdb)):
            pairs.append((pdb[i], pdb[j]))
    
    # Create a dataframe to store pairwise differences as a list for each element in the matrix
    # Initialise with zeros
    pairwise_differences = pd.DataFrame()
    pairwise_diff_sq = pd.DataFrame()

    for pair in pairs:
        corr_1 = pd.read_csv('results/heatmaps/new/corr_matrix_'+pair[0]+'_'+ff+'.csv', sep='\t')
        corr_2 = pd.read_csv('results/heatmaps/new/corr_matrix_'+pair[1]+'_'+ff+'.csv', sep='\t')

        # Ignore the first column and first row
        corr_1 = corr_1.iloc[:, 1:]
        corr_2 = corr_2.iloc[:, 1:]

        # if pair_diff is not empty, append the new pair_diff to the existing pair_diff
        if pairwise_differences.empty:
            pairwise_differences = np.abs(corr_1 - corr_2)

            pairwise_diff_sq = (corr_1 - corr_2)**2
        else:
            pairwise_differences += np.abs(corr_1 - corr_2)
            pairwise_diff_sq += (corr_1 - corr_2)**2

    print(pairwise_differences)

    pairwise_differences = pairwise_differences / len(pairs)
    pairwise_diff_sq = pairwise_diff_sq / len(pairs)
    pairwise_diff_std = np.sqrt(pairwise_diff_sq - pairwise_differences**2)

    pairwise_differences.to_csv('results/heatmaps/new/pairwise_differences_'+ff+'_avg.csv', sep='\t')
    pairwise_diff_std.to_csv('results/heatmaps/new/pairwise_differences_'+ff+'_std.csv', sep='\t')

    plt.figure(figsize=(12,12))
    plot_heatmap(pairwise_differences.T)
    
    plt.savefig('results/heatmaps/new/pairwise_differences_'+ff+'_avg.png')
    plt.close()
    
    plt.figure(figsize=(12,12))
    plot_heatmap(pairwise_diff_std.T)
    plt.savefig('results/heatmaps/new/pairwise_differences_'+ff+'_std.png')
    plt.close()
        
def get_pairwise_differences_overall():
    """
    Get the average and std dev of the pairwise differences in the correlation matrix for ALL FORCE FIELDS
    """
    ff = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']

    # Create all pairwise combinations of pdbs
    pairs = []
    for i in range(len(ff)):
        for j in range(i+1, len(ff)):
            pairs.append((ff[i], ff[j]))
    
    # Create a dataframe to store pairwise differences as a list for each element in the matrix
    # Initialise with zeros
    pairwise_differences = pd.DataFrame()
    pairwise_diff_sq = pd.DataFrame()

    for pair in pairs:
        corr_1 = pd.read_csv('results/heatmaps/new/pairwise_differences_'+pair[0]+'_avg.csv', sep='\t')
        corr_2 = pd.read_csv('results/heatmaps/new/pairwise_differences_'+pair[1]+'_avg.csv', sep='\t')

        # Ignore the first column and first row
        corr_1 = corr_1.iloc[:, 1:]
        corr_2 = corr_2.iloc[:, 1:]

        # if pair_diff is not empty, append the new pair_diff to the existing pair_diff
        if pairwise_differences.empty:
            pairwise_differences = np.abs(corr_1 - corr_2)

            pairwise_diff_sq = (corr_1 - corr_2)**2
        else:
            pairwise_differences += np.abs(corr_1 - corr_2)
            pairwise_diff_sq += (corr_1 - corr_2)**2

    print(pairwise_differences)

    pairwise_differences = pairwise_differences / len(pairs)
    pairwise_diff_sq = pairwise_diff_sq / len(pairs)
    pairwise_diff_std = np.sqrt(pairwise_diff_sq - pairwise_differences**2)

    pairwise_differences.to_csv('results/heatmaps/new/pairwise_differences_avg.csv', sep='\t')
    pairwise_diff_std.to_csv('results/heatmaps/new/pairwise_differences_std.csv', sep='\t')

    plt.figure(figsize=(12,12))
    plot_heatmap(pairwise_differences.T)
    
    plt.savefig('results/heatmaps/new/pairwise_differences_avg.png')
    plt.close()
    
    plt.figure(figsize=(12,12))
    plot_heatmap(pairwise_diff_std.T)
    plt.savefig('results/heatmaps/new/pairwise_differences_std.png')
    plt.close()

   
        
# main()
ffs = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
for ff in ffs:
    get_pairwise_differences(ff)
get_pairwise_differences_overall()
