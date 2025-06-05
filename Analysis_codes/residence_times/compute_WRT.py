"""
water_residence_times.py
Version 3
25 May 2025

Computes water residence times for the first hydration shell of a protein from a given trajectory.
Partitions these into residue-specific residence times and calculates the average residence time for each residue.
"""

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd
import sys, os

plot_dir = os.path.abspath("../initial_files")
if plot_dir not in sys.path:
    sys.path.append(plot_dir)

import get_all_ss as get_ss

plt.rcParams['font.size'] = 17
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# pdb=sys.argv[1]  # PDB ID
# ff=sys.argv[2]   # Force field used (in words)

ff_dict = {'Final_simfiles':1,
           'Simfiles':2,
           'amber_tip3p':3,
           'charmm_tip4p':4}

def get_retention_times(pdb, ff):
    if ff == 'Final_simfiles':
        topology = "/home/dynamics/akshatha/"+ff+"/"+pdb+"/"+pdb+"_md_mod.gro"
        trajectory = "/home/dynamics/akshatha/"+ff+"/"+pdb+"/"+pdb+"_md_trajmod.xtc"
    else:
        topology = "/home/dynamics/akshatha/"+ff+"/"+pdb+"/md_resolved.tpr"
        trajectory = "/home/dynamics/akshatha/"+ff+"/"+pdb+"/md_resolved.xtc"
    uni = mda.Universe(topology, trajectory)

    protein = uni.select_atoms("protein")
    protein_ca = uni.select_atoms("name CA")

    # Extract waters within 5 angstroms of the protein surface: frame 1
    hydration_shell_prev = uni.select_atoms("name OW and around 5 protein")
    print(f"Number of water molecules in the first hydration shell: {len(hydration_shell_prev)}")
    print(f"Total number of water molecules: {len(uni.select_atoms('name OW'))}")

    # Initialise a dictionary to store retention times for each water molecule.
    retention_time = {water.index: [[0] for _ in range(len(protein_ca))] for water in uni.select_atoms("name OW")}

    # iterate over trajectory frames and waters in hydration shell
    for ts in uni.trajectory:
        # Get the current positions of the hydration shell waters
        hydration_shell_current = uni.select_atoms("name OW and around 5 protein")
        
        # Precompute protein atom positions for this frame
        protein_positions = protein.positions

        # Build a mapping from water index to its position for quick lookup
        hydration_shell_current_dict = {w.index: w.position for w in hydration_shell_current}
        hydration_shell_prev_dict = {w.index: w.position for w in hydration_shell_prev}

        # Only iterate over water indices present in both previous and current hydration shell
        common_water_indices = set(hydration_shell_prev_dict.keys()) & set(hydration_shell_current_dict.keys())

        for water_idx in common_water_indices:
            water_pos = hydration_shell_current_dict[water_idx]
            # Find closest residue in current frame
            dists = np.linalg.norm(protein_positions - water_pos, axis=1)
            closest_idx = np.argmin(dists)
            closest_resid = protein.residues[protein.atoms[closest_idx].resindex].resid

            # Find closest residue in previous frame
            prev_water_pos = hydration_shell_prev_dict[water_idx]
            prev_dists = np.linalg.norm(protein_positions - prev_water_pos, axis=1)
            prev_closest_idx = np.argmin(prev_dists)
            closest_resid_prev = protein.residues[protein.atoms[prev_closest_idx].resindex].resid

            
            target_entry = retention_time[water_idx][closest_resid-1]
            if closest_resid == closest_resid_prev:
                # Add counter to last entry if the water is still next to the same residue
                target_entry[-1] += 1
            else:
                # If the water moved to a different residue in the previous frame and returned back to this residue, append a new entry
                target_entry.append(1)
                

        # Account for waters that were not in the hydration shell in the previous frame
        additional_waters = set(hydration_shell_current.indices) - set(hydration_shell_prev.indices)
        if additional_waters:
            water_indices = list(additional_waters)
            water_positions = hydration_shell_current.positions[
                [hydration_shell_current.indices.tolist().index(idx) for idx in water_indices]
            ]
            # Vectorized distance calculation
            dists = np.linalg.norm(protein.positions[None, :, :] - water_positions[:, None, :], axis=2)
            closest_indices = np.argmin(dists, axis=1)
            for i, water_idx in enumerate(water_indices):
                closest_resid = protein.residues[protein.atoms[closest_indices[i]].resindex].resid
                target_entry = retention_time[water_idx][closest_resid-1]
                target_entry.append(1)

        if ts.frame % 10 == 0:
            print(f"Frame {ts.frame}: Number of water molecules in the first hydration shell: {len(hydration_shell_current)}")
            print(f"Additional waters in hydration shell: {len(additional_waters)}")

        hydration_shell_prev = hydration_shell_current


    # Only get water indices that have non-empty lists in retention_time
    retention_time = {k: v for k, v in retention_time.items() if any(len(lst) > 0 for lst in v)}

    # So now, retention_df has water indices: and each column corresponds to a residue, with lists of retention times.
    # water index : [[], [], []]
    # Calculate average residence time for each residue

    # get the list corresponding to each water index
    for residue_lists in retention_time.values():
        # get the list within the list
        for i, lst in enumerate(residue_lists):
            if len(lst) > 0:
                residue_lists[i] = np.mean(lst)
            else:
                residue_lists[i] = 0  # or some other value to indicate no residence time

    # water index :[avg res time1, avg res time2, ...]
    retention_df = pd.DataFrame.from_dict(retention_time, orient='index', columns=[f"Residue_{i+1}" for i in range(len(protein_ca))])
    # Save the DataFrame to a CSV file
    retention_df.to_csv("retention_data/retention_times_"+pdb+"_"+str(ff_dict[ff])+".csv", index_label="Water_Index")

def label_resids(pdb):
    """
    Returns a list of residue names for a given PDB file.
    """
    file=os.path.join(os.getcwd(), '..','..','FASTA_Seqs', pdb+'.txt')
    with open(file, 'r') as file:
        lines = file.readlines()
        row = lines[0].strip()
        return list(row)
    
def return_colours(labels):
    """
    Define a colour scheme for residues depending on their secondary structure
    """
    label_dict = {1:"red", 2:"maroon", 3:"brown", 5: "blue", 4: "grey", 6:"grey", 7:"purple", 8:"pink",  9:"green", 10:"cyan", 11:"violet", 12:"black"}


    return label_dict[labels]

def plot_retention_times(pdb, f):
    """
    Plots the distribution of residence times for each residue.
    """
    retention_df = pd.read_csv("retention_data/retention_times_"+pdb+"_"+f+".csv", index_col="Water_Index")

    ss_filepath = os.path.join(os.getcwd(), '..', '..', 'mmcif_files', pdb + '_sslists.txt')
    ss_dict = get_ss.get_ss_residues(ss_filepath)
    # plt.figure(figsize=(15, 6))
    plt.subplot(5, 1, int(f)+1)
    average_value = []
    for i, column in enumerate(retention_df.columns):
        # get a subset of the column for which values are not zero
        values = retention_df[column][retention_df[column] > 0]
        if i == 0:
            print(values)
        average_value.append(values.mean())

    for residue in range(len(average_value)):
        ss_info = 12
        for key in ss_dict.keys():
            if residue+1 in ss_dict[key]:
                ss_info = key
                break
        plt.bar(residue, average_value[residue], color=return_colours(ss_info), edgecolor='black', width=0.8)

    # plt.title("Average Residence Times for Each Residue")
    # plt.xlabel("Residue Index")
    plt.grid('both')
    # plt.ylabel("Average Residence Time (frames)")

    # # # get residue names
    residues=label_resids(pdb)
    
    with open("retention_data/"+pdb+"_"+f+".txt", 'w') as file:
        file.write("Residue\tAverage Residence Time (frames)\n")
        for i, avg_time in enumerate(average_value):
            if avg_time < 50.0:
                file.write(f"{residues[i]}\t{avg_time}\n")

    plt.xticks(ticks=range(len(residues)), labels=[name for name in residues], fontsize=12)    
    # plt.tight_layout()
    # plt.savefig("residue_retention_"+pdb+"_"+ff+".png")

def plot_res_sasa(pdb):
    """
    Plots the average residue-wise SASA against the reference values in the NMR structure
    """
    
    plt.figure(figsize=(15,6))
    ffs=['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    label={'Final_simfiles':'Amber/4P-Ew', 'Simfiles':'CHARMM/3Pm', 'amber_tip3p':'Amber/3P', 'charmm_tip4p':'CHARMM/4P-Ew'}
    for ff in ffs:
        
        basepath='/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/basic_properties/'
        file=basepath+'res_sasa.xvg'
        data = pd.read_csv(file, sep='\s+', skiprows=25, names=['res', 'avg', 'std'])
        offset = data['res'].iloc[0]
        data['res'] = data['res'] - offset + 1

        # plt.subplot(2,1,1)

        ss_ranges = get_ss.get_ss_data(pdb)
        plt.fill_between(data['res']-1, data['avg'], color='grey', alpha=0.1)
        for key in ss_ranges:
            for i in range(0, len(ss_ranges[key]), 2):
                if len(ss_ranges[key]) >= 2:
                    plt.fill_between(data['res']-1, data['avg'], where=(data['res'] >= ss_ranges[key][i]) & (data['res'] <= ss_ranges[key][i+1]), color=return_colours(key), alpha=0.1)
        # plt.legend()

        plt.plot(data['res']-1, data['avg'], label=label[ff], marker='o', markersize=3, linewidth=1.5)
        plt.grid('both')
        plt.ylabel("SASA Average")
        residues=label_resids(pdb)
        plt.xticks(ticks=range(len(residues)), labels=[name for name in residues], fontsize=9,rotation=90)
        plt.yticks(rotation=90)
        plt.xlabel("Residue")
    plt.tight_layout()
    plt.savefig("sasa_residue_"+pdb+".png")

def get_colour_scale(data):
    """
    Returns colour map for each of the 20 amino acids based on sasa values
    """
    # file = pd.read_csv(f'sasa_avg_{ff}.txt', sep='\t', names=['Residue','SASA_avg',	'SASA_std'], skiprows=1)
    # Get heat map colours for each residue
    sasa = data['avg']
    residues = data['resnames']

    # convert sasa values to float
    sasa = sasa.astype(float)
    sasa_norm = (sasa - sasa.min()) / (sasa.max() - sasa.min())
    # Create a colour map based on the normalized sasa values
    # get a high contrast colour map
    cmap = plt.cm.get_cmap('plasma')
    colours = [cmap(sasa_norm[i]) for i in range(len(sasa_norm))]
    # Create a dictionary mapping residues to their corresponding colours
    colour_list = [colours[i] for i in range(len(residues))]

    # get the colour bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=sasa.min(), vmax=sasa.max()))
    sm.set_array([])  # Only needed for older versions of matplotlib to avoid warnings

    return colour_list, sm

# print(get_colour_scale('Final_simfiles'))  # Call to initialize the colour scale

def plot_all_ffs():
    """
    Plots results of all force fields for a given PDB
    """
    pdb = ['1ozi', '2gmo', '2juo', '2k0q', '2l4x', '2m68', '4bwh']
    ff = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    fig=plt.figure(figsize=(15, 10))
    for i in range(len(ff)):
        plt.subplot(4, 1, i+1)
        for protein in pdb: 
            if os.path.exists(f"retention_data/"+protein+"_"+str(i+1)+".txt"):
                data = pd.read_csv(f"retention_data/"+protein+"_"+str(i+1)+".txt", skiprows=1, names=['res', 'resid_time'], sep='\t')
                # sort the data by residue
                # data = data.sort_values(by='res')
                average_value = data['resid_time']
                residues = data['res']

                colour_dict, sm = get_colour_scale(ff[i])
                # add colour to dataframe
                colours = [colour_dict[res] for res in residues]
                plt.bar(residues, average_value, edgecolor='black', width=0.8, alpha=0.5, color=colours)
                
                plt.title(ff[i])
                plt.ylim([0,50])
                plt.grid('both')
                # plot colour bar on the right side
        # ax = plt.gca()
        # cbar = fig.colorbar(sm, orientation='vertical', ax=ax, fraction=0.025)


    # plt.tight_layout()
    plt.xlabel("Residue", fontsize=25)
    # put colour bar on right side
    # Add colorbar to the right side of the figure
    cbar = fig.colorbar(sm, ax=fig.get_axes(), orientation='vertical', fraction=0.025, pad=0.02)
    cbar.set_label("Average Accessible Surface Area", fontsize=18)
    fig.text(0.07, 0.5, "Average Retention Time (ps)", ha='center', va='center', fontsize=25, rotation=90)
    # add space between subplots
    plt.subplots_adjust(hspace=0.5)
    plt.savefig(f"retention_data/residue_retention.png")

def plot(pdb):
    fig=plt.figure(figsize=(17, 10))
    colours = ['blue', 'orange', 'green', 'red']
    plt.subplot(5, 1, 1)
    for i,force in enumerate(ff_dict.keys()):
        sasafile='/home/dynamics/akshatha/'+force+'/'+pdb+'/Analysis_equi/basic_properties/res_sasa.xvg'
        data = pd.read_csv(sasafile, sep='\s+', skiprows=25, names=['res', 'avg', 'std'])
        offset = data['res'].iloc[0]
        data['res'] = data['res'] - offset + 1
        plt.bar(data['res']-1, data['avg'], edgecolor='black', width=0.8, alpha=0.5, color=colours[i])
    plt.ylabel("SASA\nAverage", fontsize=18)
    plt.xticks(fontsize=1)
    plt.grid('both')

    for f in range(2, 5):
        plot_retention_times(pdb, str(f))
        plt.ylim([0,100])
    plt.xlabel("Residue", fontsize=20)

    fig.text(0.08, 0.4, "Average Retention Time (ps)", ha='center', va='center', fontsize=18, rotation=90)
    plt.savefig(f"retention_data/residue_retention_{pdb}.png")
    # plot_all_ffs(pdb)

def plot_final():
    """
    Plots results of all force fields for a given PDB
    """
    pdb = ['1ozi', '2gmo', '2juo', '2k0q', '2l4x', '2m68', '4bwh']
    ff = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    # pdb = ['1ozi']
    # ff = ['Final_simfiles']

    fig=plt.figure(figsize=(17, 10))
     
    for i in range(len(ff)):
        plt.subplot(4, 1, i+1)
        for protein in pdb: 
            if os.path.exists(f"retention_data/"+protein+"_"+ff[i]+".txt"):
                sasa_data = pd.read_csv(f'/home/dynamics/akshatha/'+ff[i]+'/'+protein+'/Analysis_equi/basic_properties/res_sasa.xvg', sep='\s+', skiprows=25, names=['res', 'avg', 'std'])
                offset = sasa_data['res'].iloc[0]
                sasa_data['res'] = sasa_data['res'] - offset + 1
                residues = label_resids(protein)
                sasa_data['resnames'] = residues
                retention_data = pd.read_csv(f"retention_data/"+protein+"_"+ff[i]+".txt", skiprows=1, names=['resnames', 'resid_time'], sep='\t')

                # concatenate and remove duplicate columns
                combined_data = pd.concat([sasa_data, retention_data], axis=1)
                combined_data = combined_data.loc[:, ~combined_data.columns.duplicated()]
                
                # sort the data by residue
                # combined_data = combined_data.sort_values(by='resid_time', ascending=True)
                print(combined_data)


                colour_list, sm = get_colour_scale(combined_data)
                
                for j in range(len(combined_data['resnames'])):
                    if combined_data['resid_time'].iloc[j] < 25.0:
                        plt.bar(combined_data['resnames'].iloc[j], combined_data['resid_time'].iloc[j], edgecolor='black', width=0.8, alpha=0.5, color=colour_list[j])
                plt.title(ff[i])
                plt.ylim([0,25])
                plt.grid('both')
        ax = plt.gca()
        # increase width of colour bar
        cbar = fig.colorbar(sm, orientation='vertical', ax=ax, fraction=0.025, pad=0.05)
        if i == 1:
            cbar.set_label("Time-Averaged SASA per Residue", fontsize=18, loc='top')
        # set ticks for colour bar
        cbar.set_ticks([0, 0.5, 1, 1.5, 2])

    plt.xlabel("Residue", fontsize=25)
    # put colour bar on right side
    # Add colorbar to the right side of the figure
    # cbar = fig.colorbar(sm, ax=fig.get_axes(), orientation='vertical', fraction=0.025, pad=0.02)
    
    fig.text(0.91, 0.5, "Average Retention Time (ps)", ha='center', va='center', fontsize=20, rotation=90)
    # add space between subplots
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5)
    plt.savefig(f"retention_data/residue_retention_fin.png")

plot_final()

# get_retention_times(pdb, ff)
# plot(pdb)
# plot_all_ffs()
