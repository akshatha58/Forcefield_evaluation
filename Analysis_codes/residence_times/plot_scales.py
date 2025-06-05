"""
plot_scales.py

Plots heatmaps of average SASA and water retention times for different amino acids

"""

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd
import sys, os

plot_dir = os.path.abspath("../../initial_files")
if plot_dir not in sys.path:
    sys.path.append(plot_dir)

import get_all_ss as get_ss

plt.rcParams['font.size'] = 17
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

def get_averages():
    ffs = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    for f in ffs:
        file = pd.read_csv(f'{f}_all.txt', sep='\t', index_col=0, names=['residue', 'retention'], skiprows=1)
        
        file = file[file['retention'] < 50.0]
        
        # average retention times for each residue
        file = file.groupby('residue').mean()
        # sort in ascending order of retention time
        file = file.sort_values(by='retention', ascending=True)
        with open(f'retention_50_{f}.txt', 'w') as out:
            file.to_csv(out, sep='\t', float_format='%.4f')
        
def combine_data(property):
    ff = ['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    combined_data = pd.DataFrame()
    for f in ff:
        data = pd.read_csv(f'{property}_{f}.txt', sep='\t', index_col=0)
        # delete data['SASA_std'] if it exists
        if 'SASA_std' in data.columns:
            data.drop(columns=['SASA_std'], inplace=True)
        if combined_data.empty:
            combined_data = data
        else:
            combined_data = pd.concat([combined_data, data], axis=1)
    combined_data.columns = ff

    print(combined_data)
    return combined_data

def plot_heatmap():
    """
    Plots heatmaps of average SASA and water retention times for different amino acids in a descending order
    """
    plt.subplots(1, 2, figsize=(16, 8))
    plt.subplots_adjust(wspace=0.3)

    for i, prop in enumerate(['sasa', 'retention_50']):
        data = combine_data(prop)
        if prop == 'sasa':
            title = 'Average SASA'
            ylabel = 'SASA'
        else:
            title = 'Average Water Retention Time (ps)'
            ylabel = 'Water Retention Time (ps)'

        ax = plt.subplot(1, 2, i + 1)
        # Sort the data by the first column
        # if prop == "sasa":
        data = data.sort_values(by=data.columns[0], ascending=False)
        # Reorder the index to have the same order as the first column
        data = data.reindex(data.index)
        im = ax.imshow(data, cmap='plasma', aspect='auto')
        ax.set_yticks(np.arange(len(data.index)))
        ax.set_yticklabels(data.index, rotation=45, ha='right')
        ax.set_xticks(np.arange(len(data.columns)))
        ax.set_xticklabels(np.arange(1, len(data.columns) + 1))
        ax.set_title(title)
        # ax.set_ylabel(ylabel)
        plt.colorbar(im, ax=ax)
    plt.savefig('heatmap_sorted.png', dpi=300, bbox_inches='tight')

def plot_heatmap_orders(property):
    """
    Plots heatmaps of average SASA and water retention times for different amino acids
    """
    fig, axs = plt.subplots(1,4, figsize=(5, 10))
    # plt.subplots_adjust(wspace=0.6, hspace=0.3)

    for i, ff in enumerate(['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']):
        file = pd.read_csv(f'{property}_{ff}.txt', sep='\t', index_col=0)
        if property == 'sasa':
           fig.suptitle('Average SASA')
           # SORT the data by the first column
           file = file.sort_values(by=file.columns[0], ascending=True)
           # normalise the data from 0 to 1
           file = (file - file.min()) / (file.max() - file.min())
        else:
            fig.suptitle('Average Water Retention Time (ps)')
            file = (file - file.min()) / (file.max() - file.min())

        ax = plt.subplot(1, 4, i + 1)
        # print(file)
        # get black border colours in heatmap
        im = ax.imshow(file, cmap='plasma', interpolation='nearest', aspect='equal')
        ax.set_yticks(np.arange(len(file.index)))
        ax.set_yticklabels(file.index, ha='right')
        # dont display xticks
        ax.set_xticks([])
        ax.set_xticklabels([])

        ax.set_title(i+1)

        # if property == 'sasa':
            # im.set_clim(0, 2)
    plt.colorbar(im, ax=ax, fraction=0.25)
    plt.tight_layout()
    plt.savefig(f'heatmap_{property}_norm.png')


# get_averages()
plot_heatmap_orders('sasa')
plot_heatmap_orders('retention')
