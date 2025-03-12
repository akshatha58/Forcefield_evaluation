"""
heatmap.py
7 March 2025

Constructs a heat map for statistics across proteins, force fields, and secondary structure elements
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.patches as mpatches
import sys
import os
import global_plots as gplots

datafile = sys.argv[1]
stat = os.path.splitext(os.path.basename(datafile))[0]

plt.rcParams['font.size'] = 13
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

def get_y_labels(datafile):
    y_label = []
    ffdict = {"Final_simfiles":1, "Simfiles":2, "amber_tip3p":3, "charmm_tip4p":4}

    with open(datafile, 'r') as file:
        for line in file:
            if line.startswith('pdb') or line.startswith('ff'):
                continue
            part1 = line.split()[0]  # Consider first column
            part2 = line.split()[1]  # Consider second column
            y_label.append(str(part1)+"-"+str(ffdict[part2]))
    return y_label


# Load data from text file
data = []
with open(datafile, 'r') as file:
    for line in file:
        if line.startswith('pdb') or line.startswith('ff'):
            continue
        parts = line.split()[2:]  # Skip the first two columns
        data.append([float(x) for x in parts])

data = np.array(data)

# Replace -800.000 with np.nan
data[data == -800.000] = np.nan

# Custom labels for x and y axes
x_labels = ['0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11','12']
y_labels = get_y_labels(datafile)

fig, ax = plt.figure(figsize=(8,11)), plt.gca()
# Create the heatmap
cax = ax.imshow(np.abs(data), cmap='coolwarm', interpolation='nearest')

# Add a color bar
fig.colorbar(cax)

# Add title and labels
plt.title('Heatmap of '+stat)
plt.xlabel('Secondary structures')
plt.ylabel('Protein-Force field/water model')

# Set custom x and y labels
ax.set_xticks(ticks=np.arange(len(x_labels)))
ax.set_yticks(ticks=np.arange(len(y_labels)))
ax.set_xticklabels(x_labels)
ax.set_yticklabels(y_labels)

ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
ax.grid(which='minor', color='black', linestyle='-', linewidth=1)
ax.tick_params(which='minor', size=0)

# Annotate each cell with its value
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        value = data[i, j]
        if not np.isnan(value):
            plt.text(j, i, f'{value:.2f}', ha='center', va='center', color='black', fontsize=9)


plt.savefig('comparisons/heatmaps/heatmap_plots/heatmap_'+stat+".png")
