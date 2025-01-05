"""
process_traj.py
Akshatha
6 October 2024

Get md_untag.gro and md_untag.xtc files and rename their residue and atom numbers to start from 1
Rename HID/HIP/HIE residues to HIS
Save them as md_mod.gro and md_trajmod.xtc files 

python3 process_traj.py <pdbpath> <trajfile> <topfile>
"""
import MDAnalysis as mda
from matplotlib import pyplot as plt
import numpy as np
import sys
import warnings
import os

# Inputfilepaths

pdbpath=sys.argv[1]
trajfile=sys.argv[2]
topfile=sys.argv[3]
pdb = os.path.splitext(os.path.basename(pdbpath))[0]

# suppress MDAnalysis warnings (related to PSF files)
warnings.filterwarnings('ignore')
print("Using MDAnalysis version", mda.__version__)

# Define file paths for input and output directory = os.path.dirname(path)
path = os.path.dirname(pdbpath)
results_path=path+"/Analysis/processed_trajs"
mdmod_gro=results_path+"/"+pdb+"_md_mod.gro"
mdtrajmod_xtc=results_path+"/"+pdb+"_md_trajmod.xtc"


# Load the GRO file
u = mda.Universe(topfile)

# Open a writer to save the modified .gro file (single frame)
with mda.Writer(mdmod_gro, n_atoms=u.atoms.n_atoms) as W:
    # Iterate over each frame in the trajectory
    for ts in u.trajectory:
        # Adjust residue numbers to start from 1
        for i, residue in enumerate(u.residues):
            residue.resid = i + 1
            if residue.resname in ["HIP", "HID", "HIE", "HSD", "HSP"]:
                residue.resname = "HIS"
            if residue.resname == "GLH":
                residue.resname = "GLU"
            if residue.resname == "ASH":
                residue.resname = "ASP"
        # Write the modified frame to the new GRO file
        W.write(u.atoms)

print("Modified the topology file.")

# Load the topology and trajectory files
try:
    u = mda.Universe(topfile, trajfile)
except Exception as e:
    print(f"Error loading files: {e}")
    with open(path+'chem_shift_data/input.gro', 'r') as f:
        print(f.read())
    raise

print("Successfully loaded the trajectory.")

# Open a writer to save the modified trajectory
with mda.Writer(mdtrajmod_xtc, n_atoms=u.atoms.n_atoms) as W:
    # Iterate over each frame in the trajectory
    for ts in u.trajectory:
        # print(f"Processing frame {ts.frame}")
        # Adjust residue numbers to start from 1
        for i, residue in enumerate(u.residues):
            residue.resid = i + 1
            # If the residue name is "HIP", change it to "HIS"
            if residue.resname in ["HIP", "HID", "HIE", "HSD", "HSP"]:
                residue.resname = "HIS"
            if residue.resname == "GLH":
                residue.resname = "GLU"
            if residue.resname == "ASH":
                residue.resname = "ASP"
        # Write the modified frame to the new XTC file
        W.write(u.atoms)

print("Finished processing all frames")
