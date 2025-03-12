"""
process_clusts.py

Pymol script for secondary-structure wise colouring, marking disulphide bonds, etc.
run process_clusts.py
"""

from pymol import cmd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys
import pickle

protein = "single_frame" # change name 
load protein+".pdb"

# Colour sec struct wise
util.cbss(protein, "cyan", "magenta", "salmon", _self=cmd)

# Colour disulphide bonds: specify chain and residue information in the command
# cmd.select("disulphides", "chain A and resi 6+7+11+20 or chain B and resi 7+19")
# cmd.color("tv_yellow", "disulphides")

# Change background and save the file
cmd.bg_color("white")
cmd.png(protein+".png", ray=1, width=1000, height=1000, dpi=300)
