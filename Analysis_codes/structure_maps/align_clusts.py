"""
align_clusts.py

Pymol code to align a representative cluster to the original NMR file

run align_clusts.py
"""
from pymol import cmd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys
import pickle

clust_struct = "single_frame" #change name
pdb = "1ozi" # change name

cmd.color("red", clust_struct)
cmd.color("blue", pdb)
cmd.align(clust_struct, pdb)
cmd.bg_color("white")
