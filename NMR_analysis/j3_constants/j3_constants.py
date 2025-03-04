"""
j3_constants.py

Updated code for calculating J3 constants from the trajectory
Uses mdtraj (trajectories without water)

python3 j3_constants.py <path> <offsets>
"""

import mdtraj
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import get_ss_groups as ssgroups
# import plotting_codes.plots as plotfn
import sys
import os

############################################################## LOADING TRAJECTORY AND TOPOLOGY FILES ##################################################

path = sys.argv[1]
offset = sys.argv[2]

pdb = os.path.splitext(os.path.basename(path))[0]
exp_path = path+"/Analysis/NMR_analysis/J3_constants/"+pdb+"_Jexp.txt"
combined_path = path+"/Analysis/NMR_analysis/J3_constants/"+pdb+"_Jcombined.txt"
ss_filepath = "/home/dynamics/akshatha/Final_simfiles/ss_lists/"+pdb+"_ss.txt" 
md_path = path+"/Analysis/NMR_analysis/J3_constants/"+pdb+"_Jmd.txt"

currentpath=os.getcwd()

def load_traj(path):
    """
    Loads the trajectory files into mdtraj
    """
    print("Loading trajectory...")

    pdb = os.path.splitext(os.path.basename(path))[0]
    traj_path = path +"/Analysis/processed_trajs/"+pdb+"_traj_nowater.xtc"
    topol= path+"/Analysis/processed_trajs/"+pdb+"_md_mod_nowater.gro"
    trajectory = mdtraj.load(traj_path, top=topol)

    print(trajectory, "\n")

    print("Accessing topology...")

    md_topology = trajectory.topology
    print(md_topology, "\n")
    table, bonds = md_topology.to_dataframe()

    return trajectory, md_topology

def calculate_J3(trajectory, md_topology, offset):
    """
    Calculates J3 coupling constants for the trajectory according to the coefficients derived by Bax et al. (2007)
    Returns time-averaged J3 values (mean and standard deviations)

    """
    print("Calculating HN-HA J-coupling constants from trajectory...\n")
    [indices, J] = mdtraj.compute_J3_HN_HA(trajectory, model='Bax2007')

    # Get time-averaged J coupling values for each dihedral, along with std deviation
    print("Calculating time averaged data...\n")
    mean_J = np.mean(J, axis=0)                              # Gly val val val: residue 2 (if Gly is resid 1)
    std_dev_J = np.std(J, axis=0)
    dihedrals = np.arange(1, len(indices)+1) # The zeroth residue doesn't really have a dihedral associated, because of the N terminal

    counter=-1
    with open(md_path, "w") as f:
        for index in indices:
            res=[]
            for atom in index:
                atomname = md_topology.atom(atom)
                resid = atomname.residue
                resid_name = resid.name
                res.append(resid_name)
            counter+=1

            print(res, " : ", mean_J[counter], " : ", dihedrals[counter])
            f.write(f"{int(dihedrals[counter]) + int(offset)}\t{res[-1]}\t{mean_J[counter]:.4f}\t{std_dev_J[counter]:.4f}\n")
    f.close()

def load_data(file_path):
    """
    Load the data from the csv file
    """
    data = pd.read_csv(file_path, sep="\s+", header=None)
    return data

def get_ss_info(ss_filepath):
    """
    Returns an array of helices and sheets to add to the data file
    """        
    ss_data = pd.read_csv(ss_filepath)
    helices, sheets = ssgroups.get_residues(ss_data)
    return helices, sheets

def concatenate_data(pdb, exp_filepath, md_filepath, combined_path, ss_filepath, offset):
    """
    Concatenates NMR and MD data based on common residues
    Adds an additional last column signifying secondary structure information for each residue
    """
    # Get J3 data from MD trajectories first
    trajectory, md_topology = load_traj(pdb)
    calculate_J3(trajectory, md_topology, offset)

    # Load experimental and md data data
    exp_data = load_data(exp_filepath)
    md_data = load_data(md_filepath)
    exp_res = list(exp_data[0])
    md_res = list(md_data[0])

    # Get helix, sheet information
    helices, sheets = get_ss_info(ss_filepath)

    # Get common residues to md and nmr data
    common_residues = [id for id in exp_res if id in md_res]
    print(helices, sheets)
    print(common_residues)
    # print("common residues: ", common_residues)

    # Concatenate data wrt common residues only
    with open(combined_path, "w") as f:
        f.write("ID   Res   J_exp  J_sd_exp  ID2  Res2  J_md  J_sd_md  ss_info  rel_err\n")

        for residue in common_residues:
            resd = residue - int(offset) + 1
            if resd in helices:
                ss_info = "H"
            elif resd in sheets:
                ss_info = "S"
            else:
                ss_info = "F"

            # Get the experimental data for common residue
            exp = exp_data[exp_data[0] == residue]
            exp = exp.to_numpy().flatten()

            print(exp, ss_info)

            # Get the MD data for common residue
            md = md_data[md_data[0] == residue]
            md = md.to_numpy().flatten()
            # print(md)

            rel_err = round((exp[2] - md[2])/exp[2], 4)
            # print(rel_err)
            
            #  Concatenate the data and save in file
            data = np.concatenate((exp, md))
            f.write('\t'.join(map(str, data)) + "\t" + ss_info + "\t" + str(rel_err)+'\n')
    f.close()

concatenate_data(path, exp_path, md_path, combined_path, ss_filepath, offset)

