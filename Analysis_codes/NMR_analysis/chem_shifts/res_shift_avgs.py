"""
res_shift_avgs.py
Akshatha
7 October 2024 

1. Rearranges experimental shift data to get residue-wise shifts

Import data ('resnum', 'res', 'element', 'isotope', 'shift', 'error', 'ambiguity')
Output rows: ('resnum', 'res', 'N shift', 'H shift', 'C shift', 'Nstd', 'Hstd', 'Cstd')

2. Combines md and experimental shift data into a single file, along with secondary structure information

UPDATED TO INCLUDE ALL DSSP-DEFINED SECONDARY STRUCTURES IN THE COMBINED FILE.
"""

import pandas as pd
import numpy as np
import sys
import os

plot_dir = os.path.abspath("../../process_ss")
if plot_dir not in sys.path:
    sys.path.append(plot_dir)

import get_all_ss as ssgroups

path = sys.argv[1]
offset = sys.argv[2]

pdb=path.split("/")[-1]
dir=path.split("/")[-2]

exp_shift_path=path+"/Analysis/NMR_analysis/chem_shifts/"+pdb+"_expshifts.txt"
save_path=path+"/Analysis/NMR_analysis/chem_shifts/"+pdb+"_expshifts_avg.txt"
md_path = path+"/Analysis/NMR_analysis/chem_shifts/"+pdb+"_mdshifts.txt"
combined_path = path+"/Analysis/NMR_analysis/chem_shifts/"+pdb+"_expshifts_combined.txt"
ss_filepath = path+"/../../mmcif_files/"+pdb+"_sslists.txt"

# trial files
# exp_shift_path = "1ozi_expshifts.txt"
# save_path = "1ozi_expshifts_avgd.txt"
# md_path = "1ozi_mdshifts.txt"
# combined_path = "1ozi_md_exp_combd.txt"
# ss_filepath = "/home/dynamics/akshatha/Final_simfiles/ss_lists/1ozi_ss.txt"

def load_data(file_path):
    """
    Load the data from the csv file
    """
    data = pd.read_csv(file_path, sep="\s+", header=None)
    return data
    
def get_averages(data):
    """
    Get the averages of the values in the 5th column for each unique value in the 2nd and 3rd columns
    """
    # group the data by the 2nd and 3rd columns
    grouped = data.groupby([0, 1, 2])
    print(grouped)

    # get the mean of the 5th column for each group
    averages = grouped[4].mean()

    # if one of H, N, C is missing, fill it with 0
    avg_grouped = averages.unstack().fillna(-1)
    avg_grouped = avg_grouped.to_dict(orient='index')

    return avg_grouped

def get_std_devs(data):

    """
    Workflow

    Get the standard deviations of the values (from the errors in the 6th column) 
    grouped according to the value in the 1st and 3rd columns
    
    98  VAL  C    0.223607
             H    0.141421
             N    0.223607

    convert this value to a dictionary such that
    key = (98, 'VAL')
    value = [0.223607, 0.141421, 0.223607]
     
    """
    # group the data by the 1st and 3rd columns

    grouped = data.groupby([0, 1, 2])
    std_devs = grouped[5].apply(lambda x: np.sqrt(x.mean()))
    std_dev_grouped = std_devs.unstack().fillna(-1)
    # print(std_dev_grouped)
    std_dev_grouped = std_dev_grouped.to_dict(orient='index')

    # print(std_dev_grouped)
    return std_dev_grouped

def format_data(data, offset, save_path):
    """
    Loads the experimental shifts data, processes it and saves the averages and stddevs to a file according to the following format:
    
    3	GLY	101.458	5.842	43.950	9.287	1.212	3.669	1.000
    """
    shifts = load_data(data)
    avg = get_averages(shifts)
    std_dev = get_std_devs(shifts)

    with open(save_path, "w") as f:
        # #ID	Name	C	H	N	Cstd	Hstd	Nstd
        # f.write('#ID \t Name \t N \t H \t C \t Nstd \t Hstd \t Cstd \n')

        for key in avg:
            # Output instance for each entry in std_dev
            # (3, 'GLY'): {'C': 0.22360679774997896, 'H': 0.1414213562373095, 'N': 0.22360679774997896}
            avg_formatted = '\t'.join([f"{avg[key][x]:.3f}" for x in ['N','H','C']])
            std_dev_formatted = '\t'.join([f"{std_dev[key][x]:.3f}" for x in ['N','H','C']])
            f.write(f"{int(key[0]) - int(offset)}\t{key[1]}\t{avg_formatted}\t{std_dev_formatted}\n")

    f.close()

def get_ss_info(pdb):
    """
    Returns an array of helices and sheets to add to the data file
    """        
    ss_dict = ssgroups.get_ss_residues(pdb)
    print(ss_dict)

    return ss_dict

def concatenate_data(exp, md, combined_path, pdb):
    """
    Concatenates NMR and MD data based on common residues
    Adds an additional last column signifying secondary structure information for each residue
    Filters data to only include MD predicted shifts with a coverage greater than 90%

    """
    exp_data = load_data(exp)
    md_data = load_data(md)
    md_data_filtered = md_data[md_data[8] >= 0.9]
    # left_out_resids = len(md_data) - len(md_data_filtered)

    exp_res = list(exp_data[0])
    md_res = list(md_data_filtered[0])

    # Get helix, sheet information
    ss_dict = get_ss_info(ss_filepath)

    # Get common residues to md and nmr data
    common_residues = [id for id in exp_res if id in md_res]
    # print(common_residues)

    # Concatenate data pertaining to these common residues only

    with open(combined_path, "w") as f:
        f.write("ID  Res  N_exp  H_exp  C_exp  N_sd_exp  H_sd_exp  C_sd_exp  N_md  H_md  C_md  N_sd_md  H_sd_md  C_sd_md  ss_info  rel_err_N  rel_err_H  rel_err_C \n")
        
        for residue in common_residues:
            ss_info = 12
            for key in ss_dict.keys():
                if residue - int(offset) in ss_dict[key]:
                    ss_info = key
                    break
                       
            # Get the experimental data for common residue
            exp = exp_data[exp_data[0] == residue]
            exp = exp.drop([0], axis=1)
            exp = exp.to_numpy().flatten()
            # print(exp)

            # Get the MD data for common residue
            md = md_data_filtered[md_data_filtered[0] == residue]
            md = md.drop([0, 1, 8], axis=1)
            md = md.to_numpy().flatten()
            # print(md)
            if exp[1] != -1.0:
                rel_err_N = round((exp[1] - md[0])/exp[1], 4)
            if exp[1] == -1.0:
                rel_err_N = -100
            if exp[2] != -1.0:
                rel_err_H = round((exp[2] - md[1])/exp[2], 4)
            if exp[2] == -1.0:
                rel_err_H = -100
            if exp[3] == -1.0:
                rel_err_C = -100
            if exp[3] != -1.0:
                rel_err_C = round((exp[3] - md[2])/exp[3], 4)
            # print(rel_err_N, rel_err_H, rel_err_C)

            # Concatenate the data and save in file
            data = np.concatenate((exp, md))
            datastr = ""
            for d in data:
                datastr = datastr + str(d) + "\t"
            line = str(residue) + "\t" + datastr + "\t" + str(ss_info) + "\t" + str(rel_err_N) + "\t" + str(rel_err_H) + "\t" + str(rel_err_C)  + "\n"
            
            f.write(line)

    f.close()


format_data(exp_shift_path, offset, save_path)
concatenate_data(save_path, md_path, combined_path, ss_filepath)

            














