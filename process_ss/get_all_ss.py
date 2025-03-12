"""
get_ss_groups.py 
Akshatha 
6 March 2025

Contains helper functions that obtain secondary structure information from the offset file
Useful to represent data in terms of secondary structures in plots

UPDATED VERSION TO INCLUDE ALL SECONDARY STRUCTURES CLASSIFIED BY DSSP, AND SEPARATE OUT TERMINAL HELICAL RESIDUES

"""
import pandas as pd

def get_residues(data_res):
    """
    Function to extract the helix and sheet residue numbers from the pdb_offsets.txt file
    """
    # SS_names
    alpha_helix = []
    three_10_helix = []
    pi_helix = []
    beta_bridge_isolated = []
    beta_ladder = []
    hbonded_turn = []
    bends = []
    terminal_helices_alpha = []
    terminal_helices_three10 = []
    terminal_pi_helices = []
    kappa_helices = []

    ss_dict = {1:alpha_helix, 2: three_10_helix, 3: pi_helix, 4: beta_bridge_isolated, 5: beta_ladder, 6: hbonded_turn, 7:bends, 8:kappa_helices, 9:terminal_helices_alpha, 10:terminal_helices_three10, 11:terminal_pi_helices}

    dummy = 0
    # Iterate over the row
    for index, row in data_res.iterrows():
        # split the columns of each row
        for col in row:
            dummy += 1
            # get the residue number
            res = col.split(" ")
            for r in res:
                ss_dict[dummy].append(int(r))
    # print("Original ss list")
    print(ss_dict)
    return ss_dict

def get_ranges(list):
    """
    Get starting and stopping residues for each secondary structure
    Useful for colour coding plots according to secondary structure tags
    """
    ranges = []
    if len(list) > 1:
        for i in range(len(list)-1):
            if i == 0:
                ranges.append(list[i])
            if list[i] + 1 != list[i+1]:
                ranges.append(list[i])
                ranges.append(list[i+1])
        ranges.append(list[-1])
    
    if len(list) <= 1 and 0 not in list:
        ranges = list
        
    return ranges

def get_terminal_helices(ss_dict):
    """
    Extracts the terminal alpha and three10 helices for the protein
    """
    # 1: 8, 2:9, 3:10
    for key in range(1,4):
        if 0 not in ss_dict[key]:
            ss_dict[key + 8] = get_ranges(ss_dict[key])
            ss_dict[key] = [x for x in ss_dict[key] if x not in ss_dict[key + 8]]

    # print("SS list with terminals")
    # print(ss_dict)
    return ss_dict

def get_ss_data(pdb):
    """
    Main function that returns secondary structure ranges for residues
    """
    file_res = "../../mmcif_files/"+pdb+"_sslists.txt"
         
    data_res = pd.read_csv(file_res, header=None)

    ss_dict = get_residues(data_res)
    ss_dict_updated = get_terminal_helices(ss_dict)
    ss_ranges = {}
    for key in ss_dict_updated.keys():
        if 0 not in ss_dict_updated[key]:
            ss_ranges[key] = get_ranges(ss_dict_updated[key])

    print("SS Ranges")
    print(ss_ranges)

    return ss_ranges

def get_ss_residues(file_res):
    """
    Main function that returns secondary structure ranges for residues
    """
         
    data_res = pd.read_csv(file_res, header=None)

    ss_dict = get_residues(data_res)
    ss_dict_updated = get_terminal_helices(ss_dict)


    return ss_dict_updated
