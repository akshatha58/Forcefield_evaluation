"""
get_ss_groups.py 
Akshatha 
28 August 2024

Contains helper functions that obtain secondary structure information from the offset file
Useful to represent data in terms of secondary structures in plots

"""
import pandas as pd
import os

def get_residues(data_res):
    """
    Function to extract the helix and sheet residue numbers from the pdb_offsets.txt file
    """
    helices = []
    sheets = []
    dummy = 0
    # Iterate over the row
    for index, row in data_res.iterrows():
        # split the columns of each row
        for col in row:
            dummy += 1
            # get the residue number
            res = col.split(" ")
            for r in res:
                if r.isdigit() and dummy == 1:
                    helices.append(int(r))
                if r.isdigit() and dummy == 2:
                    sheets.append(int(r))
            
    # print(helices)
    # print(sheets)
    return helices, sheets

def get_ranges(list):
    """
    Get starting and stopping residues for each secondary structure
    Useful for colour coding plots according to secondary structure tags
    """
    ranges = []
    for i in range(len(list)-1):
        if i == 0:
            ranges.append(list[i])
        if list[i] + 1 != list[i+1]:
            ranges.append(list[i])
            ranges.append(list[i+1])
    ranges.append(list[-1])
        
    return ranges

def get_ss_data(pdbpath, flag):
    """
    Uses the helper functions get_residues and get_ranges to extract the ranges for secondary structure elements
    Returns two lists: helix_ranges and sheet_ranges having the form [start1, stop1, start2, stop2...]
    Useful for colour coding plots with secondary structure tags
    """
    pdb = pdbpath.split("/")[-1].split(".")[0]
    print(pdb)

    dir = os.path.dirname(os.path.abspath(pdbpath))
    print(dir)

    if flag == 1:
        # Offsets are needed 
        file_res = dir+"/"+pdb+"/"+pdb+"_offsets.txt"
    else:
        # Offsets not needed (start from 1)
        file_res = dir+"/ss_lists/"+pdb+"_ss.txt"
        
    data_res = pd.read_csv(file_res)

    helix_ranges = []
    sheet_ranges = []

    helices, sheets = get_residues(data_res)
    if len(helices) != 0:
        helix_ranges = get_ranges(helices)
    if len(sheets) != 0:
        sheet_ranges = get_ranges(sheets)

    print(helix_ranges)
    print(sheet_ranges)

    return helix_ranges, sheet_ranges
