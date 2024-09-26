"""
ndx_groups.py 
Akshatha 
28 August 2024

Writes residue groups according to gromacs index file formatting

For example (first line, helices and second line sheets)

1r36
r 304-310 | r 323-330 | r 337-346 | r 350-353 | r 368-377 | r 386-398
r 354-356 | r 362-365 | r 410-412

Inputs: pdb id only 
Outputs: group.txt

python3 ndx_groups.py 1r36
"""

import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter('ignore', DeprecationWarning)
import get_ss_groups as get_ss
import sys

pdb = sys.argv[1]
helix_ranges, sheet_ranges = get_ss.get_ss_data(pdb, 1)

# read the first line of offsets file into offset
# with open(file_res, 'r') as f:
#     offset = f.readline().strip()

# get offsetted residue values for helix and sheet

# offset_helix = [x - int(offset) for x in helices]
# offset_sheet = [x - int(offset) for x in sheets]

# range_offset_helix = get_ranges(offset_helix)
# range_offset_sheet = get_ranges(offset_sheet)

# print(range_offset_helix)
# print(range_offset_sheet)

# Print according to the following format: r 57-61 | r 83-90
with open("/home/dynamics/akshatha/Final_simfiles/FF_analysis_codes/groups.txt", "a") as f:
    statement1 = pdb + "\n"
    f.write(statement1)

    if len(helix_ranges) != 0:
        for i in range(0, len(helix_ranges)-1, 2):
            statement2 = "r "+ str(helix_ranges[i])+ "-"+ str(helix_ranges[i+1])
            if i != len(helix_ranges)-2:
                statement2 += " | "
            
            f.write(statement2)

    f.write("\n")

    if len(sheet_ranges) != 0:
        for i in range(0, len(sheet_ranges)-1, 2):
            statement3 = "r "+ str(sheet_ranges[i])+ "-"+ str(sheet_ranges[i+1])
            if i != len(sheet_ranges)-2:
                statement3 += " | "
            
            f.write(statement3)
    f.write("\n\n")
    
    print(pdb, "groups written to groups.txt \n")

f.close()


"""
Index grouping syntaxes
r 11-15                   residue numbers 11 to 15
name 10 helices           name group 10 as helices
!10 & !11 & 1             flexible groups in the protein
del 11                    delete group
"""

