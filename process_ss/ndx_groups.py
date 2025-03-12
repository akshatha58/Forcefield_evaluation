"""
ndx_groups.py 
Akshatha 
28 August 2024

Writes residue groups according to gromacs index file formatting

For example (first line, helices and second line sheets)

1r36
r 304-310 | r 323-330 | r 337-346 | r 350-353 | r 368-377 | r 386-398
r 354-356 | r 362-365 | r 410-412

Inputs: path to directory where pdb file is located 
Outputs: group.txt

python3 ndx_groups.py 1r36
"""

import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter('ignore', DeprecationWarning)
import get_all_ss as get_ss
import sys

ff = sys.argv[1] # charmm or amber
# pdb = pdbpath.split('/')[-1].split('.')[0]
# dir= pdbpath.split('/')[0]
pdbfiles = ['1ozi', '1r36', '1z9b', '2gmo', '2juo', '2jxy', '2k0q', '2l4x', '2l6b', '2lht', '2lro', '2m68', '4bwh']

# Paths
groups="groups.txt"

for pdb in pdbfiles:
    ss_ranges = get_ss.get_ss_data(pdb)
    if ff == "amber":
        idx = 17
    elif ff == "charmm":
        idx = 18

    # Print according to the following format: r 57-61 | r 83-90
    with open(groups, "a") as f:
        index = idx
        statement1 = pdb + " " + ff+ "\n"
        f.write(statement1)

        for key in ss_ranges.keys():
            if len(ss_ranges[key]) > 1:
                statement2 = ""
                for i in range(0, len(ss_ranges[key])-1, 2):
                    if ss_ranges[key][i] != ss_ranges[key][i+1]:
                        statement2 += "r "+ str(ss_ranges[key][i])+ "-"+ str(ss_ranges[key][i+1])
                    else:
                        statement2 += "r "+ str(ss_ranges[key][i])
                    if i != len(ss_ranges[key])-2:
                        statement2 += " | "
                f.write(statement2)
                f.write("\n")
                f.write("name "+str(index) + " " + str(key) + "\n")
                index += 1
                
            elif len(ss_ranges[key]) == 1 and 0 not in ss_ranges[key]:
                statement2 = "r " + str(ss_ranges[key][0])

                f.write(statement2)
                f.write("\n")
                f.write("name "+str(index) + " " + str(key) + "\n")
                index += 1
        statement3 = ""
        for id in range(idx, index):
            statement3 += "!"+str(id) + " & "
        f.write(statement3 + " 1 \nname "+ str(index) + " flexible \nq")
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

