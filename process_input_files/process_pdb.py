"""
process_pdb.py
Akshatha
17 June 2024

Processes the PDB file, adds the last row of element information for DSSP to parse the file 
Checks for the starting residue number and if there are any discontinuous residues in the PDB file, and updates the information to PDBID_offsets.txt

Usage: python3 process_pdb.py 1r36
where 1r36.pdb is the PDB file 

"""


from Bio.PDB import PDBParser, PDBIO, Select
import pandas as pd
import sys
import os

import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', DeprecationWarning)

# Get pdb file from command line input
pdb_filepath = sys.argv[1]
pdb_file = os.path.splitext(os.path.basename(pdb_filepath))[0]

# Add the element to each line of the ATOM entry in the pdb file
class ElementAddingSelect(Select):
    def accept_atom(self, atom):
        atom.element = atom.get_id()[0]
        return True

parser = PDBParser()
structure = parser.get_structure('structure', pdb_filepath)

# Check for discontinuous residues or residues not starting with 1
# return a text file with the first line having the residue offset (if the residue number doesn't start with 1)
# and the subsequent lines having the residue number of the discontinuous residues (if present)

with open(pdb_file+"_offsets.txt", "w") as f:

    for model in structure:
        for chain in model:
            prev_residue_num = None
            first_residue_num = chain.get_unpacked_list()[0].get_id()[1]

            if first_residue_num != 1:
                print(f"First residue in chain {chain.get_id()} starts with {first_residue_num} instead of 1")
            
            f.write(str(first_residue_num - 1)+"\n")

            for residue in chain:
                hetero_flag, residue_num, insertion_code = residue.get_id()
                if prev_residue_num is not None and residue_num != prev_residue_num + 1:
                    print(f"Discontinuous residue found in chain {chain.get_id()} at residue number {residue_num}")
                    f.write(str(residue_num)+"\n")

                prev_residue_num = residue_num

            
f.close()

io = PDBIO()
io.set_structure(structure)
io.save(pdb_filepath, ElementAddingSelect())
