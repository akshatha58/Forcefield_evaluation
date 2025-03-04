"""
process_pdb.py
4 December 2024

Processes the pdb file to start residue number from 1

python3 process_pdb.py <path of pdb file>

"""
import sys
import warnings
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, PDBIO

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', DeprecationWarning)

def renumber_residues(pdb_file):
    
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)

    new_residue_id = 1
    for model in structure:
        for chain in model:
            for residue in chain:
                residue.id = (residue.id[0], new_residue_id, residue.id[2])
                new_residue_id += 1

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

if __name__ == "__main__":
    pdb_file = sys.argv[1]
    renumber_residues(pdb_file)
