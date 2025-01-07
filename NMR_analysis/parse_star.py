"""
parse_star.py
11 December 2024

Parses NMR observables from the STAR files 

Observable names : savefile names : category names : loop names

                Entries to be appended to the parsed file

                  "SHIFTS" : shift_set_1 : assigned_chemical_shifts : Atom_chem_shift
                            Seq_ID, Comp_ID, Atom_ID, Val, Val_err, Ambiguity_code

                     "J3"  : J_coupling_constants_1 : coupling_constants : Coupling_constant
                            Seq_ID_1, Comp_ID_1, Val, Val_min, Val_max, Val_err

                    "RDC"  : residual_dipolar_couplings : RDCs : RDC
                            Seq_ID_1, Comp_ID_1, Val, Val_min, Val_max, Val_err

                   "HEPF"  : H_exchange_protection_factor_2 : H_exch_protection_factors : H_exch_protection_factor
                            Seq_ID, Comp_ID, Val, Val_err

                    "NOE"  : heteronuclear_NOE : heteronucl_NOEs : Heteronucl_NOE
                            Seq_ID_1, Comp_ID_1, Val, Val_err

                    "T1"   : T1_relaxation_label : heteronucl_T1_relaxation: T1
                            Seq_ID, Comp_ID, Atom_ID, Val, Val_err

                    "T2"   : T2_relaxation_label : heteronucl_T2_relaxation : T2
                            Seq_ID, Comp_ID, Atom_ID, T2_val, T2_val_err, Rex_val, Rex_err
"""

import pynmrstar
import numpy
import pandas
import sys
import os

# Inputs: pdb file and observable to be parsed from the STAR file
pdbpath=sys.argv[1]
observable=sys.argv[2]
pdb=os.path.basename(pdbpath).split(".")[0]
dir=os.path.dirname(pdbpath)

# File paths
expdata_path=dir+"/../NMR_BMRB/all_data/"
file = expdata_path+pdb+"_data.txt"

# Dictionary of {observable: [savefile, category, loop]}

# Input entire NMR data into full_data
full_data = pynmrstar.Entry.from_file(file)
full_data.print_tree()

if observable == "SHIFTS":
    entries = []
    filepath = dir+"/"+pdb+"/Analysis/NMR_analysis/chem_shifts/"+pdb+"_expshifts.txt"
    for chemical_shift_loop in full_data.get_loops_by_category('Atom_chem_shift'):
        entries.append(chemical_shift_loop.get_tag(['Seq_ID', 'Comp_ID', 'Atom_ID', 'Atom_type','Atom_isotope_number','Val', 'Val_err', 'Ambiguity_code']))

        with open(filepath, "w") as parsedfile:
                for line in entries:
                        for residinfo in line:
                                if residinfo[2] in ["CA", "H", "N"]:
                                        parsedfile.write(f"{int(residinfo[0])}\t{residinfo[1]}\t{residinfo[3]}\t{int(residinfo[4])}\t{float(residinfo[5]):.4f}\t{float(residinfo[6]):.2f}\t{int(residinfo[7])}\n")
        parsedfile.close()
    
elif observable == "J3":
       entries = []
       filepath =  dir+"/"+pdb+"/Analysis/NMR_analysis/J3_constants/"+pdb+"_Jexp.txt"
       for j3_loop in full_data.get_loops_by_category('Coupling_constant'):
                entries.append(j3_loop.get_tag(['Seq_ID_1', 'Comp_ID_1','Val', 'Val_err']))
                with open(filepath, "w") as parsedfile:
                        for line in entries:
                                for residinfo in line:
                                        parsedfile.write(f"{int(residinfo[0])-1}\t{residinfo[1]}\t{float(residinfo[2]):.4f}\t{float(residinfo[3]):.2f}\n")
                parsedfile.close()

elif observable == "NOE":
       entries = []
       filepath =  dir+"/"+pdb+"/Analysis/NMR_analysis/NOE/"+pdb+"_NOEexp.txt"
       for NOE_loop in full_data.get_loops_by_category('Heteronucl_NOE'):
                entries.append(NOE_loop.get_tag(['Seq_ID_1', 'Comp_ID_1','Val', 'Val_err']))
                with open(filepath, "w") as parsedfile:
                        for line in entries:
                                for residinfo in line:
                                        parsedfile.write(f"{int(residinfo[0])}\t{residinfo[1]}\t{float(residinfo[2]):.4f}\t{float(residinfo[3]):.2f}\n")
                parsedfile.close()

elif observable == "T1":
       entries = []
       filepath =  dir+"/"+pdb+"/Analysis/NMR_analysis/T1/"+pdb+"_T1exp.txt"
       for T1_loop in full_data.get_loops_by_category('T1'):
                entries.append(T1_loop.get_tag(['Seq_ID_1', 'Comp_ID_1','Atom_ID','Val', 'Val_err']))
                with open(filepath, "w") as parsedfile:
                        for line in entries:
                                for residinfo in line:
                                        parsedfile.write(f"{int(residinfo[0])}\t{residinfo[1]}\t{residinfo[2]}\t{float(residinfo[3]):.4f}\t{float(residinfo[4]):.2f}\n")
                parsedfile.close()

elif observable == "T2":
       entries = []
        # Seq_ID, Comp_ID, Atom_ID, T2_val, T2_val_err, Rex_val, Rex_err
       filepath = dir+"/"+pdb+"/Analysis/NMR_analysis/T2/"+pdb+"_T2exp.txt"
       for T2_loop in full_data.get_loops_by_category('T2'):
                entries.append(T2_loop.get_tag(['Seq_ID', 'Comp_ID','Atom_ID','T2_Val', 'T2_Val_err', 'Rex_val', 'Rex_err']))
                with open(filepath, "w") as parsedfile:
                        for line in entries:
                                for residinfo in line:
                                        parsedfile.write(f"{int(residinfo[0])}\t{residinfo[1]}\t{residinfo[2]}\t{float(residinfo[3]):.4f}\t{float(residinfo[4]):.2f}\t{float(residinfo[5]):.4f}\t{float(residinfo[6]):.2f}\n")
                parsedfile.close()

elif observable == "RDC":
       entries = []
        # Seq_ID, Comp_ID, Atom_ID, T2_val, T2_val_err, Rex_val, Rex_err
       filepath =  dir+"/"+pdb+"/Analysis/NMR_analysis/RDC/"+pdb+"_RDCexp.txt"
        # Seq_ID_1, Comp_ID_1, Val, Val_min, Val_max, Val_err

       for RDC_loop in full_data.get_loops_by_category('RDC'):
                entries.append(RDC_loop.get_tag(['Seq_ID_1', 'Comp_ID_1','Val', 'Val_min', 'Val_max', 'Val_err']))
                with open(filepath, "w") as parsedfile:
                        for line in entries:
                                for residinfo in line:
                                        parsedfile.write(f"{int(residinfo[0])}\t{residinfo[1]}\t{float(residinfo[2]):.4f}\t{float(residinfo[3]):.4f}\t{float(residinfo[4]):.4f}\t{float(residinfo[5]):.2f}\n")
                parsedfile.close()









