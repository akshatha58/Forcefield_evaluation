# Modular-Fold-based-FF
Codes for implementing workflows involved in developing fold-based force field parameters for globular proteins

Folders for various workflows:

1. process_input_files
   a. Has codes to tag secondary structure elements in the topology files created by pdb2gmx
      (Needs installation of DSSP and BioPython modules)
```
./TRA.sh <pdbfilepath> <topfilepath>
```
   b. Codes for modifying the AMBER ff99SB-ILDN force field with tagged atom types and residues (also tags waters and metals, which need to be removed)
```
./modify_ffnonbonded.sh <atomtypespath> <ffnonbondedpath>
```
   c. Adds tags for bonds, angles, dihedrals (copy block of lines to an input file, copy back the modified block to the corresponding location in the ffbonded.itp file)
   
```
./modify_ffbonded.sh <inputfilepath> <outputfilepath> <paramtype>
# <paramtype> can be "BONDS" "ANGLES" or "DIHEDRALS"
``` 

2. process_trajectories
   -- Has codes to remove ss tags and renumber residues from 1 for ease of processing
   
