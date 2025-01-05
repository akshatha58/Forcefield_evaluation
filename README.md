# Modular-Fold-based-FF
Codes for implementing workflows involved in developing fold-based force field parameters for globular proteins

Folders for various workflows:

1. FF_processing
   a. Has codes to tag secondary structure elements in the topology files created by pdb2gmx
   b. Codes for modifying the AMBER ff99SB-ILDN force field with tagged atom types and residues

2. process_trajectories
   -- Has codes to remove ss tags and renumber residues from 1 for ease of processing
   
