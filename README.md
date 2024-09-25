# Fold_and_II_Struc-based_Modular_Forcefield
Codes for implementing workflows involved in developing fold-based force field parameters for globular proteins

Folders for various workflows:

1. FF_modification
   a. Has codes to tag secondary structure elements in the topology files
   b. Codes for modifying the AMBER ff99SB-ILDN force field with tagged atom types and residues

2. FF_analysis
   a. Has codes for processing ff trajectories, such as extracting a subset of frames, untagging residues and changing file formats
   b. Codes for automating various types of calculations and trajectory analyses
   c. Plotting files (matplotlib)
   d. Custom xmgrace plot formats
