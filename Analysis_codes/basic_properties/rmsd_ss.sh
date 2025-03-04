# rmsd_ss.sh
# 4 December 2024

# ./rmsd_ss.sh <directory where trajectory is stored>

path=$1
pdb=$(basename $path)
pdb=${pdb%.pdb}
GROMACS=/usr/local/gromacs/bin/gmx

process_trajpath=$path/Analysis/processed_trajs
results_path=$path/Analysis/basic_properties

# # Group residues according to assigned secondary structure
python3 ndx_groups.py $path

cd $path
$GROMACS make_ndx -f md.gro -o index.ndx

# Index grouping syntaxes (sample: type the following in the index file, and uncomment the last three lines of code)
# r 11-15                   residue numbers 11 to 15 (get groups from groups.txt)
# name 17 helices           name group 10 as helices
# name 18 sheets
# !17 & !18 & 1             flexible groups in the protein
# name 19 flexible
# q

# Crosscheck group numbers before executing
$GROMACS rms -f md.xtc -s md.tpr -o $results_path/rmsd_helix.xvg -tu ns -n index.ndx
$GROMACS rms -f md.xtc -s md.tpr -o $results_path/rmsd_sheet.xvg -tu ns -n index.ndx
$GROMACS rms -f md.xtc -s md.tpr -o $results_path/rmsd_flexible.xvg -tu ns -n index.ndx
