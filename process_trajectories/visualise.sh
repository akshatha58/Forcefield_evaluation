# visualise.sh 
# Akshatha
# 6 October 2024

# Get the first and last frames for visualisation
# Usage: ./visualise.sh <pdbfilepath>


pdbfilepath=$1
pdb=$(basename $pdbfilepath)
pdb="${pdb%.*}"

results_path=$pdbfilepath"/Analysis/processed_trajs"

cd $results_path

mdtrajmod_pdb=$results_path"/"$pdb"_md_trajmod.pdb"
mdmod_gro=$results_path"/"$pdb"_md_mod.gro"
first_frame=$results_path"/"$pdb"_firstframe.pdb"
last_frame=$results_path"/"$pdb"_lastframe.pdb"
frames=$results_path"/"$pdb"_frames.gro"

printf "1" | /usr/local/gromacs/bin/gmx trjconv -f $mdtrajmod_pdb -s $mdmod_gro -o $first_frame -dump 0
printf "1" | /usr/local/gromacs/bin/gmx trjconv -f $mdtrajmod_pdb -s $mdmod_gro -o $last_frame -dump 1000000

# # Get a trajectory of every 5 ns frames to visualise
printf "1" | /usr/local/gromacs/bin/gmx trjconv -f $mdtrajmod_pdb -s $mdmod_gro -o $frames -dt 1000



 
