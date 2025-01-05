# convert_fileformats.sh
# Akshatha
# October 6 2024

# Original traj/top file --> untag the residues --> change residue numbers and HIS names --> convert to .pdb format 
# ./convert_fileformats.sh <pdbid> <directorypath>

# Inputs
pdb=$1 # pdb id without extension
path=$2 # Path to the main directory where the trajectory is saved

# Standard paths
sourcepath=$(pwd)
results_path=$path"/Analysis/processed_trajs"

# Trajectory and topology file paths

topfile=$path/md.gro
trajfile=$path/md.xtc

# Tagged trajectory file (xtc converted to gro for processing)
mdtraj_gro=$results_path/$pdb"_md_traj.gro"

# Untagged trajectory file (after processing using mod_grofile.sh)
mdtraj_untag_gro=$results_path/$pdb"_md_traj_untag.gro"
mdtraj_untag_xtc=$results_path/$pdb"_md_traj_untag.xtc"

# Untagged and modified topology files
mduntag=$results_path"/"$pdb"_md_untag.gro"
mdmod_gro=$results_path"/"$pdb"_md_mod.gro"

# Modified trajectory files
mdtrajmod_xtc=$results_path"/"$pdb"_md_trajmod.xtc"
mdtrajmod_pdb=$results_path"/"$pdb"_md_trajmod.pdb"

echo "Converting file formats... "
echo " "

# # Untagging .gro topology file (single frame)
./modify_grofile.sh $topfile $mduntag

# # Convert md.xtc to .gro file and untag the trajectory (250 ps frames captured)
cd $path 
printf "1 0" | /usr/local/gromacs/bin/gmx trjconv -f $trajfile -s $path/md.tpr -o $mdtraj_gro -pbc whole -center -dt 500
cd $sourcepath
./modify_grofile.sh $mdtraj_gro $mdtraj_untag_gro

# Modify .gro untagged to .xtc untagged to modify residue numbers in the file
cd $results_path
printf "0" | /usr/local/gromacs/bin/gmx trjconv -f $mdtraj_untag_gro -s $mduntag -o $mdtraj_untag_xtc

#  Use MDAnalysis to modify residue numbers and return md_mod.gro and md_trajmod.xtc
cd $sourcepath
source /home/dynamics/akshatha/myenv/bin/activate
python3 process_traj.py $path/$pdb.pdb $mdtraj_untag_xtc $mdtraj_untag_gro

cd $results_path
printf "0" | /usr/local/gromacs/bin/gmx trjconv -f $mdtrajmod_xtc -s $mdmod_gro -o $mdtrajmod_pdb
# printf "1" | /usr/local/gromacs/bin/gmx trjconv -f $mdtraj_untag_gro -s $mduntag -o $mdtrajmod_pdb
echo "Files converted successfully"
