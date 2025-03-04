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
# mkdir $results_path

# Trajectory and topology file paths

topfile=$path/md.gro
trajfile=$path/md.xtc

# Tagged trajectory file
mdtraj_xtc=$results_path/$pdb"_md_traj.xtc"
mdtop_gro=$results_path/$pdb"_md_top.gro"

# Result files $mdtrajmod_xtc -s $mdmod_gro -o $mdtrajmod_pdb
mdtrajmod_xtc=$results_path/$pdb"_md_trajmod.xtc"
mdmod_gro=$results_path/$pdb"_md_mod.gro"
mdtrajmod_pdb=$results_path/$pdb"_md_trajmod.pdb"
mdtop_untag_gro=$results_path/$pdb"_mdmoduntag.gro"

echo "Converting file formats... "
echo " "

# Extract frames from trajectory file and structure file
cd $path 
# /usr/local/gromacs/bin/gmx make_ndx -n $path/index.ndx
/usr/local/gromacs/bin/gmx trjconv -f $trajfile -s $path/md.tpr -o $mdtraj_xtc -pbc whole -center -dt 500 
/usr/local/gromacs/bin/gmx trjconv -f $topfile -s $path/md.tpr -o $mdtop_gro -pbc whole -center -dt 500 

#  Use MDAnalysis to modify residue numbers and return md_mod.gro and md_trajmod.xtc
cd $sourcepath
source /home/dynamics/akshatha/myenv/bin/activate
# ./modify_grofile.sh $mdtop_gro $mdtop_untag_gro
python3 process_traj.py $path/$pdb.pdb $mdtraj_xtc $mdtop_gro

cd $results_path
/usr/local/gromacs/bin/gmx trjconv -f $mdtrajmod_xtc -s $mdmod_gro -o $mdtrajmod_pdb 
echo "Files converted successfully"

