# chem_shifts.sh
# Akshatha
# 7 October 2024

# Calculates chemical shifts from the trajectory and compares with experimental data
# ./chem_shifts.sh <path> <ff>

path=$1
ff=$2 # amber or charmm

# Get paths for top and trajfiles
pdb=$(basename $path)
traj_path=$path"/Analysis/processed_trajs"

mdmod_gro=$traj_path"/"$pdb"_md_mod.gro"
mdtrajmod_pdb=$traj_path"/"$pdb"_md_trajmod.pdb"

database="databases/db_2_2017_py3.db.gz"
shifts=$path"/Analysis/NMR_analysis/chem_shifts/"$pdb"_mdshifts.txt"

if [ "$ff" == 'amber' ] 
then
    cd $path/../../MD2NMR_amber
elif [ "$ff" == 'charmm' ]
then
     cd  $path/../../MD2NMR_charmm
fi

source /home/dynamics/akshatha/myenv/bin/activate # CHANGE PATH FOR ENVIRONMENT

./shifts.py --topo $mdmod_gro --traj $mdtrajmod_pdb -d $database >> $shifts

