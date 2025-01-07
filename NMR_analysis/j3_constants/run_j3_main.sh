# run_j3_main.sh 
# 6 January 2025

# Main code to parse, compute and plot J3 coupling constant data

path=$1
offset=$2
working_dir=$(pwd)

cd .. 

Parse experimental J3 coupling constant data
echo "Parsing experimental coupling constants..."
python3 parse_star.py $path J3

# Remove water molecules from trajectory files (format: .xtc and .gro structure file)
cd $trajpath
printf "1" | /usr/local/gromacs/bin/gmx trjconv -f $pdb"_md_trajmod.xtc" -s $pdb"_md_mod.gro" -o $pdb"_traj_nowater.xtc"
printf "1" | /usr/local/gromacs/bin/gmx trjconv -f $pdb"_md_mod.gro" -s $pdb"_md_mod.gro" -o $pdb"_md_mod_nowater.gro"

# Get J3 constant data from trajectory
echo "Getting coupling constant data from trajectory..."
cd $working_dir
python3 j3_constants.py $path $offset

# Plot deviations between NMR and MD data
echo "Plotting and acquiring statistics for J3 constants..."
cd /home/dynamics/akshatha/Final_simfiles/Analysis_optimised/plotting_codes
python3 compare_J3.py $path






