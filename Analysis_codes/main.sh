# main.sh
# Akshatha
# 6 October 2024

# Main code for executing all analyses. Calls all the other run codes.
# ./main.sh <directory> <offset> <force field>

# offset: according to the experimental data, to map the right residues (usually 0)
# ff: amber or charmm: to execute the appropriate trajectory processing and analysis

path=$1
offset=$2
ff=$3

pdb=$(basename $path)
mainpath=$(pwd)

cd $path
mkdir Analysis
cd Analysis
mkdir processed_trajs
mkdir basic_properties
mkdir NMR_analysis

# # Processing trajectories
cd $mainpath/process_trajectories
echo "Processing trajectories..."
chmod +x main_process.sh
./main_process.sh $path 

# # Compute basic properties from the trajectory
cd $mainpath/basic_properties
echo "Computing basic properties..."
chmod +x main_properties.sh
./main_properties.sh $path

# # Compute NMR properties from the modified trajectory
cd $path/Analysis_equi/NMR_analysis
mkdir chem_shifts
mkdir J3_constants

cd $mainpath/NMR_analysis/chem_shifts
echo "Computing chemical shifts..."
./run_shifts_main.sh $path $offset $ff

echo "Computing J3 constants..."
cd $mainpath/NMR_analysis/j3_constants
./run_j3_main.sh $path $offset 
