# main.sh
# Akshatha
# 6 October 2024

# Main code for executing all analyses

path=$1
offset=$2
pdb=$(basename $path)
mainpath=$(pwd)

# j_const_data=('1ozi', '1z9b', '2gmo', '2juo')

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
echo "Computing chemical shifts..."
cd $path/Analysis/NMR_analysis
mkdir chem_shifts
cd $mainpath/NMR_analysis
./nmr_main.sh $pdb $offset 

echo "Computing J3 constants..."
cd $path/Analysis/NMR_analysis
mkdir J3_constants
cd $mainpath/NMR_analysis
./nmr_main.sh $pdb $offset 2
