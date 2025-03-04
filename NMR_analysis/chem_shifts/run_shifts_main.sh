# run_shifts_main.sh 
# 6 January 2025

# Main code to parse, compute and plot chemical shift data 
# run_shifts_main.sh <path> <offset> <ff>
# offset is 0 for most proteins, 8 for 2gmo and 3 for 2l6b
# ff is either amber or charmm

path=$1
offset=$2
ff=$3

pdb=$(basename $path)
savepath=$path/Analysis_equi/NMR_analysis
trajpath=$path/Analysis_equi/processed_trajs
dir=$(dirname $path)
working_dir=$(pwd)

# # Parse experimental chemshift data
echo "Parsing experimental shifts..."
cd ..
python3 parse_star.py $path SHIFTS

cd chem_shifts
# Get MD chemshifts from trajectory
echo "Getting MD chemshifts from trajectory..."
chmod +x chem_shifts.sh
./chem_shifts.sh $path $ff
tail -n +2 $savepath/chem_shifts/$pdb"_mdshifts.txt" > temp.txt && mv temp.txt $savepath/chem_shifts/$pdb"_mdshifts.txt"

# Compile both data based on common residues, along with ss information
echo "Integrating data..."
python3 res_shift_avgs.py $path $offset

# # Plot deviations between NMR and MD data
echo "Plotting data and acquiring statistics for shifts..."
cd $working_dir/../../plotting_codes/
python3 compare_shifts.py $path
