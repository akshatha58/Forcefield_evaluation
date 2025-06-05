#!/bin/bash
# compute_property.sh
# Computes a property according to specified property name. 

# Usage: ./compute_property.sh <property_name> <input_file> <output_file>

prop=$1
path=$2

# Get pdb, ff from input_file
working_dir=$(pwd)
pdb=$(basename $path)
ff=$(basename $(dirname $path))

echo $pdb
echo $ff

# Define offset
if [ "$pdb" == "2gmo" ]; then
    offset=8
elif [ "$pdb" == "2l6b" ]; then
    offset=3
else
    offset=0
fi

# Define force field
if [ "$ff" == "Final_simfiles" ] || [ "$ff" == "amber_tip3p" ]; then
    forcefield="amber"
else    
    forcefield="charmm"
fi

# Compute property
if [ "$prop" == "Shifts" ]; then
    cd $working_dir/../NMR_analysis/chem_shifts
    ./run_shifts_main.sh $path $offset $forcefield

elif [ "$prop" == "RMSF" ]; then
    cd $path
    process_trajpath=$path/Analysis_equi/processed_trajs
    trajfile=$process_trajpath/$pdb"_md_trajmod.xtc"
    results_path=$path/Analysis_equi/basic_properties

    printf "4" | /usr/local/gromacs/bin/gmx rmsf -f $trajfile -s md.tpr -res -o $results_path/rmsf_backbone.xvg
    sed -i '1,17d' "$results_path/rmsf_backbone.xvg"

elif [ "$prop" == "SASA" ]; then
    cd $path
    process_trajpath=$path/Analysis_equi/processed_trajs
    trajfile=$process_trajpath/$pdb"_md_trajmod.xtc"
    topfile=$process_trajpath/$pdb"_md_mod.gro"
    results_path=$path/Analysis_equi/basic_properties
    printf "1" | /usr/local/gromacs/bin/gmx sasa -f $trajfile -s $topfile -o $results_path/sasa.xvg -surface -tu ns -or $results_path/res_sasa.xvg

elif [ "$prop" == "NH_Corr" ]; then

    cd $working_dir/../NMR_analysis/relaxation_data
    python3 get_correlations.py "$path"
    python3 analyse_data.py "$pdb" "0" "$ff"

elif [ "$prop" == "J3" ]; then
    cd $working_dir/../NMR_analysis/j3_constants
    # Define offset
    if [ "$pdb" == "2gmo" ]; then
        offset=8
    else
        offset=0
    fi
    ./run_j3_main.sh $path $offset
    
elif [ "$prop" == "NOE" ] || [ "$prop" == "T1" ] || [ "$prop" == "T2" ]; then 
    mkdir $path/Analysis_equi/NMR_analysis/relaxation_data
    cd $working_dir/../NMR_analysis
    python3 parse_star.py "$path" "$prop"

fi




