# run_analysis.sh
# 12 March 2025

# Automates analysis 
# Have set offset to 0, which might not give correct offsets for certain pdb files. 
# Crosscheck before running, or run files that need offset information separately.
# This automation is mainly useful for plotting already processed data. 

what_to_analyse=$1 # enter "main" or "heatmap"
mainpath=$(pwd)

# Automates main analysis (Mainly plotting NMR observables)
if [ "$what_to_analyse" == "main" ]; then
    ffs=("Simfiles" "Final_simfiles" "amber_tip3p" "charmm_tip4p")
    proteins=('1ozi' '2gmo' '2juo' '2k0q')

    # Iterate through each filename
    for ff in "${ffs[@]}"; do
        if [ "$ff" == "Simfiles" || "$ff" == "charmm_tip4p" ]; then
            forcefield="charmm"
        elif [ "$ff" == "Final_simfiles" || "$ff" == "amber_tip3p" ]; then
            forcefield="amber"
        fi
        
        for pdb in "${proteins[@]}"; do
            dir="../$ff/$pdb"
            if [ -d "$dir" ]; then 
                echo "Processing $pdb - $ff ..."
                ./main.sh $dir 0 $forcefield
            fi
        done
    done   
    
# Automates Heat map generation for all protein/force field combinations
elif ["$what_to_analyse" == "heatmap" ]; then 
    statnames=("chisq_" "relerr_" "RMSD_")
    params=("C.txt" "N.txt" "H.txt" "J3.txt")
    
    for stat in "${statnames[@]}"; do
        for param in "${params[@]}"; do
            dir="plotting_codes/comparisons/heatmaps/$stat$param"
            if [ -d "$dir" ]; then
                echo "Processing $stat - $param ..."
                cd $mainpath/plotting_codes
                python3 heatmap.py $dir 
            fi
        done
    done 
fi
    
