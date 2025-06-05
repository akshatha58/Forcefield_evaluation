# run_analysis.sh
# 12 March 2025

# Automates analysis 
# Have set offset to 0, which might not give correct offsets for certain pdb files. 
# Crosscheck before running, or run files that need offset information separately.
# This automation is mainly useful for plotting already processed data. 

what_to_analyse=$1 # enter "mapdevs" or "heatmap" or "struct"
mainpath=$(pwd)

if [ "$what_to_analyse" == "mapdevs" ]; then
    ffs=("Simfiles" "Final_simfiles" "amber_tip3p" "charmm_tip4p")
    proteins=('1ozi' '2gmo' '2juo' '2k0q' '4bwh' '2l4x' '2m68')
    j3_proteins=('1ozi' '2gmo' '2juo')

    # Iterate through each filename
    for ff in "${ffs[@]}"; do 
        echo "Processing $pdb - $ff chemical shift error mapping ..."
        ffpath=$mainpath/../$ff
        filename=$mainpath/NMR_analysis/error_maps/chemshift_err_$ff.pkl
        cd $mainpath/NMR_analysis/chem_shifts
        python3 map_devs.py "$ffpath" "$filename" 

        # echo "Processing $pdb - $ff J3 constants error mapping ..."
        # ffpath=$mainpath/../$ff/
        # filename=$mainpath/NMR_analysis/error_maps/J3_err_$ff.pkl
        # cd $mainpath/NMR_analysis/j3_constants
        # python3 map_devs.py "$ffpath" "$filename"       
    
    done   


elif [ "$what_to_analyse" == "heatmap" ]; then 
    statnames=("relerr")
    params=("C.txt" "N.txt" "H.txt")

    cd $mainpath/plotting_codes
    python3 get_all_stats.py
    cd ..
    
    for stat in "${statnames[@]}"; do
        for param in "${params[@]}"; do
            dir=$mainpath/plotting_codes/comparisons/heatmaps/$stat"_"$param
                echo "Processing $stat - $param ..."
                cd $mainpath/plotting_codes
                python3 heatmap.py $dir "shifts"
        done

        dir=$mainpath/plotting_codes/comparisons/heatmaps/$stat"_J3.txt"
        echo "Processing $stat - J3 ..."
        cd $mainpath/plotting_codes
        python3 heatmap.py $dir "J3"

    done 

elif [ "$what_to_analyse" == "struct" ]; then
    pdb=("2m68")
    for prot in "${pdb[@]}";do
        simulation_dirs=("Final_simfiles" "Simfiles" "amber_tip3p" "charmm_tip4p")
        cd plotting_codes
        python3 rmsf_modularise.py "$prot" "${simulation_dirs[@]}"
    done


elif [ "$what_to_analyse" == "clusters" ]; then
    cd NMR_analysis/chem_shifts
    pdb=('1r36' '2gmo' '2juo' '2jxy' '2k0q' '2l4x' '2l6b' '2lht' '2lro' '2m68' '4bwh')
    # pdb=('2juo' '2jxy' '2k0q' '2l4x' '2l6b' '2lht' '2lro' '2m68' '4bwh')
    # pdb=('1ozi')
    ff=('Final_simfiles' 'Simfiles' 'charmm_tip4p' 'amber_tip3p')
    for prot in "${pdb[@]}"; do
        for f in "${ff[@]}";do
                path=/home/dynamics/akshatha/$f/$prot               
                if [ -d "$path/Analysis_equi" ]; then
                    cd /home/dynamics/akshatha/Analysis_codes/basic_properties
                    ./get_clusters.sh $path $f
                fi
        done
    done

elif [ "$what_to_analyse" == "2gmo_distance" ]; then
    ff=('Final_simfiles' 'Simfiles' 'charmm_tip4p' 'amber_tip3p')
    for f in "${ff[@]}";do
        path=/home/dynamics/akshatha/$f/2gmo/Analysis_equi/processed_trajs
        # a 251 a 1063
        /usr/local/gromacs/bin/gmx make_ndx -f /home/dynamics/akshatha/$f/2gmo/Analysis_equi/processed_trajs/2gmo_nowater.gro -o index.ndx
        /usr/local/gromacs/bin/gmx pairdist -f /home/dynamics/akshatha/$f/2gmo/Analysis_equi/processed_trajs/2gmo_nowater.xtc -s /home/dynamics/akshatha/$f/2gmo/Analysis_equi/processed_trajs/2gmo_nowater.gro -o /home/dynamics/akshatha/$f/2gmo/Analysis_equi/basic_properties/2gmo_distances.xvg -n index.ndx
    done

elif [ "$what_to_analyse" == "entropy" ]; then 
    ff=('charmm_tip4p')
    pdb=('4bwh')
    for prot in "${pdb[@]}"; do
        for f in "${ff[@]}";do
                path=/home/dynamics/akshatha/$f/$prot               
                if [ -d "$path/Analysis_equi" ]; then
                    mkdir $path/Analysis_equi/entropy
                    pdbfile=/home/dynamics/akshatha/$f/$prot/Analysis_equi/processed_trajs/$prot"_md_trajmod.pdb"
                    datfile=tors_next_def.dat
                    resultsfile=/home/dynamics/akshatha/$f/$prot/Analysis_equi/entropy/sample.out

                    cd /home/dynamics/akshatha/pdb2entropy-master
                    ./pdb2entropy $pdbfile $datfile $resultsfile -c 8.0 -mi -kmi 1 -nt 250  
                    
                fi
        done
    done

elif [ "$what_to_analyse" == "nhcorr" ]; then
    cd NMR_analysis/relaxation_data
    
    pdb=('1ozi' '2gmo' '2juo' '2k0q' '2l4x' '2m68' '4bwh')
    # pdb=("1ozi")
    ff=('Final_simfiles' 'Simfiles' 'charmm_tip4p' 'amber_tip3p')
    for prot in "${pdb[@]}"; do
        for f in "${ff[@]}";do
            path=/home/dynamics/akshatha/$f/$prot
            # python3 get_correlations.py "$path"
            # python3 analyse_data.py "$prot" "0" "$f"       
        done

        python3 analyse_data.py "$prot" "1" "dummy"

    done

elif [ "$what_to_analyse" == "sims" ]; then
    ff=('Simfiles' 'charmm_tip4p' 'amber_tip3p')
    # all_ffs=('Final_simfiles' 'Simfiles' 'charmm_tip4p' 'amber_tip3p')
    pdb=('2k0q' '2l4x' '2m68' '4bwh')
    # pdb=('2gmo')
    # for f in "${all_ffs[@]}";do
        # python3 WRT_V3.py "2gmo" "$f"
    # done
    for f in "${ff[@]}"; do
        for prot in "${pdb[@]}"; do
            cd /home/dynamics/akshatha/$f
            ./gromacs_commands.sh $prot
            cd /home/dynamics/akshatha/Analysis_codes/entropy
            python3 WRT_V3.py "$pdb" "$f"
        done
    done
    
fi
