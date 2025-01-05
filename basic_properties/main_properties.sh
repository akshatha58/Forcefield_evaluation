# main_properties.sh
# Akshatha
# 6 October 2024

path=$1
pdb=$(basename $pdbpath)
echo $pdb

process_trajpath=$path/Analysis/processed_trajs
results_path=$path/Analysis/basic_properties

# Define paths for processed traj and top files
# trajfile=$process_trajpath/$pdb"_md_trajmod.xtc"
# topfile=$process_trajpath/$pdb"_firstframe.pdb"

cd $path

printf "11 0" | /usr/local/gromacs/bin/gmx energy -f md.edr -o $results_path/potential.xvg >> $results_path/metrics.txt
printf "13 0" | /usr/local/gromacs/bin/gmx energy -f md.edr -o $results_path/total_energy.xvg >> $results_path/metrics.txt
printf "15 0" | /usr/local/gromacs/bin/gmx energy -f md.edr -o $results_path/temperature.xvg >> $results_path/metrics.txt
printf "17 0" | /usr/local/gromacs/bin/gmx energy -f md.edr -o $results_path/pressure.xvg >> $results_path/metrics.txt
printf "23 0" | /usr/local/gromacs/bin/gmx energy -f md.edr -o $results_path/density.xvg >> $results_path/metrics.txt

# Ramachandran plot
/usr/local/gromacs/bin/gmx rama -f md.xtc -s md.tpr -o $results_path/rama.xvg -b 500000 -dt 5000

# RMSD
printf "4 4" | /usr/local/gromacs/bin/gmx rms -f md.xtc -s md.tpr -o $results_path/rmsd.xvg -tu ns

# Radius of gyration
printf "1" | /usr/local/gromacs/bin/gmx gyrate -f md.xtc -s md.tpr -o $results_path/rg.xvg -tu ns

# RMSF
printf "4" | /usr/local/gromacs/bin/gmx rmsf -f md.xtc -s md.tpr -res -o $results_path/rmsf_backbone.xvg

# python3 /home/dynamics/akshatha/Final_simfiles/Analysis_optimised/plotting_codes/rmsf_modularise.py $pdb

./rmsd_ss.sh $path

