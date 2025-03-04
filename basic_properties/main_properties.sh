# main_properties.sh
# Akshatha
# 6 October 2024

# ./main_properties.sh <path where trajectory is stored>

path=$1
pdb=$(basename $pdbpath)
echo $pdb

process_trajpath=$path/Analysis/processed_trajs
results_path=$path/Analysis/basic_properties
GROMACS=/usr/local/gromacs/bin/gmx
# Define paths for processed traj and top files
# trajfile=$process_trajpath/$pdb"_md_trajmod.xtc"
# topfile=$process_trajpath/$pdb"_firstframe.pdb"

cd $path

# Check indices before executing
printf "11 0" | $GROMACS energy -f md.edr -o $results_path/potential.xvg >> $results_path/metrics.txt
printf "13 0" | $GROMACS energy -f md.edr -o $results_path/total_energy.xvg >> $results_path/metrics.txt
printf "15 0" | $GROMACS energy -f md.edr -o $results_path/temperature.xvg >> $results_path/metrics.txt
printf "17 0" | $GROMACS energy -f md.edr -o $results_path/pressure.xvg >> $results_path/metrics.txt
printf "23 0" | $GROMACS energy -f md.edr -o $results_path/density.xvg >> $results_path/metrics.txt

# Ramachandran plot
$GROMACS rama -f md.xtc -s md.tpr -o $results_path/rama.xvg -b 500000 -dt 5000

# RMSD
printf "4 4" | $GROMACS rms -f md.xtc -s md.tpr -o $results_path/rmsd.xvg -tu ns
sed -i '1,18d' "$results_path/rmsd_all.xvg"
# Radius of gyration
printf "1" | $GROMACS gyrate -f md.xtc -s md.tpr -o $results_path/rg.xvg -tu ns
sed -i '1,27d' "$results_path/rg_all.xvg"

# RMSF
printf "4" | $GROMACS rmsf -f md.xtc -s md.tpr -res -o $results_path/rmsf_backbone.xvg
sed -i '1,17d' "$results_path/rmsf_backbone.xvg"

# RMSDs decomposed acc to secondary structures
./rmsd_ss.sh $path

# Clustering analysis
./get_clusters.sh $path

