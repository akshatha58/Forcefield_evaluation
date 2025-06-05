# get_clusters.sh
# 27 February 2025

# get_clusters.sh <directory path>

dir=$1
pdb=$(basename $dir)
resultsdir=$dir/Analysis_equi
mkdir clusters
GROMACS=/usr/local/gromacs/bin/gmx

trajfile=$resultsdir/processed_trajs/$pdb"_md_trajmod.xtc"
structfile=$resultsdir/processed_trajs/$pdb"_md_mod.gro"
indexfile=$dir/index.ndx

# cluster a subset of regions (flexible, rigid) together, specify groups in the index file according to the sequence generated in ndx_groups.py

$GROMACS make_ndx -f $structfile -o $indexfile
$GROMACS cluster -f $trajfile -s $structfile -method gromos -cutoff 0.25 -n $indexfile
# mv rmsd-clust.xpm $resultsdir/$pdb"_rmsd-clust.xpm"
# mv rmsd-dist.xvg clusters/rigid/$pdb"_4.xvg"
# rm rmsd-clust.xpm
# rm cluster.log
