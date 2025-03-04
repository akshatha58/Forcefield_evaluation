# get_clusters.sh
# 27 February 2025

# ./get_clusters.sh <directory path>

dir=$1
pdb=$(basename $dir)
resultsdir=$dir/Analysis
mkdir $resultsdir/clusters
GROMACS=/usr/local/gromacs/bin/gmx

trajfile=$resultsdir/processed_trajs/$pdb"_md_trajmod.xtc"
structfile=$resultsdir/processed_trajs/$pdb"_md_mod.gro"
indexfile=$dir/index.ndx

printf "1 1" | $GROMACS cluster -f $trajfile -s $structfile -method gromos -cutoff 0.25 -g $resultsdir/clusters/$pdb"_clusters.log" -cl $resultsdir/clusters/$pdb"_clusters.pdb" -n $indexfile
mv rmsd-clust.xpm $resultsdir/clusters/$pdb"_rmsd-clust.xpm"
mv rmsd-dist.xvg $resultsdir/clusters/$pdb"_rmsd-dist.xvg"
