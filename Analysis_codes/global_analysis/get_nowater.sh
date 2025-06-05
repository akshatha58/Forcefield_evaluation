#!/bin/bash
path=$1
pdb=$(basename $path)

sourcepath=$(pwd)
results_path=$path"/Analysis_equi/processed_trajs"
mdmod_gro=$results_path"/"$pdb"_md_mod.gro"
mdtrajmod_xtc=$results_path"/"$pdb"_md_trajmod.xtc"

traj_nowater=$results_path/$pdb"_nowater.xtc"
traj_nowater_gro=$results_path/$pdb"_trjnowater.gro"
struct_nowater=$results_path/$pdb"_nowater.gro"

printf "1 1" | /usr/local/gromacs/bin/gmx trjconv -f $mdtrajmod_xtc -s $mdmod_gro -o $traj_nowater -dt 200 -pbc whole -center
printf "1 1" | /usr/local/gromacs/bin/gmx trjconv -f $mdtrajmod_xtc -s $mdmod_gro -o $traj_nowater_gro -dt 200 -pbc whole -center
printf "1 1" | /usr/local/gromacs/bin/gmx trjconv -f $mdmod_gro -s $mdmod_gro -o $struct_nowater -pbc whole -center

