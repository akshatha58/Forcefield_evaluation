# Extending the simulation to 1 microsecond and concatenatic trajectories

/usr/local/gromacs/bin/gmx convert-tpr -s md.tpr -extend 500000 -o md_extend.tpr

# # first time:
taskset -c 1-128 /usr/local/gromacs/bin/gmx mdrun -v -deffnm md_extend -cpi md.cpt -noappend -nt 125
# taskset -c 129-256 /usr/local/gromacs/bin/gmx mdrun -v -deffnm md_extend -cpi md.cpt -noappend -nt 125

# # later simulations: (change the checkpoint file to md_extend.cpt)
taskset -c 129-256 /usr/local/gromacs/bin/gmx mdrun -v -deffnm md_extend -cpi md_extend.cpt -noappend -nt 125
# taskset -c 1-128 /usr/local/gromacs/bin/gmx mdrun -v -deffnm md_extend -cpi md_extend.cpt -noappend -nt 125

# Do a pairwise concatenation if the following code does not work
/usr/local/gromacs/bin/gmx trjcat -f md.xtc md_extend.part*.xtc -o md_final.xtc
/usr/local/gromacs/bin/gmx eneconv -f md.edr md_extend.part*.edr -o md_final.edr
