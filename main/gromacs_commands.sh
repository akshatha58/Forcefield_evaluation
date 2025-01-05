# gromacs_commands.sh
# Akshatha
# 24 September 2024

# GROMACS commands for running simulation in CHARMM36 + TIP3Pm
# ./gromacs_commands.sh <pdb> <ionic_strength>

# Before running this script, do the following stuff for each protein:
# 1. Double check if the PDB file has the right pH and temperature information
# 2. Copy .mdp files to the working directory
# 3. Modify the .mdp files according to the experimental conditions

pdb=$1
ionic_strength=$2

python3 /home/dynamics/akshatha/Analysis_codes/initial_files/process_pdb.py /home/dynamics/akshatha/Simfiles/$pdb/$pdb.pdb
cp -r charmm36m.ff $pdb

# -ter NH3+ and COO-

/usr/local/gromacs/bin/gmx  pdb2gmx -f $pdb/$pdb -o $pdb/processed.gro -ignh -ter

# /usr/local/gromacs/bin/gmx  pdb2gmx -f $pdb/$pdb -o $pdb/processed.gro -ignh -ter <<EOF
# 1
# 1
# 0
# 0
# EOF 

# Move these files first 
mv posre.itp $pdb
mv topol.top $pdb

/usr/local/gromacs/bin/gmx  editconf -f $pdb/processed.gro -o $pdb/newbox.gro -c -d 1.2 -bt cubic

/usr/local/gromacs/bin/gmx  solvate -cp $pdb/newbox.gro -o $pdb/solv.gro -p $pdb/topol.top

/usr/local/gromacs/bin/gmx  grompp -f ions.mdp -c $pdb/solv.gro -p $pdb/topol.top -o $pdb/ions.tpr

printf "13" | /usr/local/gromacs/bin/gmx  genion -s $pdb/ions.tpr -o $pdb/solv_ions.gro -p $pdb/topol.top -pname SOD -nname CLA -neutral -conc $ionic_strength

/usr/local/gromacs/bin/gmx  grompp -f minim.mdp -c $pdb/solv_ions.gro -p $pdb/topol.top -o $pdb/em.tpr

/usr/local/gromacs/bin/gmx  mdrun -v -deffnm $pdb/em

printf "11 0" | /usr/local/gromacs/bin/gmx  energy -f $pdb/em.edr -o $pdb/potential.xvg

/usr/local/gromacs/bin/gmx  grompp -f $pdb/nvt.mdp -c $pdb/em.gro -r $pdb/em.gro -p $pdb/topol.top -o $pdb/nvt.tpr

/usr/local/gromacs/bin/gmx  mdrun -deffnm $pdb/nvt

printf "16 0" | /usr/local/gromacs/bin/gmx  energy -f $pdb/nvt.edr -o $pdb/temperature.xvg

/usr/local/gromacs/bin/gmx  grompp -f $pdb/npt.mdp -c $pdb/nvt.gro -r $pdb/nvt.gro -t $pdb/nvt.cpt -p $pdb/topol.top -o $pdb/npt.tpr

/usr/local/gromacs/bin/gmx  mdrun -deffnm $pdb/npt

printf "17 0" | /usr/local/gromacs/bin/gmx  energy -f $pdb/npt.edr -o $pdb/pressure.xvg

printf "23 0" | /usr/local/gromacs/bin/gmx  energy -f $pdb/npt.edr -o $pdb/density.xvg

/usr/local/gromacs/bin/gmx  grompp -f $pdb/md.mdp -c $pdb/npt.gro -p $pdb/topol.top -o $pdb/md.tpr 

/usr/local/gromacs/bin/gmx  mdrun -v -deffnm $pdb/md

# taskset -c 35-162 /usr/local/gromacs/bin/gmx mdrun -v -s $pdb/md.tpr -cpi $pdb/md.cpt -deffnm $pdb/md -append -nt 128









# constraints = h-bonds
# cutoff-scheme = Verlet
# coulombtype = PME

# # Add these four lines
# vdwtype                 = cutoff
# vdw-modifier            = force-switch
# rlist                   = 1.2
# rvdw-switch             = 1.0

# Change these three quantities
# rvdw = 1.2
# rcoulomb = 1.2
# DispCorr = no
