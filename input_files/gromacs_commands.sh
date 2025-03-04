# gromacs_commands.sh
# Akshatha
# 24 September 2024

# GROMACS commands for running simulations

# Before running this script, do the following for each protein:
# 1. Double check if the PDB file has the right pH and temperature information
# 2. Copy .mdp files to the working directory
# 3. Modify the .mdp files according to the experimental conditions

# ./gromacs_commands.sh <pdbpath> <ionic strength>

pdbpath=$1
pdb=$(basename $pdbpath)
ionic_strength=$2
GROMACS=/usr/local/gromacs/bin/gmx

python3 process_pdb.py $pdbpath/$pdb.pdb

# Copy this if CHARMM36m is used
# cp -r charmm36m.ff $pdbpath

 $GROMACS pdb2gmx -f $pdbpath/$pdb.pdb -o $pdbpath/processed.gro -ignh -ter

# Move these files first 
mv posre.itp $pdbpath
mv topol.top $pdbpath

$GROMACS  editconf -f $pdbpath/processed.gro -o $pdbpath/newbox.gro -c -d 1.2 -bt cubic

$GROMACS solvate -cp $pdbpath/newbox.gro -o $pdbpath/solv.gro -p $pdbpath/topol.top

$GROMACS grompp -f ions.mdp -c $pdbpath/solv.gro -p $pdbpath/topol.top -o $pdbpath/ions.tpr

# Amber
printf "13" | $GROMACS genion -s $pdbpath/ions.tpr -o $pdbpath/solv_ions.gro -p $pdbpath/topol.top -pname  NA -nname CL -neutral -conc $ionic_strength
# CHARMM
# printf "13" | $GROMACS genion -s $pdbpath/ions.tpr -o $pdbpath/solv_ions.gro -p $pdbpath/topol.top -pname  SOD -nname CLA -neutral -conc $ionic_strength

$GROMACS  grompp -f minim.mdp -c $pdbpath/solv_ions.gro -p $pdbpath/topol.top -o $pdbpath/em.tpr

taskset -c 1-128 $GROMACS  mdrun -v -deffnm $pdbpath/em -nt 125

printf "10 0" | $GROMACS  energy -f $pdbpath/em.edr -o $pdbpath/potential.xvg

$GROMACS  grompp -f $pdbpath/nvt.mdp -c $pdbpath/em.gro -r $pdbpath/em.gro -p $pdbpath/topol.top -o $pdbpath/nvt.tpr

taskset -c 1-128 $GROMACS  mdrun -deffnm $pdbpath/nvt -nt 125

printf "16 0" | $GROMACS  energy -f $pdbpath/nvt.edr -o $pdbpath/temperature.xvg

$GROMACS  grompp -f $pdbpath/npt.mdp -c $pdbpath/nvt.gro -r $pdbpath/nvt.gro -t $pdbpath/nvt.cpt -p $pdbpath/topol.top -o $pdbpath/npt.tpr

taskset -c 1-128 $GROMACS mdrun -deffnm $pdbpath/npt -nt 125

printf "18 0" | $GROMACS energy -f $pdbpath/npt.edr -o $pdbpath/pressure.xvg

printf "24 0" | $GROMACS energy -f $pdbpath/npt.edr -o $pdbpath/density.xvg

$GROMACS grompp -f $pdbpath/md.mdp -c $pdbpath/npt.gro -p $pdbpath/topol.top -o $pdbpath/md.tpr 

taskset -c 1-128 $GROMACS  mdrun -v -deffnm $pdbpath/md -nt 125

# taskset -c 1-128 $GROMACS -v -s $pdbpath/md.tpr -cpi $pdbpath/md.cpt -deffnm $pdbpath/md -append -nt 125





# # constraints = h-bonds
# # cutoff-scheme = Verlet
# # coulombtype = PME

# # # Add these four lines
# # vdwtype                 = cutoff
# # vdw-modifier            = force-switch
# # rlist                   = 1.2
# # rvdw-switch             = 1.0

# # Change these three quantities
# # rvdw = 1.2
# # rcoulomb = 1.2
# # DispCorr = no
