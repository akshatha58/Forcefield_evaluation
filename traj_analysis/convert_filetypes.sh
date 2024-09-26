# convert_filetypes.sh
# Akshatha
# 24 September 2024

# One unified program that converts .gro to .gro without tags to .pdb files

# proteins=("1ozi" "1gmo" "1r36" "1z9b" "2ec7" "2gmo" "2juo" "2jxy" "2k0q" "2l4x" "2l6b" "2lro" "2m68" "4bwh" "7ywr")

pdb=$1

path=/home/dynamics/akshatha/Final_simfiles/$pdb
modify_grofile="modify_grofile.sh"
mod_grofile_path="/home/dynamics/akshatha/Final_simfiles/FF_analysis_codes/modify_grofile.sh"

cd $path

# Convert the xtc file to gro file format (without solvent, every 1 ns frame)
/usr/local/gromacs/bin/gmx trjconv -f md.xtc -s md.tpr -o traj.gro -pbc mol -center -dt 1000

if [ -f "$modify_grofile" ]; then
    # modify the gro file to get untagged gro file
    ./modify_grofile.sh traj.gro trajmod.gro
else
    cp $mod_grofile_path .
    chmod +x $modify_grofile
    ./modify_grofile.sh traj.gro trajmod.gro
fi   

# Get the trajectory information in pdb format
/usr/local/gromacs/bin/gmx trjconv -f trajmod.gro -s trajmod.gro -o $pdb"_trajmod.pdb"



