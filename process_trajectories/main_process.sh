# main.sh
# Akshatha 
# 6 October 2024

# Main codes that implement workflow (See README)

# give permissions (first time execution)
chmod +x convert_fileformats.sh
chmod +x convert_fileformats_charmm.sh
chmod +x modify_grofile.sh
chmod +x process_traj.py
chmod +x visualise.sh

# ./convert_fileformats.sh <directory> <ff>
# ff: amber or charmm

# Standard paths
sourcepath=$(pwd)

# Input protein and path to store data
pdbpath=$1
pdb=$(basename $pdbpath)
pdb="${pdb%.*}"
ff=$2

echo $pdb
cd $sourcepath

if [ "$ff" == 'amber' ] 
then
    ./convert_fileformats.sh $pdb $pdbpath
elif [ "$ff" == 'charmm' ]
then
    ./convert_fileformats_charmm.sh $pdb $pdbpath
fi

# ./visualise.sh $pdbpath
