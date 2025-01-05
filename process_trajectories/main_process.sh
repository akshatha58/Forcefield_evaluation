# main.sh
# Akshatha 
# 6 October 2024

# Main codes that implement workflow (See README)

# give permissions (first time execution)
chmod +x convert_fileformats.sh
chmod +x modify_grofile.sh
chmod +x process_traj.py
chmod +x visualise.sh

# ./convert_fileformats.sh <directory>

# Standard paths
sourcepath=$(pwd)

# Input protein and path to store data
pdbpath=$1
pdb=$(basename $pdbpath)
pdb="${pdb%.*}"

echo $pdb

cd $sourcepath

./convert_fileformats.sh $pdb $pdbpath
./visualise.sh $pdbpath
