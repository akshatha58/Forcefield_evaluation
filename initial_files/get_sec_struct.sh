#  get_sec_struct.sh
#  Akshatha
#  5 March 2025

# Update of TRA_v3. Updates the secondary structure files according to the original labels in the mmcif files. 

mmcif_path=/home/dynamics/akshatha/mmcif_files
pdb=$1

######################################################## EXTRACT HELIX AND SHEET RESIDUES FROM mmCIF FILE #############################################################################

declare -a alpha_helix
declare -a beta_bridge_isolated
declare -a beta_ladder
declare -a threeten_helix
declare -a kappa_helix
declare -a pi_helix
declare -a hbonded_turn
declare -a bends
declare -a helix_termini

echo "Extracting secondary structure information ..."

while IFS= read -r line
do
    alpha_helix+=($(awk 'NF == 25 {
        if ($5 == "H") {
            print $3
        }
    }' <<< "$line"))

    beta_bridge_isolated+=($(awk 'NF == 25 {
        if ($5 == "B") {
            print $3
        }
    }' <<< "$line"))

    beta_ladder+=($(awk 'NF == 25 {
        if ($5 == "E") {
            print $3
        }
    }' <<< "$line"))
    
    threeten_helix+=($(awk 'NF == 25 {
        if ($5 == "G") {
            print $3
        }
    }' <<< "$line"))

    pi_helix+=($(awk 'NF == 25 {
        if ($5 == "I") {
            print $3
        }
    }' <<< "$line"))

    hbonded_turn+=($(awk 'NF == 25 {
        if ($5 == "T") {
            print $3
        }
    }' <<< "$line"))

    bends+=($(awk 'NF == 25 {
        if ($5 == "S") {
            print $3
        }
    }' <<< "$line"))

        kappa_helix+=($(awk 'NF == 25 {
        if ($5 == "P") {
            print $3
        }
    }' <<< "$line"))

done < $mmcif_path/$pdb".mmcif"

# ################################################### SAVE SS INFORMATION ########################################################################

echo "Appending information into "$pdb"_sslists.txt ..."

echo ${alpha_helix[@]} >> $mmcif_path/$pdb"_sslists.txt"
echo ${threeten_helix[@]} >> $mmcif_path/$pdb"_sslists.txt"
echo ${pi_helix[@]} >> $mmcif_path/$pdb"_sslists.txt"
echo ${beta_bridge_isolated[@]} >> $mmcif_path/$pdb"_sslists.txt"
echo ${beta_ladder[@]} >> $mmcif_path/$pdb"_sslists.txt"
echo ${hbonded_turn[@]} >> $mmcif_path/$pdb"_sslists.txt"
echo ${bends[@]} >> $mmcif_path/$pdb"_sslists.txt"
echo ${kappa_helix[@]} >> $mmcif_path/$pdb"_sslists.txt"
# Last line: helix termini residues


