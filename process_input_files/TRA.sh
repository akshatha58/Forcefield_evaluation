#  TAG RIGID ATOMS (TRA) 
#  Akshatha
#  17 June 2024

# Gives H- or S- tags to residues and atom types present in the helix or sheet regions
# DSSP program used to analyse the corresponding PDB file for secondary structure elements 
# Tagging directly done to the topology file 

# USAGE: ./TRA.sh pdbfilepath topfilepath
######################################################### SECONDARY STRUCTURE ASSIGNMENT USING DSSP ###########################################################################################################

pdbfilepath=$1  
topfilepath=$2 

# Extracting PDBID from the path
pdbfile=$(basename $pdbfilepath)
pdbfile="${pdbfile%.*}"

# Extracting top file name from the path
topfile=$(basename $topfilepath)
topfile="${topfile%.*}"

echo " "
echo "Processing" $pdbfile "..."

# Adds element information to the last column
# Checks for discontinuous residues and residues not starting with 1, updates offset file with this information
python3 process_pdb.py $pdbfilepath

echo "Analysing secondary structure ..."

# DSSP Command for secondary structure assignment
mkdssp --mmcif-dictionary /usr/local/share/libcifpp/mmcif_pdbx_v50.dic $pdbfilepath $pdbfile".mmcif"

####################################################### EXTRACT HELIX AND SHEET RESIDUES FROM mmCIF FILE ######################################################################################################

declare -a helix_residues
declare -a sheet_residues

# for every entry in the mmcif file that has number of fields as 25, 
# if the 5th field has H,G or I, add the third field to helix_residues
# if the 5th field has E, add the third field to sheet_residues

echo "Extracting helix and sheet residues ..."

read -r offset < $pdbfile"_offsets.txt"

while IFS= read -r line
do
    helix_residues+=($(awk 'NF == 25 {
        if ($5 == "H" || $5 == "G" || $5 == "I") {
            print $3
        }
    }' <<< "$line"))
    
    sheet_residues+=($(awk 'NF == 25 {
        if ($5 == "E") {
            print $3
        }
    }' <<< "$line"))
done < $pdbfile".mmcif"

#################################################### ADDING OFFSETS TO THE RESIDUE NUMBERS #######################################################################################################################
# append Helices and corresponding residue list to the second line of PDBID_offsets.txt

echo "Helices " ${helix_residues[@]} >> $pdbfile"_ss.txt"
echo "Sheets " ${sheet_residues[@]} >> $pdbfile"_ss.txt"

echo "Adding offsets ..."

# Add offset to helix_residues
for i in "${!helix_residues[@]}"; do
    helix_residues[$i]=$((${helix_residues[$i]} + offset))
done

# Add offset to sheet_residues
for i in "${!sheet_residues[@]}"; do
    sheet_residues[$i]=$((${sheet_residues[$i]} + offset))
done


# print helix_residues onto terminal
echo " "
echo "Helices: "
echo ${helix_residues[@]}

# print helix_residues onto terminal
echo "Sheets: "
echo ${sheet_residues[@]}

# append Helices and corresponding residue list to the second line of PDBID_offsets.txt
echo "Helices " ${helix_residues[@]} >> $pdbfile"_offsets.txt"
echo "Sheets " ${sheet_residues[@]} >> $pdbfile"_offsets.txt"

###################################################### TAG HELICES AND SHEETS IN .TOP FILE ###################################################################################################################
echo " "
echo "Updating topology file ..."

awk -v s_num="$(IFS=,; echo "${sheet_residues[*]}")" -v h_num="$(IFS=,; echo "${helix_residues[*]}")" '
BEGIN {split(s_num, a, ",")} {split(h_num, b, ",")}
NF ge 8 {
    updated = 0
    for (i in a) {
        if ($3 == a[i] && $2 ~ /[a-zA-Z]/) {
            if (length($2) == 1){
                head = substr($0, 1, 14)
                mid1 = substr($0, 15, 2)
                mid = substr($0, 17, 10)
                mid2 = substr($0, 27, 2)
                end = substr($0, 29)

                # substitute mid with S-
                printf "%sS-%sS-%s\n", head, mid, end
                updated = 1
            }
            if (length($2) == 2){
                head = substr($0, 1, 13)
                mid1 = substr($0, 14, 2)
                mid = substr($0, 16, 11)
                mid2 = substr($0, 27, 2)
                end = substr($0, 29)

                # substitute mid with S-
                printf "%sS-%sS-%s\n", head, mid, end
                updated = 1
            }
        }
    }

    for (i in b) {
        if ($3 == b[i] && $2 ~ /[a-zA-Z]/) {
            if (length($2) == 1){
              head = substr($0, 1, 14)
                mid1 = substr($0, 15, 2)
                mid = substr($0, 17, 10)
                mid2 = substr($0, 27, 2)
                end = substr($0, 29)

                # substitute mid with H-
                printf "%sH-%sH-%s\n", head, mid, end
                updated = 1
            }
            if (length($2) == 2){
               head = substr($0, 1, 13)
                mid1 = substr($0, 14, 2)
                mid = substr($0, 16, 11)
                mid2 = substr($0, 27, 2)
                end = substr($0, 29)

                # substitute mid with H-
                printf "%sH-%sH-%s\n", head, mid, end
                updated = 1
            }
        }
    }

    if (updated == 0) {
        print
    }

}' $topfilepath > $topfile"_tagged.top"

# Move offset file, mmcif file and tagged top file to the same directory as the input files
mv $pdbfile"_offsets.txt" $pdbfile"_ss.txt" $pdbfile".mmcif" $topfile"_tagged.top" $(dirname $topfilepath)

echo "Done!"
