# parse_J_HNHA.sh
# Akshatha
# 3 September 2023

# Parses J coupling data from STAR format
# Saves residue number, residue name, J coupling constant value and error 

# Usage: ./parse_J_HNHA.sh 1r36 

pdb=$1
file=/home/dynamics/akshatha/NMR_BMRB/J_coupling/$pdb"_J_HNHA"
savepath=/home/dynamics/akshatha/Final_simfiles/$pdb/NMR_calc

declare -a residues  
declare -a jconstants 
declare -a err        
declare -a resname    


echo "Reading file..."


while IFS= read -r line
do
    residues+=($(awk 'NF == 39 {
        print $6
    }' <<< "$line"))

    jconstants+=($(awk 'NF == 39 {
        print $23
    }' <<< "$line"))

    err+=($(awk 'NF == 39 {
        print $26
    }' <<< "$line"))

    resname+=($(awk 'NF == 39 {
        print $8
    }' <<< "$line"))

done < $file


# Since mdtraj counts residues from zero, add an offset of one to residues
echo "Adjusting residue numbers..."

for i in "${!residues[@]}"; do
    residues[$i]=$((${residues[$i]} - 1))
done

# echo "Residues: "
# echo ${residues[@]}

# echo "J-constants: "
# echo ${jconstants[@]}

echo "Saving extracted data to file..."

# Save values in a text file

printf "%s\n" "${residues[@]}" > temp1.txt
printf "%s\n" "${jconstants[@]}" > temp2.txt
printf "%s\n" "${err[@]}" > temp3.txt
printf "%s\n" "${resname[@]}" > temp4.txt
printf "$s\n"

paste temp1.txt temp4.txt temp2.txt temp3.txt > $savepath"/j_exp.txt"
rm temp1.txt temp2.txt temp3.txt temp4.txt

echo "Done. Output saved to $savepath"/j_exp.txt""

