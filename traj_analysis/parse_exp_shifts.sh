# parse_exp_shifts.sh
# Akshatha
# 12 September 2023

# Parses chemical shift data from STAR format
# Residue number, residue name, atom type, chemical shift type (1, 13, 15), chemical shift value, error, ambiguity

# Usage: ./parse_exp_shifts.sh

# proteins=("1r36" "1ozi" "2gmo" "2juo" "2k0q" "2l4x" "2ec7" "2jxy" "2l6b" "2lro" "2m68" "4bwh")
# 7ywr: number of fields is 26 and 1z9b: number of fields is 25
proteins=("7ywr")

for pdb in "${proteins[@]}"
do
    file=/home/dynamics/akshatha/NMR_BMRB/shifts/$pdb"_shifts.txt"
    savepath=/home/dynamics/akshatha/Final_simfiles/$pdb/NMR_calc

    echo "Reading file..."

    awk '{if (NF == 26)
        # printf "%4d  %3s  %-5s  %2d  %10.4f     %2.2f     %d\n", $19, $7, $8, $10, $11, $12, $14
        printf "%4d  %3s  %-5s  %2d  %10.4f     %2.2f\n", $21, $8, $9, $11, $12, $13
    }' $file >> $savepath"/"$pdb"_shifts.txt" 

    echo "Done. File read to $savepath"/"$pdb"_shifts.txt"."
done
