# modify_gro.sh
# Akshatha
# 16 August 2024

# Modifies the gro file to remove H- and S- tags from the residues
# Usage: ./modify_gro.sh <input_gro_file> <output_gro_file>

grofile=$1
cp $grofile $2

# use sed to remove H- and S- from grofile
sed -i 's/H-//g' $2
sed -i 's/S-//g' $2

#  print all lines uniformly according to the .gro format with spacings only if the first field does NOT contain "SOL" or "CL" or "NA"

# # replace the first two lines of temp.gro with the original lines in the gro file 
# head -n 2 $2 > temp.gro 

awk '{
if ($1 !~ /SOL/ && $1 !~ /CL/ && $1 !~ /NA/ && NF == 6 && $1 !~ /Protein/) 
        printf "%8s  %5s%5d   %1.3f   %1.3f   %1.3f\n", $1, $2, $3, $4, $5, $6
else
        print $0
}' $2 >> temp.gro

# delete last line of temp.gro
# sed -i '$ d' temp.gro

# grep the line for the first occurrence of "SOL" and cat the tail of the file from that line to end of file to temp.gro
# grep -n -m 1 "SOL" $2 | cut -d: -f1 | xargs -I {} tail -n +{} $2 >> temp.gro

# rename temp.gro to mod_npt.gro and delete temp.gro
mv temp.gro $2

