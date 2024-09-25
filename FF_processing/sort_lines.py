"""
sort_lines.py
Akshatha
20 July 2024

Sorts tagged dihedral type 9 alphabetically, to keep block of similar parameters together
Outputs to sorted.txt, the corresponding lines can then be copied to the respective location in ffbonded.itp

Usage: python3 sort_lines.py inputfile.txt

"""
import sys
inputfile=sys.argv[1]

with open(inputfile, 'r') as file:
    lines = file.readlines()

# Sort the lines
lines.sort()

# Write the sorted and deduplicated lines back to a file
with open('sorted.txt', 'w') as sorted_file:
    previous_line = None
    for line in lines:
        if line != previous_line:
            sorted_file.write(line)
            previous_line = line
