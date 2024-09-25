# modify_ffnonbonded.sh
# Akshatha
# 20 July 2024

# Adds a superset of H- and S- atom type information to atomtypes.atp
# Adds tagged atom information to ffnonbonded.itp
# Note: this is primarily for AMBER force fields, where LJ interaction parameters are not explicitly included in the force field

# Usage: ./modify_ffnonbonded.sh atomtypes.atp ffnonbonded.itp

atomtypes=$1
ffnonbonded=$2

################################################### UPDATE TAGGED ATOMS TO ATOMTYPES.ATP ####################################################

# Create a backup of the original atomtypes.atp
cp $atomtypes'.atp' $atomtypes'_backup.atp'

# for each line, print H-%s and S-%s (%s= $0)
awk '{print "H-"$0}' $atomtypes'_backup.atp' >> $atomtypes'.atp'
awk '{print "S-"$0}' $atomtypes'_backup.atp' >> $atomtypes'.atp'

################################################### UPDATE TAGGED ATOMS TO FFNONBONDED.ITP ####################################################

cp $ffnonbonded'.itp' $ffnonbonded'_backup.itp'
# for each line, if the line starts with ; or [  or empty, print as is, else, print H-%s and S-%s (%s= $0)
awk '{if ($0 ~ /^;/ || $0 ~ /^\[/ || $0 ~ /^$/) {print $0} else {print $0"\nH-"$0"\nS-"$0}}' $ffnonbonded'.itp' >> $ffnonbonded'_tag.itp'

# if the line starts with "H-" or "S-", delete two spaces before the second field
awk '{if ($0 ~ /^H-/ || $0 ~ /^S-/) {sub(/  /,"",$0); print $0} else {print $0}}' $ffnonbonded'_tag.itp' > $ffnonbonded'.itp'

rm $ffnonbonded'_tag.itp'
