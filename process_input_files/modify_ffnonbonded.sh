# modify_ffnonbonded.sh
# Akshatha
# 20 July 2024

# Adds a superset of H- and S- atom type information to atomtypes.atp
# Adds tagged atom information to ffnonbonded.itp
# Note: this is primarily for AMBER force fields, where LJ interaction parameters are not explicitly included in the force field

# Usage: ./modify_ffnonbonded.sh atomtypespath ffnonbondedpath

atomtypespath=$1
ffnonbondedpath=$2

atomtypes=$(basename $atomtypespath)
ffnonbonded=$(basename $ffnonbondedpath)


################################################### UPDATE TAGGED ATOMS TO ATOMTYPES.ATP ####################################################

# Create a backup of the original atomtypes.atp
cp $atomtypespath $atomtypes'_backup.atp'

# for each line, print H-%s and S-%s (%s= $0)
awk '{print "H-"$0}' $atomtypes'_backup.atp' >> $atomtypespath
awk '{print "S-"$0}' $atomtypes'_backup.atp' >> $atomtypespath

################################################### UPDATE TAGGED ATOMS TO FFNONBONDED.ITP ####################################################

cp $ffnonbondedpath $ffnonbonded'_backup.itp'
# for each line, if the line starts with ; or [  or empty, print as is, else, print H-%s and S-%s (%s= $0)
awk '{if ($0 ~ /^;/ || $0 ~ /^\[/ || $0 ~ /^$/) {print $0} else {print $0"\nH-"$0"\nS-"$0}}' $ffnonbondedpath >> $ffnonbonded'_tag.itp'

# if the line starts with "H-" or "S-", delete two spaces before the second field
awk '{if ($0 ~ /^H-/ || $0 ~ /^S-/) {sub(/  /,"",$0); print $0} else {print $0}}' $ffnonbonded'_tag.itp' > $ffnonbondedpath

rm $ffnonbonded'_tag.itp'
mv $ffnonbonded'_backup.itp' $(dirname $ffnonbondedpath)
mv $atomtypes'_backup.atp' $(dirname $atomtypespath)