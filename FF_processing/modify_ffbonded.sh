# modify_ffbonded.sh
# Akshatha 
# 17 July 2023

# Has code for tagging combinations for bonds, angles and dihedrals

# Usage: ./modify_ffbonded.sh input_file output_file param_type

# Copy the block of parameters (BONDS, ANGLES, DIHEDRALS) to an input file
# Copy back the results into the corresponding location of the parameter blocks in the force field file 
# NOTE: Dihedral type 9 entries were sorted before adding to the output file to avoid GROMACS errors of duplicate block of parameters

#!/bin/bash

input_file=$1
output_file=$2
param_type=$3

if [ "$param_type" == "ANGLES" ]; then

  awk '{if ($0 ~ /^;/ || $0 ~ /^\[/ || $0 ~ /^$/) {print $0} else { 
      
      {printf "H-%-2s  H-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
      {printf "S-%-2s  S-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

      if ($1 != $2 && $2 != $3 && $1 != $3)
      {
        {printf "H-%-2s  %-4s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  %-4s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "H-%-2s  H-%-2s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  S-%-2s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "%-4s  %-4s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "%-4s  %-4s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "%-4s  H-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "%-4s  S-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "H-%-2s  S-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  H-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "H-%-2s  H-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  S-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
      }

      else if (($1 == $2 && $2 != $3) || ($2 == $3 && $1 != $2))
      {
        {printf "H-%-2s  %-4s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  %-4s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "H-%-2s  H-%-2s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  S-%-2s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "%-4s  %-4s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "%-4s  %-4s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "%-4s  H-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "%-4s  S-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        
        {printf "H-%-2s  S-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  H-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "H-%-2s  H-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  S-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
      } 

      else if ($1 == $2 && $2 == $3)
      {
        {printf "H-%-2s  %-4s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  %-4s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "H-%-2s  H-%-2s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  S-%-2s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        
        {printf "H-%-2s  S-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  H-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
      }

      else if ($1 == $3 && $2 != $3)
      {
        {printf "H-%-2s  %-4s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  %-4s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

        {printf "H-%-2s  H-%-2s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  S-%-2s  %-4s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
      
        {printf "H-%-2s  S-%-2s  S-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}
        {printf "S-%-2s  H-%-2s  H-%-2s     %d   %3.3f    %3.3f \n", $1, $2, $3, $4, $5, $6}

      }
      }}'  $input_file > $output_file

elif [ "$param_type" == "BONDS" ]; then
  awk '{if ($0 ~ /^;/ || $0 ~ /^\[/ || $0 ~ /^$/) {print $0} else { 
        {printf "H-%-2s H-%-2s      %d    %1.5f   %6.1f \n", $1, $2, $3, $4, $5}
        {printf "S-%-2s S-%-2s      %d    %1.5f   %6.1f \n", $1, $2, $3, $4, $5}    

        if ($1 == $2)
        {
          {printf "H-%-2s S-%-2s      %d    %1.5f   %6.1f \n", $1, $2, $3, $4, $5}
        }
        else
        {
          {printf "S-%-2s H-%-2s      %d    %1.5f   %6.1f \n", $1, $2, $3, $4, $5}
          {printf "H-%-2s S-%-2s      %d    %1.5f   %6.1f \n", $1, $2, $3, $4, $5}
        }
        }}'  $input_file > $output_file

elif [ "$param_type" == "DIHEDRALS" ]; then
  awk '{if ($0 ~ /^;/ || $0 ~ /^\[/ || $0 ~ /^$/) {print $0} else {

    {printf "H-%-2s H-%-2s H-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "S-%-2s S-%-2s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    
    if ($1 == $2 && $2 == $3 && $3 == $4)
    {
      {printf "H-%-2s %-4s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s %-4s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

      {printf "H-%-2s H-%-2s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s S-%-2s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

      {printf "H-%-2s H-%-2s H-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s S-%-2s S-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

      {printf "H-%-2s H-%-2s H-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s S-%-2s S-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "H-%-2s H-%-2s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

    }

     if ($1 == $4 && $2 == $3 && $1 != $2)
    {
      {printf "H-%-2s %-4s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s %-4s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

      {printf "H-%-2s H-%-2s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s S-%-2s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

      {printf "H-%-2s H-%-2s H-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s S-%-2s S-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

      {printf "H-%-2s H-%-2s H-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s S-%-2s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

      {printf "H-%-2s H-%-2s H-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "S-%-2s S-%-2s S-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
      {printf "H-%-2s H-%-2s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

      {printf "%-4s H-%-2s S-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8} 
    }

    if ($1 == $4 && $2 != $3 && $1 != $2)
    {
    {printf "%-4s S-%-2s H-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "%-4s H-%-2s S-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    }

    else 
    {
    {printf "H-%-2s %-4s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "S-%-2s %-4s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "%-4s %-4s %-4s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "%-4s %-4s %-4s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

    {printf "H-%-2s H-%-2s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "S-%-2s S-%-2s %-4s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "%-4s %-4s H-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "%-4s %-4s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

    {printf "H-%-2s H-%-2s H-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "S-%-2s S-%-2s S-%-2s %-4s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "%-4s H-%-2s H-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "%-4s S-%-2s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

    {printf "H-%-2s H-%-2s H-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "S-%-2s S-%-2s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}

    {printf "H-%-2s H-%-2s H-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "S-%-2s H-%-2s H-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "S-%-2s S-%-2s S-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "H-%-2s S-%-2s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "H-%-2s H-%-2s S-%-2s S-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    {printf "S-%-2s S-%-2s H-%-2s H-%-2s  %d   %6.2f    %8.5f     %d \n", $1, $2, $3, $4, $5, $6, $7, $8}
    }      
  }}'  $input_file > $output_file

  python3 sort_lines.py $output_file
fi

