#! /usr/bin/env bash

# Folders containing $already_exist will be skipped.

# Example:
# ./spectraFromHydroResults.sh ~/Downloads/test/ ./azSpectraV1.2.1.e 2

already_exist="spectra/v2data.dat"

# check argument and get directories
if [ $# -le 1 ]
then
  echo "Usage: spectraFromHydroResults.sh hydro_results_folder spectra_executable and parameters"
  exit
else
  hydroDir=$1
fi

current_dir=`pwd`

# get hydro_executable and parameters
executable=""
while [ $# -gt 1 ]
do
  arg=`echo "$2" | sed 's/ /\?/g'`
  executable="$executable $arg"
  shift
done

for sub_folder in `ls $hydroDir`
do
  if [ -e $hydroDir/$sub_folder/$already_exist ]
  then
    echo Skipping $sub_folder
  else
    echo Dealing with $sub_folder ...
    python ./record-spectra.py $hydroDir/$sub_folder $executable
  fi
done
