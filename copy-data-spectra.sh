#! /usr/bin/env bash

if [ $# -eq 0 ]
then
  echo "Usage: copy-data-spectra.sh source destination"
elif [ $# -ge 2 ]
then
  sourceDir=$1
  targetDir=$2
elif [ $# -ge 1 ]
then
  sourceDir=`pwd`
  targetDir=$1
fi

mkdir $targetDir/spectra
cp $sourceDir/*.log $targetDir/spectra
cp -R $sourceDir/results/v2data.dat $targetDir/spectra/
cp -R $sourceDir/results/v2data-inte.dat $targetDir/spectra/
cp -R $sourceDir/results/phipspectra.dat $targetDir/spectra/
cp -R $sourceDir/results/numbers.dat $targetDir/spectra/
cp -R $sourceDir/results/names.dat $targetDir/spectra/
cp -R $sourceDir/results/masses.dat $targetDir/spectra/
cp -R $sourceDir/results/momentum.dat $targetDir/spectra/


# clean up to save space...
#rm $targetDir/results/decdat2.dat
#rm $targetDir/results/decdat_mu.dat
#rm $targetDir/spectra/decdat2.dat
#rm $targetDir/spectra/decdat_mu.dat


rm *.log
