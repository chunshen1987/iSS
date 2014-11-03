#! /usr/bin/env bash

if [ $# -eq 0 ]
then
  echo "Usage: prepare-spectra.sh source [destination]"
elif [ $# -ge 2 ]
then
  sourceDir=$1
  targetDir=$2
elif [ $# -ge 1 ]
then
  sourceDir=$1
  targetDir=`pwd`
fi

rm $targetDir/results/*
cp -R $sourceDir/results/* $targetDir/results
