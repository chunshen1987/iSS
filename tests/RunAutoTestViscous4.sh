#!/usr/bin/env bash

cp tests/music_input_14 tests/music_input

./iSS.e tests/iSS_parameters.dat tests testViscousOneFluidCell4.dat

python3 tests/TestOutputFiles.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
    rm -fr ./check*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
