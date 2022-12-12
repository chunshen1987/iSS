#!/usr/bin/env bash

./iSS.e tests/iSS_parameters_ideal.dat tests testIdealOneFluidCell1.dat

python3 tests/TestOutputFiles.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
    rm -fr ./check*.dat
else
    echo "Tests FAILED :("
    exit 1
fi
