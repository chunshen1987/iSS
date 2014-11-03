#! /usr/bin/env bash
icc test.cpp arsenal.cpp RandomVariable2DArray.cpp RandomVariable.cpp TableFunction.cpp Table.cpp -o test.e -fast -B/usr/lib/i386-linux-gnu -I/usr/include/i386-linux-gnu
