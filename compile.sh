#!/usr/bin/env bash

mkdir -p build
cd build
rm -fr *
CXX=g++-12 cmake ..
make -j4
make install
cd ..
rm -fr build
