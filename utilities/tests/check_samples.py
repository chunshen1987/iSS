#!/usr/bin/env python

from numpy import *

sample_file = open("OSCAR.DAT", "r")

# head the header
for i in range(3):
    dummy = sample_file.readline()

nev = 0
E = 0.
px = 0.
py = 0.
pz = 0.
while True:
    line = sample_file.readline()
    if line == '': break
    nev += 1
    n_samples = int(line.split()[1])
    for i in range(n_samples):
        data = sample_file.readline().split()
        px += float(data[2])
        py += float(data[3])
        pz += float(data[4])
        E += float(data[5])

print("nev = {0}".format(nev))
print("<px> = {0:.4f} GeV".format(px/nev))
print("<py> = {0:.4f} GeV".format(py/nev))
print("<pz> = {0:.4f} GeV".format(pz/nev))
print("<E> = {0:.4f} GeV".format(E/nev))

