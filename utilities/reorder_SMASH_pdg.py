#!/usr/bin/env python

import numpy as np

f = open("pdg-SMASH_orig.dat")
SMASH_table = f.read().split("\n")

mass_list   = []
header_list = []
decay_list  = []

iline = 0
while iline < len(SMASH_table) - 1:
    header_list.append(SMASH_table[iline])
    mass = float(SMASH_table[iline].split()[2])
    mass_list.append(mass)
    num_decay_channels = int(SMASH_table[iline].split()[-1])
    decay_block = []
    for i in range(num_decay_channels):
        iline += 1
        decay_block.append(SMASH_table[iline])
    decay_list.append(decay_block)
    iline += 1

mass_list = np.array(mass_list)
indx_list = np.argsort(mass_list)

f = open("pdg-SMASH.dat", "w")
for idx in indx_list:
    f.write("{0}\n".format(header_list[idx]))
    for iline in range(len(decay_list[idx])):
        f.write("{0}\n".format(decay_list[idx][iline]))
f.close()
