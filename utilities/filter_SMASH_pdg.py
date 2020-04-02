#!/usr/bin/env python

import numpy as np
import sys


if len(sys.argv) < 2:
    print("{} filename".format(sys.argv[0]))
    exit()

filename = str(sys.argv[1])
f = open(filename, "r")
SMASH_table = f.read().split("\n")
f.close()

pdg_list = []
mass_list   = []
header_list = []
decay_list  = []

iline = 0
while iline < len(SMASH_table) - 1:
    pdg_list.append(SMASH_table[iline].split()[0])
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

pdg_list_wbbar = []
f = open("pdg-SMASH.dat", "w")
for idx in indx_list:
    ipart = int(pdg_list[idx])
    pdg_code_dissected = [int(d) for d in str(abs(ipart))]

    # here are filters
    if len(pdg_code_dissected) < 3:
        # photon
        continue
    if pdg_code_dissected[-2] > 3 or pdg_code_dissected[-3] > 3:
        # heavy quark mesons or baryons
        continue
    if len(pdg_code_dissected) > 3 and pdg_code_dissected[-4] > 3:
        # heavy quark baryons
        continue

    # passed all the filters, add it to the final pdg list
    pdg_list_wbbar.append(ipart)

    if len(pdg_code_dissected) > 3 and pdg_code_dissected[-4] != 0:
        # it is baryons: we need to add anti-baryons in the final pdg list
        if ipart < 0:
            print("something wrong!")
            exit(1)
        pdg_list_wbbar.append(-ipart)

    f.write("{0}\n".format(header_list[idx]))
    for iline in range(len(decay_list[idx])):
        f.write("{0}\n".format(decay_list[idx][iline]))
f.close()

f2 = open("chosen_particles_SMASH.dat", "w")
for ipart in pdg_list_wbbar:
    f2.write("{0}\n".format(ipart))
f2.close()
