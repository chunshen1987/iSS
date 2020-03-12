#!/usr/bin/env python

filename = "pdg-SMASH.dat"
pdg_list = []
with open(filename) as fp:
    line = fp.readline()
    while line:
        pdg_list.append(line.split()[0])
        line = fp.readline()
pdg_list = list(set(pdg_list))

pdg_list_final = []
for ipart in pdg_list:
    ipart = int(ipart)
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
    pdg_list_final.append(ipart)

    if len(pdg_code_dissected) > 3 and pdg_code_dissected[-4] != 0:
        # it is baryons: we need to add anti-baryons in the final pdg list
        if ipart < 0:
            print("something wrong!")
            exit(1)
        pdg_list_final.append(-ipart)


for ipart in pdg_list_final:
    print(ipart)
