#!/usr/bin/env python

sign = lambda a: (a>0) - (a<0)

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
    temp = int(ipart)
    pdg_code_dissected = [int(d) for d in str(abs(temp))]
    if len(pdg_code_dissected) < 3:
        # photon
        continue
    if pdg_code_dissected[-2] > 3 or pdg_code_dissected[-3] > 3:
        # heavy quark mesons or baryons
        continue
    if len(pdg_code_dissected) > 3 and pdg_code_dissected[-4] > 3:
        # heavy quark baryons
        continue
    if len(pdg_code_dissected) > 3 and pdg_code_dissected[-4] != 0:
        # baryons
        pdg_list_final.append(-temp)
    pdg_list_final.append(temp)

for ipart in pdg_list_final:
    print(ipart)
