#!/usr/bin/env python3

import numpy as np

data = np.loadtxt("checkReconstructedTmunu.dat")
ResSum = np.sum(abs(data[:, 1]))

Nfailed = 0
if abs(ResSum) > 0.002:
    print("Diff: ", ResSum)
    Nfailed = 1

exit(Nfailed)
