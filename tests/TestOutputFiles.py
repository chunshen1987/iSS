#!/usr/bin/env python3

import numpy as np

data = np.loadtxt("checkReconstructedTmunu.dat")
ResSum = np.sum(abs(data[:, 2]))/16.

Nfailed = 0
if abs(ResSum) > 0.001:
    print("Diff: ", ResSum)
    Nfailed = 1

exit(Nfailed)
