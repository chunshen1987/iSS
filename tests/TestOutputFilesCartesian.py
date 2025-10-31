#!/usr/bin/env python3

import numpy as np
Nfailed = 0

data1 = np.loadtxt("checkReconstructedTmunu_test1.dat")
ResSum1 = np.mean(abs(data1[:, 2]))

data2 = np.loadtxt("checkReconstructedTmunu_test2.dat")
ResSum2 = np.mean(abs(data2[:, 2]))

data3 = np.loadtxt("checkReconstructedTmunu_test3.dat")
ResSum3 = np.mean(abs(data3[:, 2]))

data4 = np.loadtxt("checkReconstructedTmunu_test4.dat")
ResSum4 = np.mean(abs(data4[:, 2]))

data1_spectra = np.loadtxt("check_211_spectra_test1.dat")
data2_spectra = np.loadtxt("check_211_spectra_test2.dat")
data3_spectra = np.loadtxt("check_211_spectra_test3.dat")

# Check that the spectra are similar within 0.1%
tolerance = 0.001
for i in range(len(data1_spectra)):
    val1 = data1_spectra[i,1]
    val2 = data2_spectra[i,1]
    val3 = data3_spectra[i,1]
    if val1 > 1e-8 and abs(val1 - val2) / val1 > tolerance:
        print(f"Spectra test 1 and 2 differ more than {tolerance*100}% at pT={data1_spectra[i,0]}")
        Nfailed = 1
    if val1 > 1e-8 and abs(val1 - val3) / val1 > tolerance:
        print(f"Spectra test 1 and 3 differ more than {tolerance*100}% at pT={data1_spectra[i,0]}")
        Nfailed = 1

if abs(ResSum1) > 0.001:
    print("Diff 1: ", ResSum1)
    Nfailed = 1
if abs(ResSum2) > 0.001:
    print("Diff 2: ", ResSum2)
    Nfailed = 1
if abs(ResSum3) > 0.001:
    print("Diff 3: ", ResSum3)
    Nfailed = 1
if abs(ResSum4) > 0.001:
    print("Diff 4: ", ResSum4)
    Nfailed = 1

exit(Nfailed)
