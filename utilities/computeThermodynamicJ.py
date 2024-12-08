#!/urs/bin/env python3

"""
    This script compute thermodynamic integrals with special functions
"""

import mpmath
import numpy as np

mpmath.mp.dpf = 50
mpmath.mp.pretty = True


def computeThermodynamicJ(nu: float, m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_nu(m)
        Eq. (9) in https://www.overleaf.com/project/674f23fdff4ede2888cab91e
    """
    res = 0.
    if m < 20.:
        truncOrder = 5
    else:
        truncOrder = 1
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder + 1):
        arg = (k * m)**2 / 4.
        res += sign**(k-1) * k * (
            1./(k**(nu + 2)) * mpmath.gamma(nu + 2)
            * mpmath.hyp1f2(-0.5, -0.5 - nu/2., -nu/2., arg)
            + m**(nu + 2)/4.*np.sqrt(np.pi)*(
                mpmath.gamma(-1 - nu/2.) / mpmath.gamma(0.5 - nu/2.)
                * mpmath.hyp1f2(0.5 + nu/2., 0.5, 2 + nu/2., arg)
                - k * m * mpmath.gamma(-1.5 - nu/2.) / mpmath.gamma(- nu/2.)
                * mpmath.hyp1f2(1 + nu/2., 1.5, 2.5 + nu/2., arg))
        )
    return res


def thermodynamicJ0(m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_0(m)
        Eq. (10) in https://www.overleaf.com/project/674f23fdff4ede2888cab91e
    """
    res = 0.
    truncOrder = 5
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder):
        arg = (k * m)
        res += sign**(k-1) * m * mpmath.besselk(1, arg)
    return res


def thermodynamicJ1(m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_1(m)
        Eq. (11) in https://www.overleaf.com/project/674f23fdff4ede2888cab91e
    """
    res = 0.
    truncOrder = 5
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder):
        arg = (k * m)
        res += sign**(k-1) * m**2 * mpmath.besselk(2, arg)
    return res


def thermodynamicJ2(m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_2(m)
        Eq. (12) in https://www.overleaf.com/project/674f23fdff4ede2888cab91e
    """
    res = 0.
    truncOrder = 5
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder):
        arg = (k * m)
        res += sign**(k-1) * m**2 / k * (
            arg * mpmath.besselk(1, arg) + 3. * mpmath.besselk(2, arg)
        )
    return res


def thermodynamicJ3(m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_3(m)
        Eq. (13) in https://www.overleaf.com/project/674f23fdff4ede2888cab91e
    """
    res = 0.
    truncOrder = 5
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder):
        arg = (k * m)
        res += sign**(k-1) * m**3 / k * (
            arg * mpmath.besselk(2, arg) + 3. * mpmath.besselk(3, arg)
        )
    return res


def thermodynamicJ0p5(m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_0p5(m)
    """
    res = 0.
    if m < 5.:
        truncOrder = 5
    else:
        truncOrder = 2
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder):
        arg = (k * m) / 2.
        res += sign**(k-1) * k * np.sqrt(np.pi) / (16. * k**2.5) * (
                np.sqrt(2.) * k**2 * m**2 * np.pi * mpmath.besseli(-0.25, arg)**2
                + np.sqrt(2.) * np.pi * (
                    (-3 + k**2 * m**2) * mpmath.besseli(0.25, arg)**2
                    + 3 * k * m * mpmath.besseli(0.25, arg) * mpmath.besseli(0.75, arg)
                    + 3 * k * m * mpmath.besseli(-0.75, arg) * (
                            mpmath.besseli(0.25, arg)
                            - 2 * k * m * mpmath.besseli(0.75, arg))
                    + 3 * k**2 * m**2 * (
                            mpmath.besseli(0.75, arg)**2
                            + mpmath.besseli(1.25, arg)**2)
                )
                - 2 * k * m * mpmath.besseli(-0.25, arg) * (
                    np.sqrt(2) * k * m * np.pi * mpmath.besseli(0.25, arg)
                    - 3 * mpmath.besselk(-0.75, arg)
                )
        )
    return res


def thermodynamicJ1p5(m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_1p5(m)
    """
    res = 0.
    if m < 5.:
        truncOrder = 5
    else:
        truncOrder = 2
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder):
        arg = (k * m) / 2.
        res += sign**(k-1) * k * np.sqrt(np.pi) / (32. * k**3.5) * (
                5 * np.sqrt(2.) * k**2 * m**2 * np.pi * mpmath.besseli(-0.25, arg)**2
                + np.sqrt(2.) * np.pi * (
                    5 * (-3 + k**2 * m**2) * mpmath.besseli(0.25, arg)**2
                    + k * m * (15 + 8 * k**2 * m**2) * mpmath.besseli(0.25, arg)
                            * mpmath.besseli(0.75, arg)
                    + k * m * mpmath.besseli(-0.75, arg) * (
                            (15 - 8 * k**2 * m**2) * mpmath.besseli(0.25, arg)
                            - 30 * k * m * mpmath.besseli(0.75, arg))
                    + 15 * k**2 * m**2 * (
                            mpmath.besseli(0.75, arg)**2
                            + mpmath.besseli(1.25, arg)**2)
                )
                + 2 * k * m * mpmath.besseli(-0.25, arg) * (
                    - 5 * np.sqrt(2) * k * m * np.pi * mpmath.besseli(0.25, arg)
                    + (15 + 8 * k**2 * m**2) * mpmath.besselk(-0.75, arg)
                )
        )
    return res


def thermodynamicJ2p5(m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_2p5(m)
    """
    res = 0.
    if m < 5.:
        truncOrder = 5
    else:
        truncOrder = 2
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder):
        arg = (k * m) / 2.
        res += sign**(k-1) * k * np.sqrt(np.pi) / (64 * k**4.5) * (
                np.sqrt(2.) * k**2 * m**2 * np.pi * (35 + 8 * k**2 * m**2)
                    * mpmath.besseli(-0.25, arg)**2
                + np.sqrt(2.) * np.pi * (
                    (-105 + 27 * k**2 * m**2 + 8 * k**4 * m**4)
                        * mpmath.besseli(0.25, arg)**2
                    + k * m * (105 + 64 * k**2 * m**2) * mpmath.besseli(0.25, arg)
                            * mpmath.besseli(0.75, arg)
                    - k * m * mpmath.besseli(-0.75, arg) * (
                            3 * (-35 + 16 * k**2 * m**2) * mpmath.besseli(0.25, arg)
                            + 2 * k * m * (105 + 8 * k**2 * m**2) * mpmath.besseli(0.75, arg))
                    + k**2 * m**2 * (105 + 8 * k**2 * m**2) * (
                            mpmath.besseli(0.75, arg)**2
                            + mpmath.besseli(1.25, arg)**2)
                )
                - 2 * k * m * mpmath.besseli(-0.25, arg) * (
                    np.sqrt(2) * k * m * (35 + 8 * k**2 * m**2) * np.pi * mpmath.besseli(0.25, arg)
                    - (105 + 64 * k**2 * m**2) * mpmath.besselk(-0.75, arg)
                )
        )
    return res


def thermodynamicJ3p5(m: float, bosonFlag: bool) -> float:
    """
        Compute thermodynamic integral J_3p5(m)
    """
    res = 0.
    if m < 5.:
        truncOrder = 5
    else:
        truncOrder = 2
    sign = 1.
    if bosonFlag:
        sign = -1.
    for k in range(1, truncOrder):
        arg = (k * m) / 2.
        res += sign**(k-1) * k * np.sqrt(np.pi) / (128 * k**5.5) * (
                5 * np.sqrt(2.) * k**2 * m**2 * np.pi * (63 + 16 * k**2 * m**2)
                    * mpmath.besseli(-0.25, arg)**2
                + np.sqrt(2.) * np.pi * (
                    (-945 + 219 * k**2 * m**2 + 80 * k**4 * m**4)
                        * mpmath.besseli(0.25, arg)**2
                    + k * m * (945 + 600 * k**2 * m**2 + 32 * k**4 * m**4)
                            * mpmath.besseli(0.25, arg) * mpmath.besseli(0.75, arg)
                    - k * m * mpmath.besseli(-0.75, arg) * (
                            (-945 + 408 * k**2 * m**2 + 32 * k**4 * m**4) * mpmath.besseli(0.25, arg)
                            + 6 * k * m * (315 + 32 * k**2 * m**2) * mpmath.besseli(0.75, arg))
                    + 3 * k**2 * m**2 * (315 + 32 * k**2 * m**2) * (
                            mpmath.besseli(0.75, arg)**2
                            + mpmath.besseli(1.25, arg)**2)
                )
                - 2 * k * m * mpmath.besseli(-0.25, arg) * (
                    5 * np.sqrt(2) * k * m * (63 + 16 * k**2 * m**2) * np.pi * mpmath.besseli(0.25, arg)
                    - (945 + 600 * k**2 * m**2 + 32 * k**4 * m**4) * mpmath.besselk(-0.75, arg)
                )
        )
    return res

#print(computeThermodynamicJ(0.01, 0.5, False), thermodynamicJ0(0.5, False))
#print(computeThermodynamicJ(0.01, 0.5, True), thermodynamicJ0(0.5, True))
#print(computeThermodynamicJ(1.01, 0.5, False), thermodynamicJ1(0.5, False))
#print(computeThermodynamicJ(1.01, 0.5, True), thermodynamicJ1(0.5, True))
#print(computeThermodynamicJ(2.01, 0.5, False), thermodynamicJ2(0.5, False))
#print(computeThermodynamicJ(2.01, 0.5, True), thermodynamicJ2(0.5, True))
#print(computeThermodynamicJ(3.01, 0.5, False), thermodynamicJ3(0.5, False))
#print(computeThermodynamicJ(3.01, 0.5, True), thermodynamicJ3(0.5, True))

m_min = 0.5
m_max = 15.
dm = 0.05
n_points = int((m_max - m_min) / dm) + 1
mList = np.linspace(m_min, m_max, n_points)

#for m_i in mList:
#    print(f"{m_i:.6e}, {float(thermodynamicJ3p5(m_i, False)):.6e}, {float(computeThermodynamicJ(3.5, m_i, False)):.6e}")
#exit(0)

gammaList = [0.0, 0.5, 1.0, 1.5]

for bosonFlag in [True, False]:
    if bosonFlag:
        f1 = open("J_gamma_BE.dat", "w")
        f2 = open("J_gammaplus2_BE.dat", "w")
        f3 = open("J_2_BE.dat", "w")
    else:
        f1 = open("J_gamma_FD.dat", "w")
        f2 = open("J_gammaplus2_FD.dat", "w")
        f3 = open("J_2_FD.dat", "w")
    f1.write(f"# m/T  J_gamma  (gamma = {gammaList})\n")
    f2.write(f"# m/T  J_gamma+2  (gamma = {gammaList})\n")
    f3.write(f"# m/T  J_2\n")
    for m_i in mList:
        J_2 = thermodynamicJ2(m_i, bosonFlag)
        f3.write(f"{m_i:.8e}  {float(J_2):.8e}\n")
        f1.write(f"{m_i:.8e}")
        f2.write(f"{m_i:.8e}")
        for gamma_i in gammaList:
            if abs(gamma_i - 0) < 1e-10:
                J_gamma = thermodynamicJ0(m_i, bosonFlag)
                J_gammaplus2 = thermodynamicJ2(m_i, bosonFlag)
            elif abs(gamma_i - 0.5) < 1e-10:
                J_gamma = thermodynamicJ0p5(m_i, bosonFlag)
                J_gammaplus2 = thermodynamicJ2p5(m_i, bosonFlag)
            elif abs(gamma_i - 1) < 1e-10:
                J_gamma = thermodynamicJ1(m_i, bosonFlag)
                J_gammaplus2 = thermodynamicJ3(m_i, bosonFlag)
            elif abs(gamma_i - 1.5) < 1e-10:
                J_gamma = thermodynamicJ1p5(m_i, bosonFlag)
                J_gammaplus2 = thermodynamicJ3p5(m_i, bosonFlag)
            else:
                J_gamma = computeThermodynamicJ(gamma_i, m_i, bosonFlag)
                J_gammaplus2 = computeThermodynamicJ(gamma_i + 2, m_i, bosonFlag)
            f1.write(f"  {float(J_gamma):.8e}")
            f2.write(f"  {float(J_gammaplus2):.8e}")
        f1.write("\n")
        f2.write("\n")
    f1.close()
    f2.close()
    f3.close()
