#!/usr/bin/env python3
"""Write a synthetic constant energy density MUSIC 3+1D binary surface.

This is the 300K-cell surface of doc/FSSW_performance.md. It writes
<work_dir>/surface.dat (36 floats per cell, see doc/surface_text_format.md)
and a matching <work_dir>/music_input (EoS 9, bulk on, rho_B off).

The surface is tau = tau_f(r) on an (x, y, eta_s) grid, with
e = 0.26 GeV/fm^3 everywhere, radial and elliptic flow, shear stress and a
negative bulk pressure. Cells with |grad tau_f| > 1 are spacelike.

Usage: make_benchmark_surface.py WORK_DIR [--seed N]
"""
import argparse
import os

import numpy as np

HBARC = 0.197327053  # GeV fm

E_DEC = 0.26         # GeV/fm^3
T_GUESS = 0.150      # GeV, recomputed by iSS with regulateEOS = 1
P_GUESS = 0.060      # GeV/fm^3, recomputed by iSS with regulateEOS = 1
RADIUS = 8.0         # fm
DX = 0.25            # fm
DETA = 0.1
ETA_MAX = 4.6
TAU_CENTER, TAU_EDGE = 9.0, 2.0  # fm
MAX_GRADIENT = 2.0   # cap on |d tau_f / dr|

MUSIC_INPUT = """EOS_to_use  9
Include_Bulk_Visc_Yes_1_No_0  1
Include_Rhob_Yes_1_No_0  0
turn_on_baryon_diffusion  0
freeze_surface_in_binary  1
output_vorticity  0
EndOfData
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work_dir")
    parser.add_argument("--seed", type=int, default=12345,
                        help="seed of the per-cell fluctuations")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)

    xs = np.arange(-RADIUS + DX/2, RADIUS, DX)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    r = np.hypot(X, Y)
    inside = r < RADIUS
    X, Y, r = X[inside], Y[inside], r[inside]
    phi = np.arctan2(Y, X)
    etas = np.arange(-ETA_MAX, ETA_MAX + 1e-9, DETA)

    # tau_f(r) = tau_edge + (tau_c - tau_edge) sqrt(1 - (r/R)^2)
    s = np.sqrt(np.clip(1 - (r/RADIUS)**2, 1e-4, None))
    tau_f = TAU_EDGE + (TAU_CENTER - TAU_EDGE)*s
    dF_dr = np.clip(-(TAU_CENTER - TAU_EDGE)*(r/RADIUS**2)/s,
                    -MAX_GRADIENT, 0.)
    dF_dx = dF_dr*X/np.maximum(r, 1e-12)
    dF_dy = dF_dr*Y/np.maximum(r, 1e-12)

    n = X.size
    volume = DX*DX*DETA
    slices = []
    for eta in etas:
        a = np.zeros((n, 36))
        a[:, 0], a[:, 1], a[:, 2], a[:, 3] = tau_f, X, Y, eta
        # MUSIC stores dsigma_mu without the factor tau in components 0-2
        a[:, 4] = volume
        a[:, 5] = -dF_dx*volume
        a[:, 6] = -dF_dy*volume
        # flow: radial and elliptic, Bjorken in eta_s (u^eta = 0)
        ur = (0.9*(r/RADIUS)*(1 + 0.1*np.cos(2*phi))
              *np.exp(-0.5*(eta/3.0)**2))
        ux = ur*np.cos(phi) + 0.02*rng.standard_normal(n)
        uy = ur*np.sin(phi) + 0.02*rng.standard_normal(n)
        a[:, 8] = np.sqrt(1 + ux**2 + uy**2)
        a[:, 9], a[:, 10] = ux, uy
        a[:, 12] = E_DEC/HBARC                              # fm^-4
        a[:, 13] = T_GUESS/HBARC                            # fm^-1
        a[:, 17] = (E_DEC + P_GUESS)/T_GUESS/HBARC**3       # fm^-3
        # shear: iSS (regulate_Wmunu) makes it transverse and traceless
        A = 0.08*(E_DEC + P_GUESS)*(1 + 0.3*rng.standard_normal(n))
        pixx = A*(1 + 0.2*np.cos(2*phi)) + 0.01*rng.standard_normal(n)
        piyy = A*(1 - 0.2*np.cos(2*phi)) + 0.01*rng.standard_normal(n)
        a[:, 22] = pixx/HBARC
        a[:, 23] = 0.01*rng.standard_normal(n)/HBARC
        a[:, 25] = piyy/HBARC
        a[:, 27] = -(pixx + piyy)/HBARC
        # bulk pressure
        a[:, 28] = -0.015*(1 + 0.3*rng.standard_normal(n))/HBARC
        slices.append(a.astype(np.float32))
    surface = np.concatenate(slices)

    os.makedirs(args.work_dir, exist_ok=True)
    surface.tofile(os.path.join(args.work_dir, "surface.dat"))
    with open(os.path.join(args.work_dir, "music_input"), "w") as f:
        f.write(MUSIC_INPUT)

    tau = surface[:, 0]
    d0, d1, d2 = tau*surface[:, 4], tau*surface[:, 5], tau*surface[:, 6]
    flux = surface[:, 8]*d0 + surface[:, 9]*d1 + surface[:, 10]*d2
    spacelike = d0**2 - d1**2 - d2**2 < 0
    print(f"{surface.shape[0]} cells ({n} per slice, {etas.size} slices), "
          f"sum u.dsigma = {flux.sum():.1f} fm^3")
    print(f"spacelike: {100*spacelike.mean():.1f}% of cells, "
          f"{100*flux[spacelike].sum()/flux.sum():.1f}% of u.dsigma")


if __name__ == "__main__":
    main()
