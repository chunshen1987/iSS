#!/usr/bin/env python3
"""Write the spacelike, flowing single-cell test of doc/FSSW_performance.md.

It starts from tests/testViscousOneFluidCell2.dat and sets
dsigma_mu = tau (1000, 3000, 0, 0) (lab |dsigma^i|/dsigma^0 = 3) and
u = (0.5, 0.3, 0). It writes <work_dir>/spacelike.dat (text) and a
matching <work_dir>/music_input (EoS 9, bulk on, rho_B off).

Usage: make_spacelike_cell.py WORK_DIR   (run from the iSS repo root)

Sample it with, for example,
  ./iSS.e iSS_parameters.dat WORK_DIR spacelike.dat \
      number_of_repeated_sampling=1500 use_binary_format=1
and compare two codes or seeds with sample_stats.py. The sampler's flux
acceptance matters most here, because half of the directions have
negative flux.
"""
import argparse
import os

MUSIC_INPUT = """EOS_to_use  9
Include_Bulk_Visc_Yes_1_No_0  1
Include_Rhob_Yes_1_No_0  0
turn_on_baryon_diffusion  0
freeze_surface_in_binary  0
EndOfData
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work_dir")
    parser.add_argument("--template",
                        default="tests/testViscousOneFluidCell2.dat")
    args = parser.parse_args()

    with open(args.template) as f:
        cell = [float(x) for x in f.read().split()]
    cell[4:8] = [1000., 3000., 0., 0.]    # dsigma_mu without the factor tau
    cell[8:12] = [1.1576, 0.5, 0.3, 0.]   # u^mu (iSS recomputes u^0)

    os.makedirs(args.work_dir, exist_ok=True)
    with open(os.path.join(args.work_dir, "spacelike.dat"), "w") as f:
        f.write("  ".join(f"{x:.10g}" for x in cell) + "\n")
    with open(os.path.join(args.work_dir, "music_input"), "w") as f:
        f.write(MUSIC_INPUT)


if __name__ == "__main__":
    main()
