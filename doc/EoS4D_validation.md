# Validation of the 4D EoS path (EOS_to_use = 20)

The CI tests do not cover EoS 20, so these runs were done by hand with
the downloaded tables:

- `iSS_tables/EOS_tables/HRG4D/HRG_{p,t,mub,muq,mus}_b.dat`, from
  `download_HRG4D.sh` (UrQMD)
- `iSS_tables/deltaf_tables/urqmd/NEoS4D_{22mom,CE}_deltafCoeff.bin`,
  from `download_NEoS4D_deltafCoeffs.sh`

## Table edges (review finding #1)

A small program called `EOS_4D::getThermalVariables` at points inside
and beyond the table, built once with `eos_4D.cpp` before the fix and
once after.

| Input (e in GeV/fm³; nB, nQ, nS in 1/fm³) | Before the fix | After the fix |
|---|---|---|
| e = 0.3, all n = 0 (inside) | T = 0.1530, P = 0.04976 | identical |
| e = 1000 (above the T range) | throws `out_of_range` | T = 0.170 GeV (table top) |
| e = 0.3, nB = 50 | **μB = −0.70** (read from the neighbouring row) | μB = +0.70 (table edge) |
| e = 0.3, nQ = 50 | μQ = −0.13 | μQ = +0.13 |
| e = 0.3, nS = 50 | μB = −0.70 | μB = +0.70 |
| e = 1000, all n = 500 | throws `out_of_range` | table-top values |

## End-to-end sampling

All runs use the BQS viscous CI cell `tests/testViscousOneFluidCell3.dat`
with `tests/music_input_14` changed to `EOS_to_use 20` and the CI
parameter files (`regulateEOS = 1`). "Residual" is the mean
|ΔT^{μν}| that `tests/TestOutputFiles.py` checks against 1e-3.

| Case | Settings | Result |
|---|---|---|
| A | 22-moment δf (kind 20) | runs; residual 1.04e-3 |
| B | Chapman-Enskog δf (kind 21) | runs; residual 5.2e-4 |
| C | RTA δf (kind 1), 22-moment file removed | runs (exited before fix #4); residual 9.2e-4 |
| D | `regulateEOS = 0`, kind 20 | stops with the error from fix #3 |
| E | no δf | runs; residual 7.8e-4 |
| F | kind 20, cell with Qn = 0.4 nB and Sn = 0 added after Bn | runs; residual 1.73e-3 |

## Open issue: the HRG4D table does not match the sampled hadron gas

With EoS 20, the sampled hadrons do not reproduce the cell's own energy
density and charges. EoS 14 reproduces the same cell to about 0.1%.

| Run | T⁰⁰ sampled / surface | nB sampled / surface | nQ sampled / surface |
|---|---|---|---|
| EoS 14 (CI test Viscous3) | 0.3009 / 0.3005 (+0.1%) | 5012 / 5000 (+0.2%) | — |
| EoS 20, case E (no δf), nQ = 0 | 0.2897 / 0.3005 (−3.6%) | 4914 / 5000 (−1.7%) | — |
| EoS 20, case F, nQ = 0.4 nB | 0.2828 / 0.3005 (−5.9%) | 4796 / 5000 (−4.1%) | 1901 / 2000 (−4.9%) |

Turning δf off does not change the gap, and the edge fix leaves values
inside the table unchanged. So the gap comes from the table. The T and
μ it returns for this (e, nB, nQ, nS) do not reproduce e and n when the
UrQMD particle list (`pdg-urqmd_v3.3+.dat`) is sampled.

Possible causes:

- the table was generated from a particle list different from
  `pdg-urqmd_v3.3+.dat`;
- the grid in the tilde variables is too coarse around this point;
- the e, n → tilde-variable map in `EOS_4D::get_tilde_variables` does
  not match the one used to build the table.

This should be checked against how the HRG4D tables were made, and on
a real MUSIC EoS 20 binary surface, before EoS 20 is used for
production runs.
