# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

iSS (iSpectraSampler) is a C++ Monte Carlo sampler that turns a hydrodynamic freeze-out hypersurface (usually from MUSIC) into hadron samples via the Cooper-Frye formula, with optional delta f corrections, resonance decays and spin polarization.

## Build

Requires GSL and zlib (OpenMP optional). CMake installs into the repo root: `iSS.e` goes in the root and `libiSS.so` goes in `lib/`.

```bash
mkdir -p build && cd build && cmake .. && make -j8 && make install && cd ..
# or: ./compile.sh   (same, but deletes build/ afterwards)
```

- Compiled with `-std=c++11` and many warnings (`-Wshadow -Wold-style-cast ...`), so keep new code warning-free.
- `cmake -Dunittest=ON ..` builds the doctest unit tests (`unittest_Spectators.e`, installed to `tests/`). Without it `DOCTEST_CONFIG_DISABLE` is defined, so `TEST_CASE` blocks in sources compile away.
- The `GNUmakefile` and `src/GNUmakefile` are a legacy make path. Use CMake.

## Running and tests

```bash
./iSS.e [param_file] [work_dir] [surface_file] [key=value ...]
```

- The defaults are `iSS_parameters.dat`, `results`, and `surface.dat`.
- Any `key=value` argument overrides a parameter from the file (`ParameterReader::readFromArguments`).
- **Run from the repo root.** Table paths are relative (`iSS_tables/...`, and `eos_4D.cpp` hard-codes `./iSS_tables/EOS_tables/HRG4D/` and `./iSS_tables/deltaf_tables/urqmd/`).
- `work_dir` must contain a `music_input` file. iSS reads hydro settings from it: `EOS_to_use`, `Include_Bulk_Visc_Yes_1_No_0`, `Include_Rhob_Yes_1_No_0`, `turn_on_baryon_diffusion`, `freeze_surface_in_binary`, `output_vorticity`.

Tests are end-to-end sampling checks on a single fluid cell (the same ones CI runs in `.github/workflows/BuildTest.yml`). Run them from the repo root after `make install`:

```bash
bash tests/RunAutoTestIdeal1.sh      # one script per test: RunAutoTest{Ideal,Viscous}{1..4}.sh
```

- Each script copies a `tests/music_input_<EOS>` to `tests/music_input` and runs `iSS.e` on a `tests/test*OneFluidCell*.dat` surface. Then `tests/TestOutputFiles.py` (needs numpy) checks that the mean residual in `checkReconstructedTmunu.dat` is below 1e-3.
- The test parameter files set `perform_checks = 1`, which compares T^{μν} and charge currents rebuilt from the samples against the surface.
- All sampling tests use `MC_sampling = 4` (FSSW) and `perform_decays = 0`. CI also builds and runs `tests/TestResonanceDecayAtRest.cpp` (`g++ -std=c++11 -Isrc tests/TestResonanceDecayAtRest.cpp src/particle_decay.cpp src/Random.cpp`, run from the repo root), which checks four-momentum conservation of single decays at rest and at γ up to 300. Sampling followed by decays, and the older `EmissionFunctionArray` path, are not covered, so check changes to those by hand (e.g. `./iSS.e tests/iSS_parameters.dat tests <surface> perform_decays=1`).
- The Viscous tests take several minutes each and use about 3.5 GB of RAM.

## Architecture

The flow is `main.cpp` → `iSS::shell()` (`iSS.cpp`):

1. **Surface input**: `read_FOdata` (`readindata.cpp`) reads the surface according to `hydro_mode` (VISH2+1, or MUSIC boost-invariant or 3+1D, text or binary) and the `music_input` flags.
   - It sets `iEOS_MUSIC` and picks the afterburner, which also selects the particle table. EOS 9 forces UrQMD and 91 forces SMASH.
   - With `regulateEOS = 1`, it recomputes T and μ from an HRG EoS table, only for EoS 9, 91, 12, 14 and 20 (others are left as they are). `EOS_to_use = 20` uses the 4D (μB, μQ, μS) `EOS_4D` class instead of the 2D tables.
   - `iSS` then transforms the cells to the local rest frame (`FO_surf_LRF`).
2. **Sampling** (`MC_sampling`):
   - `4` → `FSSW` (`FSSW.cpp`): the fast, default sampler. It samples hadron counts per cell, then momenta in the LRF via `MomentumSamplerBase`, `BosonMomentumSampler` and `FermionMomentumSampler`, then boosts them back to the lab frame.
   - `1`–`3` → the older `EmissionFunctionArray` (`emissionfunction.cpp`), which builds spectra and samples them with the `pT/phi/eta` Gauss tables in `iSS_tables/bin_tables`.
3. **Delta f**: `bulk_deltaf_kind` selects the model.
   - `0`/`1`/`11` are 14-moment, RTA and OSU.
   - `20`/`21` are the WSU 22-moment and Chapman-Enskog coefficients for NEoS-BQS. These come from `iSS_tables/deltaf_tables/<afterburner>/`, or from the 4D binary coefficient files via `EOS_4D` when the 4D EoS is used.
4. **Decays**: with `perform_decays = 1`, `particle_decay` uses `pdg-urqmd_v3.3+.dat`, `pdg-s95pv1.dat` or `pdg-SMASH.dat`. Its callers, the feed-down loops in `FSSW.cpp` and `emissionfunction.cpp`, re-queue every daughter that `check_particle_stable` reports as unstable. FSSW decays events in parallel, one random stream per event, through the `perform_decays` overload that takes the generator. FSSW skips decays for SMASH.
5. **Other options**: `calculate_polarization` → `SpinPolarization` (needs vorticity arrays in the surface); `include_spectators` → `Spectators`; `local_charge_conservation`.

`iSS` can also be used as a library (`libiSS.so`) by frameworks such as JETSCAPE/XSCAPE. They call `shell()` and read results through `get_number_of_sampled_events`, `get_hadron` and the related getters.

## Data tables

- `iSS_tables/` holds the particle lists (`pdg-*.dat`, `chosen_particles_*.dat`), the HRG EoS tables and the delta f coefficients.
- The 4D EoS and its delta f coefficients are not in git. Download them with:
  - `iSS_tables/EOS_tables/download_HRG4D.sh` (creates `HRG4D/`), run from `iSS_tables/EOS_tables/`;
  - `iSS_tables/deltaf_tables/urqmd/download_NEoS4D_deltafCoeffs.sh`, run from that same `urqmd/` folder.

## Notes

`doc/` holds development reports: merge history, review fixes, 4D EoS validation (including an open mismatch between the HRG4D table and the sampled gas) and the MUSIC surface column layout (`doc/surface_text_format.md`).

## Branches

- `main` is the integration branch.
- `XSCAPE` and `JETSCAPE` hold the framework-integration variants. PRs from outside contributors often target `XSCAPE`.
