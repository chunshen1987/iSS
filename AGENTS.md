# Scientific programming in iSS

Work as a scientific programmer specializing in theoretical relativistic
heavy-ion physics. iSS (iSpectraSampler) turns a hydrodynamic freeze-out
hypersurface into hadron spectra and Monte Carlo events using the Cooper–Frye
prescription. Use the equations, physical limits, and numerical evidence to
evaluate changes alongside software correctness.

## Physics and numerical methods

- For a change to a physics kernel, identify the governing equation, its
  approximations, and how its symbols map to code. Check the original model
  or table documentation when a coefficient or convention is unclear. Give
  the reference and equation number when available.
- Trace units, signs, tensor indices, and frames across the input, calculation,
  and output. In particular, check the metric and index placement in
  `p^mu dSigma_mu`, the Milne and Cartesian representations, and the boosts
  between the fluid local rest frame and the lab. Use `hbarC` from
  `src/data_struct.h` for the existing GeV and fm conversions.
- Check how the Cooper–Frye flux, equilibrium distribution, viscous and
  diffusion corrections, and regulators enter both the integrated yield and
  the sampled momentum distribution. Preserve each mode's treatment of
  signed flux, quantum statistics, and charge conjugation unless the change
  calls for a new prescription. Check shared physics across sampler paths
  when applicable.
- For equation-of-state or coefficient-table changes, verify the axis and
  field ordering, units, interpolation, and behavior at table boundaries.
  Check that the hadron species, degeneracies, and baryon, electric, and
  strangeness charges match the selected table and afterburner.
- Prefer an analytic limit or an independent numerical calculation as a
  reference. Useful limits include an ideal static cell, vanishing viscous
  corrections, and charge conjugation at finite chemical potentials. Check
  mass shells and four-momentum conservation for decay changes. For sampled
  events, test ensemble means and uncertainties unless the selected mode
  imposes an event-by-event constraint.
- Justify numerical tolerances with the scale of the quantity, convergence,
  or Monte Carlo uncertainty. When changing an approximation, cutoff, or
  fallback, explain the affected regime and quantify its effect on relevant
  observables with a reproducible case.

## Code paths and compatibility

- `src/iSS.cpp` orchestrates surface loading, frame transformations, and
  dispatch. `src/readindata.cpp` implements hydro input and EoS regulation.
  Consult [the MUSIC surface layout](doc/surface_text_format.md) for text and
  binary formats, and [the 4D EoS validation notes](doc/EoS4D_validation.md)
  for the reported mismatch between the HRG4D table and the sampled gas.
- `src/FSSW.cpp` and the momentum sampler classes implement the default fast
  sampler (`MC_sampling = 4`). `src/emissionfunction.cpp` implements smooth
  spectra and sampling modes 1–3. `src/eos_4D.cpp` handles EoS 20 and its 4D
  viscous coefficient tables. Resonance feed-down uses `src/particle_decay.cpp`
  from both sampler paths. The [README](README.md) describes the decay model.
- Preserve documented CLI and parameter behavior, MUSIC and afterburner file
  formats, and the public `iSS` library interface used by downstream
  frameworks. Check both CLI and library behavior when changing initialization,
  event storage, or repeated calls to `shell()`.
- Keep production code C++11-compatible and use the existing RNG and OpenMP
  patterns. When modifying sampling or threading, track random draw order
  and ownership. Preserve existing fixed-seed guarantees. For intentional RNG
  changes, validate the resulting distributions and record the change in
  reproducibility.

## Build and validate

Use CMake with GSL and zlib. OpenMP is optional. From the repository root,
build and run a focused check without installing over local binaries.

```sh
cmake -S . -B build
cmake --build build --parallel 2
python3 tests/TestSamplerConsistency.py build/src/iSS.e
```

- For normal runs, invoke `./iSS.e [param_file] [work_dir] [surface_file]
  [key=value ...]` from the repository root so relative `iSS_tables/` paths
  resolve. For MUSIC input, put `music_input` and the surface in `work_dir`.
  Use `iSS_parameters.dat` to check defaults and command-line override names.
- Choose checks for the affected feature. `tests/` includes focused checks
  for MUSIC parsers, 4D tables, charge conjugation, spectra, and oversampling.
  `tests/Test4DEOSTables.py` uses synthetic tables. Physical EoS 20 runs
  require the external tables listed in the 4D EoS validation notes.
- CI runs the single-cell ideal and viscous scripts and the standalone
  `tests/TestResonanceDecayAtRest.cpp` test. The single-cell scripts expect
  `./iSS.e` after installation, require NumPy, and overwrite the tracked
  `tests/music_input`. Run them in a disposable checkout or isolate their
  inputs and outputs. These scripts exercise FSSW reconstruction with
  decays disabled. Check other sampling modes and sampling with feed-down
  directly when they are affected.
- For numerical or stochastic comparisons, record the input and table
  choices, parameter overrides, `randomSeed`, event count, and OpenMP thread
  count as applicable. Report the observable, reference, and residual or
  uncertainty. For performance work, compare wall time and physics observables
  on the same workload using [the benchmark utilities](utilities/benchmark/README.md).
