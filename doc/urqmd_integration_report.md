# iSS–UrQMD integration report

## Summary

iSS now supports two paths from Cooper–Frye samples to UrQMD. Both work with `MC_sampling = 1` through `4` and use the UrQMD particle table selected by `afterburner_type = 1`.

| Mode | Setting | Event handoff |
| --- | --- | --- |
| Direct restart file | `use_UrQMD_format = 1` | iSS writes `fort.14`. Pass it to standalone UrQMD as `ftn10`. |
| In-memory transport | `run_urqmd_in_memory = 1` | iSS passes each sampled event to the linked UrQMD 3.4 cascade and exposes the transported particles through its event getters. |

## Direct restart file

The C++ writer follows the active `osc2u` conversion path. It maps PDG IDs to UrQMD types, orders baryons before mesons, and back-propagates hadrons to time zero while preserving their original hadronization times. Empty samples are omitted. The third line of `iSS_tables/OSCAR_header.txt` supplies the projectile, target, reference frame, and lab kinetic energy per nucleon. The shipped metadata describe Au+Au at 100 GeV per nucleon in the equal-speed frame.

## In-memory transport

Build with GNU Fortran and `patch` available. CMake uses `urqmd_code/urqmd/src` by default. Set `ISS_URQMD_SOURCE_DIR` to use a different local copy of the UrQMD 3.4 sources.

```sh
cmake -S . -B build -DISS_ENABLE_URQMD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Set `run_urqmd_in_memory = 1` and `afterburner_type = 1` in `iSS_parameters.dat`. This takes precedence over restart-file and OSCAR output. iSS holds the sampled events in memory, transports them one at a time, and replaces each event with its final particles. Enabled sampler checks run before transport.

The linked cascade uses UrQMD's default physics options with `ctp 22 1.0` from the supplied burner file. The `urqmd_duration` and `urqmd_max_step` settings are in fm/c and both default to `8000`. By default, the first transport event uses iSS's resolved `randomSeed`. With `randomSeed = -1`, both modules share one seed from the local system clock. A positive `urqmd_seed` overrides the transport seed.

The event handoff skips OSCAR, `fort.14`, and per-species sample files. The existing gzip and binary options can save final particles. Set both to `0` to keep final events in memory. UrQMD may still read or create its initialization cache, `tables.dat`.

C++ callers can link the UrQMD-enabled `libiSS` and use `UrqmdTransport::evolve` with a `std::vector<iSS_Hadron>`. The result contains the final particles and the number of collision attempts. The wrapper serializes access to UrQMD's shared Fortran state. [C++ interface](../src/urqmd_transport.h)

## Verification

- Fixed-seed restart files matched the compiled `osc2u` byte for byte in sampling modes 2 and 4. Modes 1 and 3 matched within iSS's in-memory particle precision. Standalone UrQMD successfully processed direct restart output.
- The in-memory C++ tests exercised charged-rho decay, proton–proton scattering, empty input, event isolation, and seeded replay. [C++ test](../tests/TestUrQMDTransport.cpp)
- Command-line tests passed for all four sampling modes and checked final binary and gzip output, memory-only operation, and the absence of intermediate event files. The build with the Fortran interface disabled also passed the existing restart-output and OSCAR2013 regressions. [CLI test](../tests/TestUrQMDTransport.py)
