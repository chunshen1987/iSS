iSS
===

Monte Carlo sampler for particle distribution from Cooper-Frye freeze-out procedure

Resonance decays (`perform_decays = 1`)
---------------------------------------

A decay channel is drawn with the table's branching ratios, normalised over
the channels the mother can reach at its pole mass. A channel that is closed
at the daughters' pole masses is reached by giving its broad daughters (width
above 0.1 MeV) masses from their Breit-Wigner distributions, each no lower
than the mass at which all of its own channels stay open, so that the
daughter keeps the table's branching ratios; a channel that cannot be reached
that way (e.g. a0(980), f0(980) -> K Kbar in `pdg-s95pv1.dat`) is skipped and
its share goes to the others. A particle that is below all of its channels
is kept in the final state, and the number of such particles is reported at
the end of the run; with the shipped tables this does not happen, because a
daughter is never given a mass at which it cannot decay. A negative number
of daughters (how
`pdg-urqmd_v3.3+.dat` and `pdg-s95pv1.dat` mark closed channels) is read as
its absolute value. Two- and three-body channels use the original routines and
random-number sequence; four- and five-body ones are uniform in phase space
(Raubold-Lynch).

Direct UrQMD output
-------------------

Set `use_UrQMD_format = 1` to write sampled events in UrQMD's restart
format to `fort.14` in the working directory. Pass that file to UrQMD as
`ftn10` in place of the file normally produced by `osc2u`. The option takes
priority over OSCAR, gzip, and binary event output and works with
`MC_sampling = 1` through `4`. It holds sampled events in memory while
writing and omits empty samples.

Use `afterburner_type = 1` to select the UrQMD particle table. The writer
reports an error if a sampled PDG ID has no UrQMD mapping. Update the
third line of `iSS_tables/OSCAR_header.txt` for the projectile, target,
reference frame, and lab kinetic energy per nucleon. The shipped header
describes Au+Au at 100 GeV per nucleon in the equal-speed frame. As in
`osc2u`, the impact parameter is zero and the frame setting describes the
sampled momenta. The writer back-propagates the hadrons to time zero and
records their original hadronization times for UrQMD.

In-memory UrQMD transport
------------------------

Build the optional UrQMD 3.4 interface with GNU Fortran and `patch`.
CMake looks for local sources in `urqmd_code/urqmd/src` by default.
Set `ISS_URQMD_SOURCE_DIR` at configure time to use a different local copy.
If the sources are missing, CMake uses Git to clone
`Chunshen1987/urqmd_afterburner` from Bitbucket into
`build/_deps/urqmd_afterburner` and reuses that checkout on later runs.

```sh
cmake -S . -B build -DISS_ENABLE_URQMD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Set `run_urqmd_in_memory = 1` and `afterburner_type = 1` in
`iSS_parameters.dat`. This works with `MC_sampling = 1` through `4`.
iSS samples one event, runs UrQMD on it, and writes the final particles
before sampling the next event. The run uses `number_of_repeated_sampling`
or the existing adaptive event count in sampling modes 2 and 4. Enabled
sampler checks accumulate statistics from each sample before transport.

The mode skips OSCAR, file 14, and per-species sample files. The selected
`use_gzip_format` or `use_binary_format` output contains one record per event,
including empty events. After a disk-output run, the iSS event getters expose
only the last transported event at index 0. Set both format switches to
`0` to retain all transported events in memory instead. UrQMD can still
read or create its initialization table cache (`tables.dat`).

`urqmd_duration` and `urqmd_max_step` are in fm/c and default to `8000`.
`urqmd_seed = -1` (the default) uses iSS's resolved `randomSeed` for the
first transport event. When `randomSeed = -1`, both modules use the same
seed taken from the local system clock. Set a positive `urqmd_seed` to
override the transport seed. UrQMD transport accepts seeds from 1 through
2,147,483,398. The linked cascade uses UrQMD's default physics options and
the supplied burner file's `ctp 22 1.0` setting. It evolves particles in the
frame of the input samples.

C++ callers can link the UrQMD-enabled `libiSS` and include the event wrapper
from `src`. The wrapper serializes access to UrQMD's shared Fortran state.

```cpp
#include "urqmd_transport.h"

UrqmdTransport transport(20.0, 5.0, 123987);
auto result = transport.evolve(sample);  // sample is a std::vector<iSS_Hadron>
// result.particles contains the final event.
```
