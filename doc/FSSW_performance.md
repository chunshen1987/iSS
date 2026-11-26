# FSSW performance: benchmark and optimizations

A benchmark of `MC_sampling = 4` (FSSW) with `bulk_deltaf_kind` 20 and 21
on a 300K-cell surface, the profile of where the time went, and four
commits on `optimization` that make the run **5.5 times faster**:

| | Before (`39c2e26`) | After (`b66c806`), 12 threads | After, 1 thread |
|---|---|---|---|
| `bulk_deltaf_kind = 20` | 10.90 s | 1.96 s | 6.82 s |
| `bulk_deltaf_kind = 21` | 10.26 s | 1.90 s | — |
| Peak memory (kind 20) | 364 MB | 399 MB | 392 MB |

Machine: Intel i7-9750H laptop (6 cores, 12 threads), 6 GB RAM.

The surface generators, the timing script and the statistical comparison
are in [utilities/benchmark/](../utilities/benchmark/README.md).

## Benchmark

### Surface

The repo has no surface this large, so the benchmark uses a synthetic
3+1D MUSIC binary surface (36 floats per cell,
[surface_text_format.md](surface_text_format.md)), written by
`utilities/benchmark/make_benchmark_surface.py`:

- 300,204 cells: a disc of radius 8 fm on a 0.25 fm grid (3228 cells),
  times 93 η_s slices from −4.6 to 4.6 (Δη_s = 0.1).
- τ_f(r) = 2 + 7 √(1 − (r/R)²) fm, with |dτ_f/dr| capped at 2. Cells
  with |∇τ_f| > 1 are spacelike: 43% of the cells and 44% of the total
  u·dσ (18,831 fm³). That is more than a typical MUSIC surface, so it is
  a hard case for the flux sampler.
- **ε = 0.26 GeV/fm³ in every cell**, μ_B = 0. With EoS 9 and
  `regulateEOS = 1`, every cell then has the same T.
- Flow u^r = 0.9 (r/R)(1 + 0.1 cos 2φ) exp(−η_s²/18) with Gaussian noise
  of 0.02 per component; u^η = 0.
- π^{xx} ≈ π^{yy} ≈ 0.08 (ε + P) with a cos 2φ modulation and 30%
  per-cell noise, τ²π^{ηη} = −(π^{xx} + π^{yy}); `regulate_Wmunu` makes
  it transverse and traceless. Π ≈ −0.015 GeV/fm³ with 30% noise.

`music_input`: EoS 9, bulk on, ρ_B off, diffusion off, binary surface.

### Run

```bash
./iSS.e iSS_parameters.dat <work_dir> surface.dat \
    bulk_deltaf_kind=20 number_of_repeated_sampling=1000 \
    randomSeed=1 use_binary_format=1
```

Each event has about 5,400 hadrons (5.4M per run, 215 MB of binary
output). The timings above are 3 interleaved repeats of each binary; the
spread was below 1%. "Before" is a build of `39c2e26` linked without
`libiSS.so` (`-Dlink_with_lib=OFF`), so that it cannot load a newer
installed library. It reproduces the original output byte for byte.

## Where the time went (before, kind 20)

`perf` is not available on this machine, so a copy of the code was
instrumented with `rdtsc` timers and counters:

| Phase | Time | Notes |
|---|---|---|
| Per-species yields (321 species × 300K cells) | 4.4 s | Does not shrink with the number of events |
| Momentum sampling and δf rejection | 4.9 s | 3.6 tries per accepted hadron, ~250 ns per try |
| Picking the cell (binary search of the CDF) | 0.8 s | |
| Building the hadron records | 0.6 s | |
| Reading the surface, regulation, LRF transform | 0.5 s | |
| CDF construction | 0.1 s | |

Of the 19.4M momentum tries, 9.2M (47%) were rejected by the flux factor,
4.8M (25%) by the δf factor and 5.4M accepted. Each rejection redrew |p|
from the table, which costs 130 ns.

## Changes

| Commit | Change | Output | How it was checked |
|---|---|---|---|
| `74ba878` | `Stopwatch` uses `steady_clock`. It used `clock()`, the CPU time of the whole process, which adds up the time of all threads | Unchanged | Build; output identical to `39c2e26` |
| `b64c2fa` | Yield loop: reuse `calculate_dN_analytic` while (T, μ) repeat; read a compact 72-byte copy of each cell instead of ~250 bytes and a heap vector; OpenMP over cells | **Bit-identical** | Output identical to `39c2e26` at 1 and 12 threads |
| `48edc78` | Exact flux sampler (below); δf bound of 1 when no δf term is on; skip NaN \|p\| | New random stream | Statistical comparisons (below); CI tests |
| `b66c806` | Parallel events: one RNG stream per event, one momentum sampler per thread, explicit-generator overloads, exact cache key | New random stream, same for any thread count | Same output at 1, 2, 4, 6 and 12 threads; statistical comparisons; CI tests |

### Yields (`b64c2fa`)

`calculate_dN_analytic` depends only on the species, T and μ. At constant
ε and μ_B = 0 every cell matches, and 96M calls become a few hundred.
After that the loop is limited by memory traffic, about 11 GB/s on this
laptop, so it now reads a compact array built once in the constructor.
The diffusion, PCE and bulk-quadrature paths still read the full cell.
On a surface where T varies, the cache helps less, but there the loop is
limited by computation and gains more from the OpenMP threads.

### Exact flux sampler (`48edc78`)

The old acceptance was p·dσ/(p⁰(|dσ⁰| + |dσ⃗|)) times (1 + δf)/2. The flux
bound does not depend on v = |p|/p⁰, so slow hadrons (most of the heavy
ones) in flowing cells were rejected often, and every rejection redrew
|p|.

With a = dσ⁰ and b = |dσ⃗| in the LRF, the cut flux max(0, a + v b cos θ)
averages over directions to

- h(v) = a for v b ≤ a,
- h(v) = (a + v b)²/(4 v b) otherwise,

which grows with v up to h_max = h(1). The sampler now keeps |p| with
probability h(v)/h_max, which is always the case for an LRF-timelike
cell, and then draws only the direction against the bound a + v b. This
samples the same distribution exactly. The δf factor still restarts from
|p|, because its angular average depends on |p|.

Tries per hadron fell from 3.6 to 1.9, and |p| draws from 19.4M to
10.2M. With every δf term off, the acceptance had still been divided by
1 + `deltaf_max_ratio`, so ideal runs threw away half of all tries; the
bound is now 1 in that case.

### Parallel events (`b66c806`)

Within each species, the events are sampled in parallel.

- Each event has its own `mt19937`, seeded with `std::seed_seq` from two
  draws of the main generator and the event index. A fixed `randomSeed`
  gives the same samples for any number of threads, and repeated calls
  still continue the main stream.
- The hadron counts are drawn serially, in event order, from the same
  generators as before.
- Each thread has its own `MomentumSamplerShell`, because the samplers
  cache the last (m − μ)/T. The cache is now reused only for an
  identical key. With the old 1e-8 tolerance, a sample could depend on
  which particle the thread drew before.
- `MomentumSamplerShell`, `MomentumSamplerBase` and
  `RandomVariable1DArray` have overloads that take the generator. The
  old signatures forward to them.

Scaling, kind 20 (same output at every thread count):

| Threads | 1 | 2 | 4 | 6 | 12 |
|---|---|---|---|---|---|
| Wall time | 6.82 s | 4.01 s | 2.61 s | 2.19 s | 1.98 s |

## Validation

- **Bit-identity.** `b64c2fa` reproduces the output of `39c2e26`.
  `b66c806` gives the same output at 1–12 threads.
- **Statistics on the benchmark surface.** For 12 species
  (π^±, π⁰, K^±, p, p̄, Λ, Ξ, Ω, φ, ω), the comparison covered the
  multiplicity, ⟨p_T⟩, ⟨p_T²⟩, v₁, v₂ (|y| < 1), ⟨y²⟩, ⟨E⟩, ⟨p_x⟩,
  ⟨p_y⟩ and ⟨p_x²⟩, 121 observables per comparison
  (`utilities/benchmark/sample_stats.py`). Old vs new agree:
  for kind 20, max |z| = 2.3 and rms z = 0.92. For kind 21, pairs with
  independent seeds give rms z = 0.79 and 0.92. Two baseline runs with
  different seeds give rms z = 1.00, the noise level. One kind-21 pion
  ⟨p_x⟩ value came out at z = 3.3; pooled over 3 old and 4 new seeds the
  difference is 2.2σ. It was singled out for being the largest, so that
  is consistent with noise.
- **Spacelike single cell.** One cell with lab-frame |dσ⃗|/dσ⁰ = 3 and
  u = (0.5, 0.3, 0) (`utilities/benchmark/make_spacelike_cell.py`),
  with and without δf, sampled with 1500 events
  (1.3e7 hadrons per run) by the old and new code. The first version of
  the new sampler failed this test: kaon ⟨p_T⟩ was 5.5σ low, because it
  skipped the |p| acceptance when v b ≤ a < b. With that fixed, every
  observable agrees within errors.
- **CI tests.** All 8 pass.

## Surfaces with varying T and chemical potentials

The constant-ε surface is the easy case for the caches, so the same
geometry, flow and viscous fields were also run with three thermodynamic
profiles (`make_benchmark_surface.py --profile`). In each, ε follows a
smooth profile from about 0.15 to 0.47 GeV/fm³ (clipped to 0.10–0.60)
with 10% noise per cell, so no two neighbouring cells match exactly.

| Profile | EoS | Densities | What `regulateEOS = 1` does per cell |
|---|---|---|---|
| `varying-e` | 9 | none | T and P from the 2D HRG table |
| `nB` | 14 | n_B up to 0.16 fm⁻³, two humps at η_s = ±2.5 (baryon stopping) | T, P and μ_B, μ_S, μ_Q from the 2D NEoS-BQS table |
| `BQS` | 20 | as `nB`, with n_Q = 0.4 n_B (5% noise) and n_S around 0 (±0.01 fm⁻³) | T, P, three μ's and the δf coefficients from the 4D HRG tables: 11 4D interpolations |

The `BQS` surface carries total B = 1048, Q = 419 and S = 0.009. Its
samples have p/p̄ ≈ 107/16 and K⁺/K⁻ ≈ 183/101, so the chemical potentials
from the tables reach the sampler. All tilde variables stay inside the
4D grid.

### Timings (1000 events)

| Profile | δf | Before | After, 12 threads | After, 1 thread |
|---|---|---|---|---|
| `varying-e` | 20 | 11.71 s | 2.25 s | 10.54 s |
| `varying-e` | 21 | 11.05 s | 2.15 s | 9.91 s |
| `nB` | 20 | 11.65 s | 2.29 s | 10.50 s |
| `nB` | 21 | 11.03 s | 2.19 s | 9.79 s |
| `BQS` | 20 | 11.42 s | 2.45 s | 10.32 s |
| `BQS` | 21 | 11.53 s | 2.37 s | 10.28 s |

With 12 threads the speedup is 4.7–5.2 times, a little less than the 5.5
times at constant ε. On one thread it is only about 1.1 times (1.6 times
at constant ε). The 12-thread and 1-thread outputs are identical for
every profile.

Where the single-thread time goes (δf = 20; a 1-event run gives setup
plus yields):

| | Constant ε | `varying-e` | `BQS` |
|---|---|---|---|
| Setup: reading, `regulateEOS`, LRF transform | 0.38 s | 0.39 s | 0.61 s |
| Yields of 321 species | 1.07 s | 3.96 s | 4.16 s |
| Sampling the hadrons | 5.4 s | 6.2 s | 5.4 s |

- The yield cache stops helping when T varies, and the yield loop is back
  to about 4 s of computation. That part runs in parallel (about 0.7 s
  on 12 threads), so the parallel speedup holds.
- The momentum samplers recompute their cached CDF terms for nearly
  every hadron, but that costs little: sampling takes about as long as
  at constant ε.
- The 4D lookups add about 0.2 s to the serial setup. With 12 threads the
  setup (0.6–0.8 s) is the largest single part of a 2.4 s run.

### Validation

Before vs after, the same 121 observables as above, for each profile and
δf kind:

| Profile | δf 20: max \|z\|, rms z | δf 21: max \|z\|, rms z |
|---|---|---|
| `varying-e` | 2.54, 1.13 | 2.46, 0.98 |
| `nB` | 2.31, 0.87 | 3.22, 1.04 |
| `BQS` | 2.97, 1.13 | 3.36, 1.14 |

Two of the 726 observables have |z| > 3, about what chance gives (≈2):
Ξ ⟨y²⟩ for `nB`, δf 21, and π⁰ ⟨p_y⟩ for `BQS`, δf 21. They are
different species, observables and surfaces. This does not test the open
mismatch between the HRG4D table and the sampled gas
([EoS4D_validation.md](EoS4D_validation.md)), which is the same before
and after these changes.

## Things to know

- **Same seed, different events.** A fixed `randomSeed` still gives
  reproducible results, but not the same events as `39c2e26`. To keep
  existing outputs byte for byte, stop at `b64c2fa` (about 8 s on this
  benchmark).
- **Threads.** OpenMP uses every core by default. Set `OMP_NUM_THREADS`
  on shared nodes, or when iSS runs as a library inside JETSCAPE/XSCAPE.
  Without OpenMP the code runs serially. Compilers without OpenMP warn
  about the pragmas, as `spin_polarization.cpp` already does (CMake skips
  OpenMP for AppleClang).
- **Constant ε is the easy case.** When T varies from cell to cell, the
  yield cache rarely hits, and single-thread runs gain little (see the
  section above). The parallel speedup holds.

## Open items

- **Resonance decays** (`perform_decays = 1`) still run serially after
  sampling, and the CI tests do not cover them. They could now take most
  of the run time. The same per-event streams would let them run in
  parallel.
- **The |p| draw** (`MomentumSamplerBase::inverse_CDF`, about 130 ns) is
  a binary search over a table that combines three columns on every
  probe. Precombining the table per (m/T, μ/T) would only pay off when T
  is constant.
- **Setup** (0.4–0.8 s) is serial and now a quarter to a third of the
  run. The cells are independent, so `regulate_surface_cells` could run
  in parallel. For EoS 20, the 11 `FourDLInterp` calls per cell each
  recompute the same tilde variables and the 16 weights, and read the
  table through bounds-checked `at()`. One call that interpolates all
  fields at once would do that work once.
- **Yields when T varies** (about 4 s of CPU per run). Per-cell factors
  shared by all species, such as e^{μ_B/T}, e^{μ_S/T}, e^{μ_Q/T} and
  1/T, could be computed once instead of once per species.
- **Minor.** `getCENEOSBQSCoefficients` and
  `get22momNEOSBQSCoefficients` multiply the energy index by the
  energy-axis length instead of the n_B-axis length. Both are 200, so
  the results are not affected.
