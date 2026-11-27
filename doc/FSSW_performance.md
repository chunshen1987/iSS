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

### Follow-up: the setup and the yields when T varies

Two follow-up commits address the cost found above:

- **EoS setup.** `regulate_surface_cells` runs over the cells with
  OpenMP. The out-of-range warnings of the 2D lookup are collected and
  printed after the loop, in cell order. For EoS 20 one call,
  `EOS_4D::getThermalVariablesAndDeltafCoeffs`, interpolates all 11
  fields with one set of tilde variables, corner indices and weights,
  without `at()`. It uses the same weights and summation order as
  `FourDLInterp`. On 10⁶ random points, inside and beyond the table, it
  matches the separate calls bit for bit for both δf tables.
- **Yields.** The cell keeps 1/T, and `calculate_dN_analytic` and
  `interpolate_bessel` no longer divide: 1/n comes from a table, the
  table position from a stored 1/Δx, and the post-processing multiplies
  by T instead of dividing by β.

The idea in the previous open item, computing the fugacities e^{μ/T}
once per cell, was also tried. It made no measurable difference
(10.46 → 10.42 s on one thread for `nB`), so it was left out, and
e^{μ/T} is computed as before.

Setup and yields of a 1-event run, δf = 20, best of three for the setup:

| | Setup, before → after | Yields, 1 thread | Yields, 12 threads |
|---|---|---|---|
| Constant ε | 0.38 → 0.38 s | 1.07 → 0.96 s | 0.62 → 0.52 s |
| `varying-e` | 0.38 → 0.38 s | 4.00 → 3.23 s | 0.79 → 0.75 s |
| `nB` | 0.45 → 0.45 s | 4.20 → 3.34 s | 0.83 → 0.80 s |
| `BQS` | 0.62 → 0.49 s (12 threads) | 4.17 → 3.33 s | 0.84 → 0.78 s |

1000 events, before → after the follow-up:

| Profile | δf | 12 threads | 1 thread |
|---|---|---|---|
| constant ε | 20 | 1.94 → 1.84 s | 6.70 → 6.55 s |
| constant ε | 21 | 1.90 → 1.76 s | 6.23 → 6.13 s |
| `varying-e` | 20 | 2.30 → 2.20 s | 10.45 → 9.66 s |
| `varying-e` | 21 | 2.15 → 2.10 s | 9.87 → 9.08 s |
| `nB` | 20 | 2.31 → 2.27 s | 10.63 → 9.74 s |
| `nB` | 21 | 2.21 → 2.16 s | 9.83 → 9.11 s |
| `BQS` | 20 | 2.48 → 2.24 s | 10.51 → 9.44 s |
| `BQS` | 21 | 2.40 → 2.21 s | 10.35 → 9.57 s |

- The yields lose about 20% of their single-thread time, and runs on one
  thread get 7–10% faster when T varies. With 12 threads the yield loop
  was already short, and the gain is 2–10%.
- The setup is not dominated by `regulateEOS`. Running it in parallel
  changes nothing for EoS 9 and 14. For EoS 20 it saves 0.13 s; the rest
  is reading the surface and the 176 MB of 4D tables, and the LRF
  transform.
- The total yields of all 321 species agree with the previous code to
  every printed digit. Removing the divisions changes the rounding, so
  the samples are not byte-identical. For
  each profile the output is the same on 1 and 12 threads.

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
different species, observables and surfaces. The follow-up commits,
against the code before them (seeds 3 and 5), on all four surfaces:
max |z| 1.91–3.03, rms z 0.83–1.11, one of 968 observables above 3.
This does not test the open
mismatch between the HRG4D table and the sampled gas
([EoS4D_validation.md](EoS4D_validation.md)), which is the same before
and after these changes.

## Resonance decays

With `perform_decays = 1`, FSSW decays the sampled resonances after
sampling, in `perform_resonance_feed_down`. This step used to run
serially. On the constant-ε benchmark it took about 1.3 s, 40% of a
12-thread run, and it doubled the peak memory.

The events now decay in parallel:

- `particle_decay` has overloads of `perform_decays` and of the two-,
  three- and n-body routines that take the random generator as an
  argument. The old signatures forward to the shared generator, so the
  `EmissionFunctionArray` path (`MC_sampling` 1–3) is unchanged. The two
  warning counters are atomic.
- Each event decays with its own stream, `Random(key_hi, key_lo, iev)`,
  keyed from two draws of the main generator after sampling, as in the
  sampling itself. The output does not depend on the number of threads.
- The event's list moves into a per-thread queue instead of being
  copied. The daughter list and the buffers of `perform_decays` are
  reused instead of being allocated for every particle.

1000 events, δf = 20, with decays (5.4k primaries and 9.1k final hadrons
per event):

| Surface | Threads | Before | After | Decay step after |
|---|---|---|---|---|
| constant ε | 12 | 3.18 s, 812 MB | 2.07 s, 537 MB | 0.14 s |
| constant ε | 1 | 7.84 s, 640 MB | 7.51 s, 492 MB | 0.83 s |
| `BQS` | 12 | 3.59 s, 811 MB | 2.50 s, 540 MB | 0.14 s |
| `BQS` | 1 | 10.72 s, 646 MB | 10.27 s, 489 MB | 0.81 s |

On 1, 2, 6 and 12 threads the decay step takes 0.90, 0.48, 0.19 and
0.16 s (before reusing the buffers), and the output is the same.

### Validation

- **Statistics.** Final-state samples, before vs after (seeds 3 and 5):
  101 observables of the species that survive the decays, max |z| 1.94
  (constant ε) and 2.44 (`BQS`), rms z 0.86 and 1.05, none above 3.
- **Old path.** `MC_sampling = 2` with decays on a CI cell (20 events,
  518K final hadrons) gives the same file before and after.
- **Buffers.** Reusing the buffers does not change the output.
- **Decays at rest.** `tests/TestResonanceDecayAtRest.cpp`, which CI does
  not run, passes.
- **CI tests.** All 8 pass (they run without decays).
- **Four-momentum.** With the same seed, a run without decays gives
  exactly the primaries that the run with decays decays. Per event, the
  total E changed by 0.1 GeV (median) and at most 0.47 GeV, out of
  43 TeV, and p_z by up to 0.63 GeV. The code before these changes gave
  the same size (0.34 and 0.49 GeV at most). The next section finds the
  cause and fixes it.

### Four-momentum conservation in decays

Each decay routine built the daughters in the mother's rest frame with
total energy `mother->mass`, then boosted them with v = p/E and
γ = 1/√(1 − v²) from the mother's stored E and p. Both steps lose
precision for a fast mother, because `iSS_Hadron` stores E and p as
floats:

- **The velocity was computed in float.** `mother->px/mother->E` divides
  two floats. At large γ, 1 − v is small, and γ = 1/√(1 − v²) inherits
  the rounding of v amplified by about γ². A φ exactly on shell at
  E = 64 GeV lost 0.002 GeV in every decay.
- **The stored mother is off its mass shell.** Rounding E and p to float
  moves √(E² − p²) away from the mass by about γ² times the float
  precision (10⁻³ at γ = 100–300). Daughters that share `mother->mass` in
  the rest frame then carry (M/√(E² − p²))·(E, p) instead of (E, p).

Both errors grow as γ². In the benchmark, resonances at |y| up to about
5 reach γ = 50–100 and more.

The fix (in `particle_decay.cpp`, for the two-, three- and n-body
routines):

- Decay at the mass of the stored four-momentum, `rest_frame_mass`:
  √(E² − p²) in double, which is exact for float inputs. Right at a
  threshold, where it may fall below the daughters' masses, the mother's
  mass is kept. Channel selection still uses `mother->mass`.
- Boost with P/M instead of the velocity, `boost_from_rest_frame`:
  E′ = (E_P·E* + P·p*)/M and p′ = p* + P·[E* + P·p*/(E_P + M)]/M.
  Daughters with total (M, 0) then add up to P, with no 1 − v
  cancellation.

The random draws are unchanged; only the kinematics move, at the γ²·ε
level.

Per decay, ρ, K*, Δ⁺⁺, ω and φ with random directions, 20000 decays of
each species per γ range (100000 per row):

| γ | Mean \|ΔE\|/E before | after | max \|ΔE\| before | after |
|---|---|---|---|---|
| 1–2 | 5.4·10⁻⁸ | 1.1·10⁻⁸ | 9.5·10⁻⁷ GeV | 1.2·10⁻⁷ GeV |
| 10–30 | 1.1·10⁻⁵ | 1.2·10⁻⁸ | 3.9·10⁻³ GeV | 1.9·10⁻⁶ GeV |
| 30–100 | 1.1·10⁻⁴ | 1.3·10⁻⁸ | 0.11 GeV | 4.1·10⁻⁶ GeV |
| 100–300 | 1.1·10⁻³ | 1.3·10⁻⁸ | 3.1 GeV | 1.5·10⁻⁵ GeV |

After the fix, what is left is the float rounding of the stored
daughters, the same at every γ.

Per event on the `BQS` surface (100 events, same seed with and without
decays), the change of the total E went from a median of 0.1 GeV
(at most 0.47 GeV) to 1.5·10⁻⁵ GeV (at most 5.8·10⁻⁵ GeV), and of p_z
from at most 0.63 GeV to 6·10⁻⁵ GeV. That is below the float rounding
of an event, about 0.003 GeV. The final-state samples agree statistically
with the code before the fix (101 observables, max |z| 2.80–2.86, none
above 3), and the decay step still takes 0.14 s.

`tests/TestResonanceDecayAtRest.cpp` now also decays ρ, ω, Δ⁺⁺ and K*
at γ = 3, 30 and 300 and requires the daughters to add up to the
mother's four-momentum to 10⁻⁶ of E. The old code fails it (ρ at
γ = 30: 23.0709 against 23.07 GeV).

About 0.4% of the final particles are still off their mass shell by more
than 1 MeV. Those are high-energy particles and photons, for which
E² − p² computed from the stored floats loses the mass. Their E and p
are each accurate to the float precision; only double-precision storage
in `iSS_Hadron` would change this.

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

- **Decays in CI.** The CI tests run without decays, and
  `tests/TestResonanceDecayAtRest.cpp` is not built by CMake or run in
  CI.
- **The |p| draw** (`MomentumSamplerBase::inverse_CDF`, about 130 ns) is
  a binary search over a table that combines three columns on every
  probe. Precombining the table per (m/T, μ/T) would only pay off when T
  is constant.
- **Setup** (0.4–0.5 s) is serial reading: the surface, and for EoS 20
  the 176 MB of 4D tables. The δf tables are read one float at a time
  into six arrays.
- **Yields when T varies** (3.2 s of CPU). The loop runs once per species
  over all cells. With 12 threads it takes 0.75 s, of which 0.14 s is
  the serial inverse-CDF setup (`RandomVariable1DArray`) per species.
  Flattening the `double**` Bessel table made no difference.
- **GPU.** A CUDA prototype of the yield loop and the case against a
  CUDA or Metal backend for now: [GPU_feasibility.md](GPU_feasibility.md).
- **Minor.** `getCENEOSBQSCoefficients` and
  `get22momNEOSBQSCoefficients` multiply the energy index by the
  energy-axis length instead of the n_B-axis length. Both are 200, so
  the results are not affected.
