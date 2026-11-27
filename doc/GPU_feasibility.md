# Is a GPU port of FSSW worthwhile?

**Not now.** On the 300K-cell benchmark of
[FSSW_performance.md](FSSW_performance.md), even a complete port of the
FSSW yields and sampling would take a 1000-event run from 2.2 s to about
0.8 s. It would mean a second implementation of most of the sampler,
which CI cannot test. The part that ports easily, the yields, saves at
most about 0.4 s. No CUDA or Metal backend was written.

Machine: Intel i7-9750H (6 cores, 12 threads) and an NVIDIA GTX 1650
Mobile (Turing, sm_75, 16 SMs, 4 GB), CUDA 13.5.

## What a GPU could replace

A 12-thread run on the `varying-e` surface, `bulk_deltaf_kind = 20`,
1000 events, 2.20 s in all (code at `3fc92b5`):

| Part | Time | On a GPU? |
|---|---|---|
| Setup: reading the surface and tables, `regulateEOS`, LRF transform | 0.38 s | No: disk and parsing |
| Yields of 321 species in 300K cells | 0.75 s, of which 0.14 s is the serial inverse-CDF setup | Yes, a regular loop |
| Sampling 5.4M hadrons (2.0 µs each on one thread) | 0.91 s | Possible, with divergent rejection loops |
| Writing the samples, exit | 0.15 s | No |

Making the yields and the sampling free would leave about 0.55 s.

## Yield prototype

[`utilities/benchmark/gpu_yield_bench.cu`](../utilities/benchmark/gpu_yield_bench.cu)
is a standalone prototype of the yield loop for the 22-moment δf
(`calculate_dN_analytic`, three Bessel functions, up to 10 terms for the
light species). It runs the same function with OpenMP and as a CUDA
kernel (one thread per cell, a loop over species, species in constant
memory). It uses 300K random cells spanning the benchmark's T and μ, and
369 species (iSS samples 321 of them).

| | Time |
|---|---|
| CPU, 12 threads | 0.42–0.45 s |
| GPU FP64 kernel | 0.19 s |
| GPU FP32 kernel | 0.023 s |
| Copying all yields back, FP64 (886 MB) / FP32 | 0.11 s / 0.055 s |
| CUDA start-up (`cudaFree(0)`), once per process | 0.14–0.17 s |
| Uploading the cells and the Bessel table | 0.003 s |

Agreement with the CPU, maximum relative difference:

| | Per cell | Per species total |
|---|---|---|
| FP64 | 1.0e-15 | 3.6e-16 |
| FP32 | 2.8e-6 | 6.3e-7 |

- **FP64 gains nothing.** The kernel is only 2.3 times faster than 12
  CPU threads, because this GPU computes FP64 at 1/32 of its FP32 rate.
  The copy back and the start-up use up the gain.
- **FP32 is 20 times faster** than the CPU in the kernel. This only pays
  off if the per-species inverse CDF and the choice of cells for each
  hadron also run on the GPU, so that the 300K yields of each species
  never leave it. The net gain is then about 0.6 s minus the 0.15 s
  start-up, so about 0.4 s of 2.2 s.
- FP32 precision would be enough for sampling, but the CDF sums over
  300K cells would need FP64 or compensated accumulation.

## Why not port the sampler too

The sampling is the other half of the time that a GPU could remove. A
port would have to reimplement, for each backend:

- the exact flux sampler and the boson and fermion momentum samplers
  (`MomentumSamplerBase` and its CDF tables);
- the δf kinds 0, 1, 11, 20 and 21, the bulk quadrature, diffusion δf
  and PCE;
- the boosts and the η_s and y ranges;
- a counter-based random number generator. The GPU events could not
  reproduce the per-event `mt19937` streams of the CPU, so the two
  backends would only agree statistically.

The best case is about 2.2 s → 0.8 s on this benchmark. Against that:

- iSS is a small step. If the samples go on to UrQMD or SMASH, that
  step usually costs much more.
- Used as a library (JETSCAPE/XSCAPE), every process pays the CUDA
  start-up, which is about as long as the 12-thread yield loop.
- **Metal has no FP64.** A Metal backend would differ numerically from
  the CPU path, and it can only be built and tested on macOS.
- The CI runners have no GPU, so a GPU backend would go untested in CI,
  and every physics change would have to be made in three places.
- Resonance decays (`perform_decays = 1`) run on the CPU. They now run
  in parallel over events (0.14 s of a 2.1 s run on 12 threads).

## When to revisit

- **Throughput runs**, with one single-thread iSS job per core. A run
  then costs about 9.7 s per core, and a GPU would gain far more
  relative to that.
- **Much larger surfaces or event counts**, with no afterburner.
- **GPU nodes** that are a standard part of the workflow.

If so, start with one optional CUDA piece behind a build flag: the FP32
yields with the inverse CDF and cell selection on the GPU. The momentum
sampling stays on the CPU, and the CPU path remains the reference.
