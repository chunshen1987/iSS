# Development notes

Reports from the November 2026 work that merged the 4D EoS branch and
cleaned up the issues found in review.

| Report | What it covers |
|---|---|
| [merge_4DEoS_main_dev.md](merge_4DEoS_main_dev.md) | Merging `origin/4DEoS` into `main` and `main` into `dev`: conflicts and how they were resolved |
| [code_review_fixes.md](code_review_fixes.md) | Findings from the code review of the merged `dev`, the fixes and how each was checked |
| [EoS4D_validation.md](EoS4D_validation.md) | End-to-end tests of `EOS_to_use = 20` and the open mismatch between the HRG4D table and the sampled hadron gas |
| [surface_text_format.md](surface_text_format.md) | Column layout of MUSIC freeze-out surfaces (binary and text) and the matching MUSIC change for Qn, Sn |
| [PR6_review.md](PR6_review.md) | Review of pull request #6 (n-body decays), which targets `XSCAPE` |
| [FSSW_performance.md](FSSW_performance.md) | Benchmark of FSSW on a 300K-cell surface (δf 20, 21), the profile, and the four optimization commits on `optimization` |
| [GPU_feasibility.md](GPU_feasibility.md) | Whether to port FSSW to CUDA or Metal: what a GPU could replace, a CUDA prototype of the yield loop, and why no backend was written |

Unless a report says otherwise, "the CI tests" means the 8 scripts
`tests/RunAutoTest{Ideal1-4,Viscous1-4}.sh`, run from the repo root
after `make install`. They all use `MC_sampling = 4` (FSSW),
`perform_decays = 0` and text surfaces with EoS 9, 91, 12 or 14.
