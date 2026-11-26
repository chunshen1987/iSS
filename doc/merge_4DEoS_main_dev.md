# Merging 4DEoS into main and main into dev

## 1. `origin/4DEoS` → `merge4D` → `main`

`merge4D` branched from `main`. `origin/4DEoS` added 95 commits. These
brought the 4D EoS (`src/eos_4D.*`, `EOS_to_use = 20`), spin
polarization, the NEOS3D/HRG4D EoS and 4D delta f download scripts, and
changes to the emission function, readers and samplers.

Merge commit `a648273`. One conflict:

- **`src/FSSW.cpp`, `sample_momemtum_from_a_fluid_cell`: rapidity of a
  sampled hadron.** 4DEoS used `y = 0.5*log((E + pz)/(E - pz))` and
  resampled when the result was NaN or Inf. `merge4D` (`8b3e43b`) used
  `y = asinh(pz/mT)`, which cannot be NaN or Inf. Kept `asinh`. All
  other 4DEoS changes to the file merged without conflict.

`main` was then fast-forwarded to `merge4D`. `CLAUDE.md` was added on
`main` as `2725216`.

Build: clean. CI tests: all 8 pass.

## 2. `main` → `dev`

`dev` had 8 commits that were not on `main`, and `main` had 32 that were
not on `dev`. Merge commit `f65385a`. Conflicts in three files:

| File | Conflict | Resolution |
|---|---|---|
| `src/FSSW.h`, `src/FSSW.cpp` (about 25 hunks) | `dev` renamed the flags (`bUSE_*` → `USE_*` as `int`, `bINCLUDE_*_DELTAF` → `INCLUDE_*_DELTAF`), and `main` edited the same lines | Kept `dev`'s names and renamed the `b*` names coming from `main`, including the new `bUSE_OSCAR2013` → `USE_OSCAR2013` |
| `src/FSSW.cpp`, constructor | `dev` reads parameters with defaults and limits `deltaf_max_ratio` to [0, 1]; `main` added `use_OSCAR2013` | Kept `dev`'s constructor and added `USE_OSCAR2013 = getVal("use_OSCAR2013", 0)` |
| `src/FSSW.cpp`, delta f coefficients | `main` uses `surf->visCoeffs` for EoS 20 and `|nB|` for NEoS-BQS, and renamed `bulkvisCoefficients` → `visCoefficients` | Took `main`'s version; code that merged without conflict already used the new name |
| `src/FSSW.cpp`, rapidity | Both sides use `asinh`; `dev` also keeps a NaN check | Kept `dev`'s |
| `src/readindata.cpp`, valid EoS ids | `dev` added 42 (EoS Generator 1D), `main` added 20 | Kept both. `main` labelled 20 "BEST" by mistake; relabelled it "4D NEoS (HRG4D)" |

One resolved hunk dropped a `{` that `dev`'s surrounding code needed (the
`INCLUDE_DIFFUSION_DELTAF` check before `get_deltaf_qmu_coeff` in the
sampling loop). The first build caught it, and the brace was restored.

To check that nothing from `main` was lost, the merged `FSSW.cpp` was
compared with `main`'s, ignoring the flag renames. Every remaining
difference comes from one of `dev`'s own commits: parameter defaults,
the `deltaf_max_ratio` limit and acceptance factor, the general `bulkPi`
fallback, whitespace, and the NaN check.

Build: clean, no warnings. CI tests: all 8 pass.
