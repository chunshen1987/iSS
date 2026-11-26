# Code review of the merged dev branch: findings and fixes

A `/code-review` of `dev` after the `main` merge (`f65385a`) reported 10
findings and 3 side notes. One finding (9) was already in `dev`'s code.
The others came with the 4D EoS work from `main`. None came from the
conflict resolution. The review read the code without running it, and
none of the findings were covered by the CI tests.

All of them are fixed on `dev`.

## Findings

| # | Where | Problem | Fix | Commit | How it was checked |
|---|---|---|---|---|---|
| 1 | `eos_4D.cpp`, `FourDLInterp` | Each grid index was limited to `N-1` and its neighbour `N` read. At the T edge this threw `std::out_of_range`. At the μB edge it read the next row of the table and returned **μB = −0.70 for large positive nB** | Limit the lower index to `N-2`; the weight clamps to 1 at the edge | `1865f55` | Probe program querying the table before and after the fix: the throws and the wrong sign are gone, and values inside the table are identical ([EoS4D_validation.md](EoS4D_validation.md)) |
| 2 | `readindata.cpp`, `regulate_surface_cells` | A text surface with `output_vorticity 1` read `vorticity_arr[24..33]` from an empty vector | Ignore `output_vorticity` for text surfaces, with a warning; only binary readers fill vorticity | `1865f55` | Text cell with `output_vorticity 1`: warning, normal exit |
| 3 | `FSSW.cpp`, both delta f coefficient sites | EoS 20 always took `surf->visCoeffs`, which is empty unless `regulateEOS = 1` and `bulk_deltaf_kind` is 20 or 21 | Use the per-cell 4D coefficients only for kind 20 or 21; exit with a clear error when EoS 20 and kind 20 or 21 run without `regulateEOS = 1` | `1865f55` | EoS 20 with `regulateEOS = 0`: "EOS_to_use = 20 with bulk_deltaf_kind = 20 requires regulateEOS = 1" |
| 4 | `FSSW.cpp`, constructor | EoS 20 always loaded the 4D delta f file, and exited when it was missing, even for kinds that never use it; otherwise it kept a second, unused copy of a 90 MB table | Removed `FSSW::eos_4d_` and the load (the coefficients come from `read_FOdata`) | `1865f55` | EoS 20 with RTA delta f (kind 1) and the 22-moment file removed: runs |
| 5 | `readindata.cpp` | `regulateEOS = 1` is the default and regulated every EoS. Anything other than 12 or 14 used the μ = 0 `HRGEOS_PST` table, so finite-μB EoS such as BEST (17) lost their chemical potentials | Regulate only EoS 9, 91, 12, 14 and 20; warn and keep the surface as it is for any other | `acf71c2` | EoS 17 with `regulateEOS = 1`: warning, μ kept |
| 6 | `iSS.cpp`, `shell` | Since `6674945` (2022), `MC_sampling = 0` skipped everything, so `calculate_vn = 1` no longer wrote the smooth spectra and v_n | Call `generate_samples()` when `MC_sampling != 0` or `calculate_vn == 1`; `EmissionFunctionArray` does not sample when `MC_sampling = 0` | `e207b5d` | Flowing ideal cell: `dN_pTdpTdphidy.dat` and the `thermal_*_vndata.dat` files written, no particle samples |
| 7 | `readindata.cpp`, surface readers | Binary surfaces need 36 floats per cell (Qn, Sn at 30–31), so 34-float files from older MUSIC no longer load; text surfaces set Qn = Sn = 0 | Older binary layout: dropped on purpose, no compatibility needed. Text: read Qn, Sn after Bn for EoS 20 ([surface_text_format.md](surface_text_format.md)) | `c699c7e` | EoS 20 text cell with Qn = 0.4 nB: the surface's nQ comes through (2000 instead of 0) |
| 8 | `readindata.cpp`, `read_in_HRG_EOS` | EoS table path built from `particle_table_path_` instead of `table_path_` | Use `table_path_`, like the other EoS tables | `1865f55` | Build and CI tests (both paths are `iSS_tables` in `main.cpp`) |
| 9 | `FSSW.cpp`, constructor | `deltaf_max_ratio` was read without a default, and `getVal(name)` calls `exit(1)` when the key is missing (JETSCAPE/XSCAPE and older parameter files) | Default 1.0, which also keeps the old acceptance factor (1 + δf)/2 | `1865f55` | Build and CI tests |
| 10 | `iSS.cpp`, `perform_checks` | The sampled charges were summed as (B, S, Q) and the surface totals as (B, Q, S), so the check printed Q next to S | Sum the sampled charges as (B, Q, S) | `1865f55` | BQS CI test (Ideal4) log |

Fixing #10 exposed the text-reader half of #7. With the order corrected,
the BQS CI cell showed surface Q = 0 against sampled Q = +2005. The old
order had compared surface Q (0) with sampled S (−0.11) and looked
fine.

## Side notes

All fixed in `cb31b03`.

| Where | Problem | Fix |
|---|---|---|
| `readindata.cpp`, both MUSIC readers | Cells with u·dσ < 0 were dropped without a message | Still dropped (a sampler cannot draw a negative number of particles), but counted and reported: "dropped N cells with u.dsigma < 0 (x% of the total \|u.dsigma\|)" |
| `FSSW.cpp`, `getParticleInfo` | Returned index 0 for an unknown particle id, so its charges were counted as particle 0's | Returns −1. `computeAvgTotalEnergyMomentum` skips such particles and warns with the count |
| `FSSW.cpp`, `computeAvgTotalEnergyMomentum` | `sqrt` of a variance that rounding can make slightly negative gave NaN | Variance clamped at 0; also returns early when there are no events |

Checked with a two-cell surface whose second cell has its normal vector
flipped: "dropped 1 cells with u.dsigma < 0 (50% of the total
|u.dsigma|)", and finite charge averages and errors.

## Verification summary

- Build: clean, no warnings, after every change.
- CI tests: all 8 pass after every change.
- The 4D EoS paths (#1, #3, #4, #7) were run with the HRG4D tables
  downloaded by `iSS_tables/EOS_tables/download_HRG4D.sh` and
  `iSS_tables/deltaf_tables/urqmd/download_NEoS4D_deltafCoeffs.sh`. Those
  tables are listed in `.gitignore` (`7ec0e14`).
