# MUSIC freeze-out surface layout

iSS reads MUSIC surfaces with `read_FOdata::read_FOsurfdat_MUSIC`
(`hydro_mode = 2`, full 3+1D) and
`read_FOdata::read_FOsurfdat_MUSIC_boost_invariant` (`hydro_mode = 1`).
The format is selected by `freeze_surface_in_binary` in `music_input`.

## Binary (`freeze_surface_in_binary 1`)

36 floats per cell, plus 38 when `output_vorticity 1`. Older MUSIC
versions wrote 34 floats, without Qn and Sn; iSS does not read those
files any more.

| Index | Quantity |
|---|---|
| 0–3 | τ, x, y, η_s |
| 4–7 | dσ_μ |
| 8–11 | u^μ |
| 12–17 | e, T, μB, μS, μQ, (e + P)/T (in fm⁻¹ units) |
| 18–27 | π^{μν} (10 components) |
| 28 | Π (bulk) |
| 29 | nB |
| 30 | nQ |
| 31 | nS |
| 32–35 | q^μ |
| 36–73 | vorticity and related arrays (only with `output_vorticity 1`) |

MUSIC's default branch (`main`) still writes 34 floats. The 36-float
layout is written by the `main_4D` and `chun_dev_4D` branches.

## Text (`freeze_surface_in_binary 0`)

One cell per line. Columns 0–27 are always present, in the same order
as binary indices 0–27. The rest depend on `music_input`:

| Columns | Present when |
|---|---|
| Π | `Include_Bulk_Visc_Yes_1_No_0 1` |
| nB | `Include_Rhob_Yes_1_No_0 1` |
| nQ, nS | `Include_Rhob_Yes_1_No_0 1` **and** `EOS_to_use 20` (since `c699c7e`) |
| q^μ (4 columns) | `turn_on_baryon_diffusion 1` |

For every EoS other than 20, Qn and Sn are not read and are set to 0,
which matches the text files MUSIC writes today. The text readers do
not read vorticity: `output_vorticity 1` with a text surface is ignored
with a warning.

## Matching change in MUSIC

MUSIC's text writer does not output ρQ and ρS on any branch, including
`main_4D` and `chun_dev_4D`. For EoS 20 text surfaces to carry them,
MUSIC `main_4D`'s `src/evolve.cpp` needs the change below, in all three
places that write a text surface. It follows MUSIC's existing
`turn_on_rhob` checks and uses the same variables as the binary path a
few lines above each one. It has not been compiled against MUSIC.

```diff
--- a/src/evolve.cpp
+++ b/src/evolve.cpp
@@ -737,6 +737,10 @@
                                << fluid_center.Wmunu[ii] << " ";
                     if (DATA.turn_on_bulk) s_file << fluid_center.pi_b << " ";
                     if (DATA.turn_on_rhob) s_file << fluid_center.rhob << " ";
+                    // iSS reads rhoq and rhos after rhob for the 4D EoS
+                    if (DATA.turn_on_rhob && DATA.whichEOS == 20)
+                        s_file << fluid_center.rhoq << " "
+                               << fluid_center.rhos << " ";
                     if (DATA.turn_on_diff)
                         for (int ii = 10; ii < 14; ii++)
                             s_file << std::scientific << std::setprecision(10)
@@ -1037,6 +1041,9 @@
                        << Wetaeta_center << " ";
                 if (DATA.turn_on_bulk) s_file << pi_b_center << " ";
                 if (DATA.turn_on_rhob) s_file << rhob_center << " ";
+                // iSS reads rhoq and rhos after rhob for the 4D EoS
+                if (DATA.turn_on_rhob && DATA.whichEOS == 20)
+                    s_file << rhoq_center << " " << rhos_center << " ";
                 if (DATA.turn_on_diff)
                     s_file << qtau_center << " " << qx_center << " "
                            << qy_center << " " << qeta_center << " ";
@@ -1339,6 +1346,10 @@
                             s_file << fluid_center.pi_b << " ";
                         if (DATA.turn_on_rhob)
                             s_file << fluid_center.rhob << " ";
+                        // iSS reads rhoq and rhos after rhob for the 4D EoS
+                        if (DATA.turn_on_rhob && DATA.whichEOS == 20)
+                            s_file << fluid_center.rhoq << " "
+                                   << fluid_center.rhos << " ";
                         if (DATA.turn_on_diff)
                             for (int ii = 10; ii < 14; ii++)
                                 s_file << std::scientific
```
