# FSSW benchmark scripts

The scripts behind [doc/FSSW_performance.md](../../doc/FSSW_performance.md).
They need Python 3 with numpy. Run them from the iSS repo root.

| Script | What it does |
|---|---|
| `make_benchmark_surface.py WORK_DIR [--profile P]` | Writes the 300K-cell MUSIC binary surface and its `music_input`. `P` is `constant` (default: ε = 0.26 GeV/fm³, EoS 9), `varying-e` (EoS 9), `nB` (EoS 14) or `BQS` (n_B, n_Q, n_S; EoS 20, needs the downloaded 4D tables) |
| `make_spacelike_cell.py WORK_DIR` | Writes the spacelike, flowing single-cell test surface and its `music_input` |
| `sample_stats.py summarize BIN JSON` | Reduces `particle_samples.bin` to per-species means and errors |
| `sample_stats.py compare A.json B.json` | Prints the z-score of each observable, and the max and rms |
| `run_benchmark.sh WORK_DIR OUT_DIR EXE...` | Times several executables, interleaved and repeated, and keeps a summary of each run |

## Example

```bash
python3 utilities/benchmark/make_benchmark_surface.py /path/to/bench

# build an executable that does not load lib/libiSS.so
mkdir -p /path/to/new/build && cd /path/to/new/build
cmake /path/to/iSS -Dlink_with_lib=OFF && make -j iSS.e && cd /path/to/iSS

utilities/benchmark/run_benchmark.sh /path/to/bench /path/to/out \
    /path/to/old/build/src/iSS.e /path/to/new/build/src/iSS.e
python3 utilities/benchmark/sample_stats.py compare \
    /path/to/out/old_df20_1.json /path/to/out/new_df20_1.json
```

`run_benchmark.sh` labels each executable by the folder three levels up
(`/path/to/old/build/src/iSS.e` → `old`). The samples are written to
`particle_samples.bin` in the repo root and deleted after each run; a
1000-event run of the benchmark surface is about 215 MB. Keep in mind
that `/tmp` is often a RAM-backed tmpfs.
