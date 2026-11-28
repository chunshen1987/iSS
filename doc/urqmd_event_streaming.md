# Event-by-event iSS and UrQMD transport

With `run_urqmd_in_memory = 1`, iSS now samples one event, runs UrQMD transport, and writes and flushes the final particle record before sampling the next event. This removes the need to hold the full batch of sampled or transported particles in memory during a disk-output run.

The loop supports sampling modes 1–4. It uses `number_of_repeated_sampling` unless the existing adaptive particle target is enabled in mode 2 or 4. It keeps one UrQMD transport instance for the run, preserves empty events as zero-count records, and accumulates enabled sampler diagnostics from the particles before transport.

The output filenames and event record layouts remain `particle_samples.bin` and `particle_samples.gz`. The iSS getters expose the last transported event at index 0 in a disk-output run. When both output switches are disabled, all final events remain available through the getters.

## Reproducibility and trade-offs

The event-first loop changes how random draws are assigned relative to the old batch workflow. Runs with the same seed and settings can be replayed. Sampling tables that depend on the particle species are recalculated for each event to avoid keeping one set per species for the full surface. This trades more sampling work for a smaller memory footprint.

## Verification

Build with `ISS_ENABLE_URQMD=ON` and run `ctest --test-dir build --output-on-failure`. The command-line test uses `randomSeed = 1`, a one-cell surface, and a 20 fm/c UrQMD duration. It checks binary and gzip records, the same first event in one-event and three-event runs, zero-count records, the adaptive count for a zero-yield surface, and pre-transport diagnostics.

Both UrQMD CTest cases passed. A C++ smoke test confirmed that the last of three events matched between disk-output and memory-only runs and that the getters exposed one and three events, respectively. The direct UrQMD output, sampler consistency, and OSCAR2013 checks passed with transport disabled.
