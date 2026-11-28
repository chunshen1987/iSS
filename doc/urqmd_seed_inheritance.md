# Shared iSS and UrQMD random seed

In-memory UrQMD now inherits the resolved iSS `randomSeed` for its first transport event. This is the default when `urqmd_seed` is omitted or set to `-1`. A positive `urqmd_seed` overrides the inherited value.

With `randomSeed = -1`, iSS derives a positive seed from the local system clock and passes that same value to UrQMD. The coupled run logs both seeds. To replay a clock-seeded run, set `randomSeed` to the logged iSS seed and keep the other run settings fixed. UrQMD's supported seed range remains 1 through 2,147,483,398. The direct `fort.14` writer retains its `osc2u` conversion convention.

## Verification

- The coupled tests verified fixed-seed inheritance for all four sampling modes, an explicit UrQMD override, and inheritance from an older parameter file that omits `urqmd_seed`. Replaying a clock-seeded run with its logged seed produced byte-identical particle output. Run with `ctest --test-dir build --output-on-failure` after building with `ISS_ENABLE_URQMD=ON`.
- All three RNG unit tests passed, including replay from a clock-derived seed. On this system, the bundled doctest required `-DDOCTEST_CONFIG_NO_POSIX_SIGNALS` when building the test executable.
- The direct UrQMD output and sampler consistency regression tests passed with `ISS_ENABLE_URQMD=OFF`.
