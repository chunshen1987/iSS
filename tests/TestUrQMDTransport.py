"""Exercise the in-memory UrQMD mode and check that it skips event files."""

import gzip
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
RECORD = struct.Struct("=i9f")


def run(args, work, count=3):
    result = subprocess.run(
        args, cwd=work, capture_output=True, text=True, timeout=90,
        env=dict(os.environ, OMP_NUM_THREADS="2", URQMD_TAB="tables.dat"),
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"UrQMD processed {count} samples" in result.stdout, result.stdout
    seeds = re.search(r"iSS random seed = (\d+), UrQMD transport seed = (\d+)",
                      result.stdout)
    assert seeds, result.stdout
    return tuple(map(int, seeds.groups()))


def check_binary_events(data, expected_events):
    offset = 0
    for _ in range(expected_events):
        count, = struct.unpack_from("=i", data, offset)
        offset += 4
        assert count > 0
        for _ in range(count):
            pid, *values = RECORD.unpack_from(data, offset)
            offset += RECORD.size
            assert pid != 0
            assert all(map(math.isfinite, values))
            assert values[1] == 20., values
    assert offset == len(data)


def check_gzip_events(path, expected_events):
    with gzip.open(path, "rt") as stream:
        for _ in range(expected_events):
            count = int(next(stream))
            assert count > 0
            for _ in range(count):
                pid, *values = next(stream).split()
                assert int(pid) != 0 and len(values) == 9
                assert all(math.isfinite(float(value)) for value in values)
                assert float(values[1]) == 20.
        assert next(stream, None) is None


with TemporaryDirectory() as workdir:
    work = Path(workdir)
    (work / "music_input").symlink_to(ROOT / "tests/music_input_14")
    tables = work / "iSS_tables"
    tables.mkdir()
    for name in ("bin_tables", "pdg-urqmd_v3.3+.dat", "z_exp_m_z.dat"):
        (tables / name).symlink_to(ROOT / "iSS_tables" / name)
    (tables / "chosen_particles_urqmd_v3.3+.dat").write_text("211\n213\n")
    surface = (ROOT / "tests/testViscousOneFluidCell3.dat").read_text().split()
    surface[4] = "40"
    (work / "surface.dat").write_text(" ".join(surface) + "\n")

    # Legacy samplers normally write per-species samples. Enable every
    # intermediate format to verify that in-memory transport takes priority.
    for mode in (1, 2, 3, 4):
        args = [
            str(EXE), str(ROOT / "iSS_parameters.dat"), ".", "surface.dat",
            "hydro_mode=1", "regulateEOS=0", "bulk_deltaf_kind=1",
            "include_deltaf_bulk=0", "include_deltaf_shear=0", "randomSeed=1",
            "y_LB=-0.5", "y_RB=0.5", "sample_upto_desired_particle_number=0",
            "number_of_repeated_sampling=3", f"MC_sampling={mode}",
            "run_urqmd_in_memory=1", "urqmd_duration=20", "urqmd_max_step=5",
            "output_samples_into_files=1", "store_samples_in_memory=0",
            "use_UrQMD_format=1", "use_OSCAR_format=1", "use_binary_format=1",
            f"use_gzip_format={int(mode == 2)}",
        ]
        assert run(args, work) == (1, 1)
        assert not list(work.glob("fort.*"))
        assert not list(work.glob("samples_*.dat"))
        assert not (work / "OSCAR.DAT").exists()
        if mode == 2:
            final = work / "particle_samples.gz"
            check_gzip_events(final, 3)
            lines = gzip.decompress(final.read_bytes()).splitlines(keepends=True)
            first_event = b"".join(lines[:1 + int(lines[0])])
        else:
            final = work / "particle_samples.bin"
            data = final.read_bytes()
            check_binary_events(data, 3)
            count, = struct.unpack_from("=i", data)
            first_event = data[:4 + count*RECORD.size]
        # The first event must be independent of how many events follow it.
        assert run([*args, "number_of_repeated_sampling=1"], work, count=1) == (1, 1)
        actual = final.read_bytes()
        if mode == 2:
            actual = gzip.decompress(actual)
        assert actual == first_event
        final.unlink()

    assert run([*args, "use_binary_format=0", "use_gzip_format=0"], work) == (1, 1)
    assert not list(work.glob("particle_samples.*"))

    # Parameter files older than this setting still inherit the iSS seed
    # when they omit it. An explicit positive seed overrides inheritance.
    parameters = Path(args[1]).read_text().splitlines()
    (work / "parameters.dat").write_text(
        "\n".join(line for line in parameters if not line.startswith("urqmd_seed")) + "\n")
    args[1] = "parameters.dat"
    assert run([*args, "urqmd_seed=17"], work) == (1, 17)
    iss_seed, urqmd_seed = run([*args, "randomSeed=-1"], work)
    assert iss_seed == urqmd_seed and 1 <= iss_seed <= (1 << 30)
    output = work / "particle_samples.bin"
    clock_seeded_events = output.read_bytes()
    assert run([*args, f"randomSeed={iss_seed}"], work) == (iss_seed, iss_seed)
    assert output.read_bytes() == clock_seeded_events

    assert run([*args, "perform_checks=1"], work) == (1, 1)
    assert (work / "checkReconstructedTmunu.dat").is_file()
    assert (work / "check_211_spectra.dat").is_file()

    # Empty events still occupy one record each. The adaptive count in FSSW
    # falls back to one event when the surface has zero thermal yield.
    surface[4] = "0"
    (work / "surface.dat").write_text(" ".join(surface) + "\n")
    assert run(args, work) == (1, 1)
    assert output.read_bytes() == struct.pack("=3i", 0, 0, 0)
    assert run([*args, "sample_upto_desired_particle_number=1"], work, count=1) == (1, 1)
    assert output.read_bytes() == struct.pack("=i", 0)
    assert run([*args, "MC_sampling=2", "sample_upto_desired_particle_number=1"],
               work, count=1) == (1, 1)
    assert output.read_bytes() == struct.pack("=i", 0)

print("UrQMD event streaming and shared seed checks passed.")
