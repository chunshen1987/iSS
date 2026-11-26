"""Check unambiguous and colliding MUSIC binary record sizes."""

import math
import re
import struct
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
CELL = list(map(float, (ROOT / "tests/testViscousOneFluidCell3.dat")
                .read_text().split()))


def pack_record(values):
    return struct.pack(f"={len(values)}f", *values)


def run_surface(directory, mode, vorticity, payload, *options):
    (directory / "music_input").write_text(
        "EOS_to_use 20\n"
        "Include_Bulk_Visc_Yes_1_No_0 1\n"
        "Include_Rhob_Yes_1_No_0 1\n"
        "turn_on_baryon_diffusion 0\n"
        "freeze_surface_in_binary 1\n"
        f"output_vorticity {vorticity}\n\n"
    )
    (directory / "surface.dat").write_bytes(payload)
    return subprocess.run(
        [str(EXE), "iSS_parameters.dat", str(directory), "surface.dat",
         f"hydro_mode={mode}", "regulateEOS=0", "MC_sampling=0", *options],
        cwd=ROOT, capture_output=True, text=True, timeout=10,
    )


def check_cells(result, count):
    assert result.returncode == 0, result.stdout + result.stderr
    counts = re.findall(r"total number of cells: (\d+)", result.stdout)
    assert counts == [str(count)] * 2, result.stdout
    charges = re.findall(r"nQ\[\d\] = ([-+\d.eE]+)", result.stdout)
    assert list(map(float, charges)) == [5000.*count, 2000.*count, -1000.*count]


def check_error(result, message):
    output = result.stdout + result.stderr
    assert result.returncode == 1, output
    assert message in output, output
    assert "total number of cells:" not in output, output


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    for mode in (1, 2):
        for vorticity in (0, 1):
            suffix = [0.] * (38*vorticity)
            legacy = pack_record(CELL + [0.]*4 + suffix)
            current = pack_record(CELL + [0.04, -0.02] + [0.]*4 + suffix)
            gcd = math.gcd(len(legacy), len(current))
            legacy_collision = legacy * (len(current)//gcd)
            current_count = len(legacy)//gcd
            current_collision = current * current_count

            check_cells(run_surface(directory, mode, vorticity, current), 1)
            check_error(run_surface(directory, mode, vorticity, legacy),
                        "matches the legacy")
            for payload in (legacy_collision, current_collision):
                check_error(run_surface(directory, mode, vorticity, payload),
                            "Ambiguous MUSIC binary surface")
            check_cells(run_surface(directory, mode, vorticity, current_collision,
                                    "music_binary_surface_columns=36"),
                        current_count)
            check_error(run_surface(directory, mode, vorticity, current + b"\0"),
                        "expected whole")
            check_error(run_surface(directory, mode, vorticity, current,
                                    "music_binary_surface_columns=36.5"),
                        "music_binary_surface_columns must be")

print("MUSIC binary surface checks passed for both hydro modes.")
