"""Check sample diagnostics, OSCAR2013 events and VISH2+1 charge totals."""

import math
import re
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
CELL = (ROOT / "tests/testViscousOneFluidCell3.dat").read_text().split()


def run(directory, *options):
    return subprocess.run(
        [str(EXE), str(ROOT / "iSS_parameters.dat"), ".", "surface.dat", *options],
        cwd=directory, capture_output=True, text=True, timeout=15,
    )


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    result = run(directory, "MC_sampling=0", "calculate_vn=1",
                 "bulk_deltaf_kind=1", "perform_checks=1")
    assert result.returncode == 1, result.stdout + result.stderr
    assert "perform_checks = 1 requires MC_sampling = 4" in result.stdout

    (directory / "iSS_tables").symlink_to(ROOT / "iSS_tables")
    (directory / "music_input").symlink_to(ROOT / "tests/music_input_14")
    cell = CELL.copy()
    cell[4] = "1e-8"
    surface = directory / "surface.dat"
    surface.write_text(" ".join(cell) + "\n")
    result = run(directory, "number_of_repeated_sampling=3", "randomSeed=1",
                 "use_OSCAR_format=1", "use_OSCAR2013=1")
    assert result.returncode == 0, result.stdout + result.stderr
    lines = (directory / "OSCAR.DAT").read_text().splitlines()
    assert lines[2:] == [
        line
        for event in range(1, 4)
        for line in (f"# event {event} out 0", f"# event {event} end 0 impact 0.000")
    ], lines

    cell[4] = "100"
    surface.write_text(" ".join(cell) + "\n")
    result = run(directory, "number_of_repeated_sampling=0",
                 "perform_checks=1", "use_binary_format=0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Sample checks require a surface and events" in result.stdout
    assert not (directory / "checkReconstructedTmunu.dat").exists()

    result = run(directory, "number_of_repeated_sampling=1", "randomSeed=1",
                 "perform_checks=1", "use_binary_format=0")
    assert result.returncode == 0, result.stdout + result.stderr
    reconstruction = (directory / "checkReconstructedTmunu.dat").read_text()
    values = [float(value) for value in " ".join(reconstruction.splitlines()[1:]).split()]
    assert len(values) == 51 and all(map(math.isfinite, values)), reconstruction

    cell[4] = "0"
    surface.write_text(" ".join(cell) + "\n")
    result = run(directory, "number_of_repeated_sampling=1",
                 "perform_checks=1", "use_binary_format=0")
    assert result.returncode == 1, result.stdout + result.stderr
    assert "Sample reconstruction requires a finite, positive surface volume" in result.stdout

    (directory / "checkReconstructedTmunu.dat").write_text("0 0 nan\n" * 17)
    result = subprocess.run([sys.executable, str(ROOT / "tests/TestOutputFiles.py")],
                            cwd=directory, capture_output=True, text=True, timeout=10)
    assert result.returncode == 1 and "nonfinite reconstruction data" in result.stderr, result

    # VISH2+1 rows contain nB, muB and muS, without nQ, nS or muQ.
    (directory / "decdat2.dat").write_text(
        "10 2 0 0 0 0 0.3 0.1 0.15 0 0 0.05 0 0 0 0 0 0 0\n"
    )
    surface.write_text("0 0 0 0\n")
    result = run(directory, "hydro_mode=0", "MC_sampling=0", "regulateEOS=0")
    assert result.returncode == 0, result.stdout + result.stderr
    charges = re.findall(r"nQ\[\d\] = (\S+)", result.stdout)
    assert list(map(float, charges)) == [2., 0., 0.], result.stdout

print("Surface charge, sample check and OSCAR2013 regressions passed.")
