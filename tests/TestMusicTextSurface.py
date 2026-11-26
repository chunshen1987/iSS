"""Check MUSIC text input with a built iSS executable and no 4D tables."""

import re
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
CELL = (ROOT / "tests/testViscousOneFluidCell3.dat").read_text().split()
CHARGED_CELL = CELL + ["0.04", "-0.02"]


def run_surface(directory, mode, rows, eos=20, diffusion=False):
    (directory / "music_input").write_text(
        f"EOS_to_use {eos}\n"
        "Include_Bulk_Visc_Yes_1_No_0 1\n"
        "Include_Rhob_Yes_1_No_0 1\n"
        f"turn_on_baryon_diffusion {int(diffusion)}\n"
        "freeze_surface_in_binary 0\n"
    )
    (directory / "surface.dat").write_text(
        "\n".join(" ".join(row) for row in rows)
    )
    return subprocess.run(
        [str(EXE), "iSS_parameters.dat", str(directory), "surface.dat",
         f"hydro_mode={mode}", "regulateEOS=0", "MC_sampling=0"],
        cwd=ROOT, capture_output=True, text=True, timeout=10,
    )


def check_charges(result, expected):
    assert result.returncode == 0, result.stdout + result.stderr
    charges = re.findall(r"nQ\[\d\] = ([-+\d.eE]+)", result.stdout)
    assert list(map(float, charges)) == expected, result.stdout


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    for mode in (1, 2):
        check_charges(run_surface(directory, mode, [CHARGED_CELL, CHARGED_CELL]),
                      [10000., 4000., -2000.])
        check_charges(run_surface(directory, mode, [CELL], eos=14),
                      [5000., 0., 0.])
        diffusion = ["0.01", "0.02", "0.03", "0.04"]
        check_charges(
            run_surface(directory, mode, [CHARGED_CELL + diffusion],
                        diffusion=True),
            [5000., 2000., -1000.],
        )

        for malformed in (CELL, CELL + ["0.04"], CELL + ["bad", "0"],
                          CHARGED_CELL + ["0"]):
            result = run_surface(directory, mode,
                                 [CHARGED_CELL, malformed, CHARGED_CELL])
            output = result.stdout + result.stderr
            assert result.returncode == 1, output
            assert "surface.dat at line 2: expected 32 numeric columns" in output
            assert "requires nQ and nS after nB" in output

        result = run_surface(directory, mode, [CELL + diffusion],
                             diffusion=True)
        output = result.stdout + result.stderr
        assert result.returncode == 1, output
        assert "surface.dat at line 1: expected 36 numeric columns" in output

print("MUSIC text surface checks passed for both hydro modes.")
