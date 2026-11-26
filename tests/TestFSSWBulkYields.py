"""Compare FSSW bulk yield shifts with smooth Cooper-Frye integration."""

import math
import re
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()


def calculate_yield(directory, mode, kind, ratio=1, bulk=True):
    result = subprocess.run(
        [str(EXE), str(ROOT / "iSS_parameters.dat"), ".",
         "testViscousOneFluidCell3.dat", f"MC_sampling={mode}",
         f"calculate_vn={int(mode == 0)}", f"bulk_deltaf_kind={kind}",
         f"include_deltaf_bulk={int(bulk)}", f"deltaf_max_ratio={ratio}",
         "hydro_mode=1", "regulateEOS=0", "include_deltaf_shear=0",
         "include_deltaf_diffusion=0",
         "number_of_repeated_sampling=0", "use_binary_format=0"],
        cwd=directory, capture_output=True, text=True, timeout=15,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    if mode == 0:
        integrated = directory / "thermal_211_integrated__vndata.dat"
        return float(integrated.read_text().splitlines()[0].split()[1])
    match = re.search(r"Sampling using dN_dy=([^,]+),", result.stdout)
    assert match, result.stdout
    return float(match[1])


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    for source, name in (("music_input_14", "music_input"),
                         ("testViscousOneFluidCell3.dat",
                          "testViscousOneFluidCell3.dat")):
        (directory / name).symlink_to(ROOT / "tests" / source)
    tables = directory / "iSS_tables"
    tables.mkdir()
    for name in ("bin_tables", "deltaf_tables", "pdg-urqmd_v3.3+.dat"):
        (tables / name).symlink_to(ROOT / "iSS_tables" / name)
    (tables / "chosen_particles_urqmd_v3.3+.dat").write_text("211\n")

    fssw_equilibrium = calculate_yield(directory, 4, 0, bulk=False)
    smooth_equilibrium = calculate_yield(directory, 0, 0, bulk=False)
    for kind in (0, 2, 3, 4):
        for ratio in (1, 0.2):
            fssw_shift = calculate_yield(directory, 4, kind, ratio) - fssw_equilibrium
            smooth_shift = (
                calculate_yield(directory, 0, kind, ratio) - smooth_equilibrium
            )
            # The smooth reference integrates pT and rapidity on a separate grid.
            assert math.isclose(fssw_shift, smooth_shift, rel_tol=0.02), (
                kind, ratio, fssw_shift, smooth_shift
            )
        assert calculate_yield(directory, 4, kind, 0) == fssw_equilibrium, kind

print("FSSW bulk yield shifts agree with smooth integration for kinds 0, 2, 3 and 4.")
