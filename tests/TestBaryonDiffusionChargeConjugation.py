"""Check conjugate diffusion corrections in spectra and FSSW sampling."""

import re
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
CELL = (ROOT / "tests/testViscousOneFluidCell3.dat").read_text().split()


def run(directory, mode, sign, diffusion):
    pid = sign*2212
    chosen = directory / "iSS_tables/chosen_particles_urqmd_v3.3+.dat"
    chosen.write_text(f"{pid}\n")
    cell = CELL.copy()
    cell[5] = "1500"  # spatial surface normal along the diffusion current
    cell[14] = str(sign*1.5)  # muB in fm^-1
    cell[-1] = str(sign*0.1)  # baryon density
    cell.extend(["0", str(sign*0.02), "0", "0"])
    (directory / "surface.dat").write_text(" ".join(cell) + "\n")
    result = subprocess.run(
        [str(EXE), str(ROOT / "iSS_parameters.dat"), ".", "surface.dat",
         f"MC_sampling={mode}", f"calculate_vn={int(mode == 0)}",
         "regulateEOS=0", "bulk_deltaf_kind=1", "include_deltaf_bulk=0",
         "include_deltaf_shear=0", f"include_deltaf_diffusion={int(diffusion)}",
         "use_OSCAR_format=1", "use_OSCAR2013=1",
         "number_of_repeated_sampling=1", "randomSeed=1"],
        cwd=directory, capture_output=True, text=True, timeout=15,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    if mode == 0:
        return (directory / f"thermal_{pid}_vndata.dat").read_bytes()

    yield_match = re.search(r"Sampling using dN_dy=([^,]+),", result.stdout)
    assert yield_match, result.stdout
    momenta = []
    for line in (directory / "OSCAR.DAT").read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split()
        assert int(fields[9]) == pid, line
        momenta.append(tuple(map(float, fields[5:9])))
    assert momenta, result.stdout
    return float(yield_match[1]), momenta


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    music_input = (ROOT / "tests/music_input_14").read_text()
    (directory / "music_input").write_text(
        music_input.replace("turn_on_baryon_diffusion  0",
                            "turn_on_baryon_diffusion  1")
    )
    tables = directory / "iSS_tables"
    tables.mkdir()
    for name in ("bin_tables", "deltaf_tables", "pdg-urqmd_v3.3+.dat"):
        (tables / name).symlink_to(ROOT / "iSS_tables" / name)

    for mode in (0, 4):
        proton = run(directory, mode, 1, True)
        antiproton = run(directory, mode, -1, True)
        baseline = run(directory, mode, 1, False)
        assert proton == antiproton, f"Diffusion breaks charge conjugation in mode {mode}"
        if mode == 4:
            assert proton[0] != baseline[0], "Diffusion did not change the yield"
            assert proton[1] != baseline[1], "Diffusion did not change the samples"
        else:
            assert proton != baseline, "Diffusion did not change the spectrum"

print("Diffusion changes both spectra and FSSW yields symmetrically under charge conjugation.")
