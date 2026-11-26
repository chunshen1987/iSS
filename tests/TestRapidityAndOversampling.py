"""Check 3+1D midrapidity spectra and the FSSW event cap."""

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
        [str(EXE), str(ROOT / "iSS_parameters.dat"), ".", "surface.dat",
         "regulateEOS=0", "bulk_deltaf_kind=1", "include_deltaf_bulk=0",
         "include_deltaf_shear=0", "use_binary_format=0", *options],
        cwd=directory, capture_output=True, text=True, timeout=15,
    )


def read_spectrum(directory, eta, hydro_mode=2):
    cell = CELL.copy()
    cell[3] = str(eta)
    (directory / "surface.dat").write_text(" ".join(cell) + "\n")
    result = run(directory, f"hydro_mode={hydro_mode}", "MC_sampling=0",
                 "calculate_vn=1")
    assert result.returncode == 0, result.stdout + result.stderr
    output = (directory / "thermal_211_vndata.dat").read_text().splitlines()
    return [tuple(map(float, row.split()[:3])) for row in output]


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    (directory / "music_input").symlink_to(ROOT / "tests/music_input_14")
    tables = directory / "iSS_tables"
    tables.mkdir()
    for name in ("bin_tables", "pdg-urqmd_v3.3+.dat"):
        (tables / name).symlink_to(ROOT / "iSS_tables" / name)
    (tables / "chosen_particles_urqmd_v3.3+.dat").write_text("211\n")

    central = read_spectrum(directory, 0)
    shifted = read_spectrum(directory, 2)
    assert shifted == read_spectrum(directory, -2)
    assert read_spectrum(directory, 0, 1) == read_spectrum(directory, 2, 1)

    # The cell is static, muQ = 0, and its normal is purely temporal.
    # The spectrum at y = 0 scales as cosh(eta_s) f0(mT*cosh(eta_s)/T).
    particle_table = ROOT / "iSS_tables/pdg-urqmd_v3.3+.dat"
    pion = next(row.split() for row in particle_table.read_text().splitlines()
                if row.split()[:2] == ["211", "Pion(+)"])
    pion_mass = float(pion[2])
    temperature = float(CELL[13])*0.197327053
    cosh_eta = math.cosh(2.)
    for (pT, _, unshifted_yield), (_, _, shifted_yield) in zip(central[:8], shifted[:8]):
        mt_over_t = math.hypot(pT, pion_mass)/temperature
        ratio = cosh_eta*math.expm1(mt_over_t)/math.expm1(mt_over_t*cosh_eta)
        assert math.isclose(shifted_yield/unshifted_yield, ratio, rel_tol=2e-5)

    result = run(directory, "MC_sampling=3")
    assert result.returncode == 1, result.stdout + result.stderr
    assert "MC_sampling = 3 requires a boost-invariant surface" in result.stdout

    cell = CELL.copy()
    for surface_weight, expected_events in (("0.001", 100), ("0.000001", 100),
                                            ("0", 1)):
        cell[4] = surface_weight
        (directory / "surface.dat").write_text(" ".join(cell) + "\n")
        result = run(directory, "MC_sampling=4",
                     "sample_upto_desired_particle_number=1",
                     "number_of_particles_needed=100000",
                     "maximum_sampling_events=100", "randomSeed=1")
        assert result.returncode == 0, result.stdout + result.stderr
        match = re.search(r"number of repeated sampling = (\d+)", result.stdout)
        assert match, result.stdout
        assert int(match[1]) == expected_events, (surface_weight, match[1])

print("3+1D midrapidity spectra and low-yield oversampling checks passed.")
