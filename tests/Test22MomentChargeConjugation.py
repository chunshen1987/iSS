"""Check conjugate 22-moment yields and sampled momenta with a fixed seed."""

import re
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
CELL = (ROOT / "tests/testViscousOneFluidCell3.dat").read_text().split()


def sample(directory, particle_id, density):
    tables = directory / "iSS_tables"
    (tables / "chosen_particles_urqmd_v3.3+.dat").write_text(f"{particle_id}\n")
    (directory / "surface.dat").write_text(
        " ".join(CELL[:-1] + [str(density)]) + "\n"
    )
    result = subprocess.run(
        [str(EXE), str(ROOT / "iSS_parameters.dat"), ".", "surface.dat",
         # Use the fixture's exact zero chemical potentials at nB = 0.
         "MC_sampling=4", f"regulateEOS={int(density != 0.)}", "bulk_deltaf_kind=20",
         "include_deltaf_bulk=1", "include_deltaf_shear=0",
         "use_OSCAR_format=1", "use_OSCAR2013=1",
         "number_of_repeated_sampling=1", "randomSeed=1"],
        cwd=directory, capture_output=True, text=True, timeout=10,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    yield_match = re.search(r"Sampling using dN_dy=([^,]+),", result.stdout)
    assert yield_match, result.stdout
    momenta = []
    for line in (directory / "OSCAR.DAT").read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split()
        assert int(fields[9]) == particle_id, line
        momenta.append(tuple(map(float, fields[5:9])))
    assert momenta, result.stdout
    return float(yield_match[1]), momenta


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    (directory / "music_input").symlink_to(ROOT / "tests/music_input_14")
    tables = directory / "iSS_tables"
    tables.mkdir()
    for name in ("EOS_tables", "deltaf_tables", "pdg-urqmd_v3.3+.dat"):
        (tables / name).symlink_to(ROOT / "iSS_tables" / name)

    for particle_id, density in ((2212, 0.1), (2212, -0.1), (2212, 0.0),
                                 (2112, 0.1), (3122, 0.1)):
        yield_, momenta = sample(directory, particle_id, density)
        anti_yield, anti_momenta = sample(directory, -particle_id, -density)
        assert yield_ == anti_yield, (particle_id, density, yield_, anti_yield)
        assert momenta == anti_momenta, (particle_id, density)

print("22-moment yields and sampled momenta match under charge conjugation.")
