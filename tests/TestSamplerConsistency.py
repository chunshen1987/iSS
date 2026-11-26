"""Exercise yield and spectrum consistency, reconstruction, and I/O errors."""

import math
import re
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
CELL = (ROOT / "tests/testViscousOneFluidCell3.dat").read_text().split()


def setup(directory, name):
    work = directory / name
    tables = work / "iSS_tables"
    tables.mkdir(parents=True)
    for table in ("pdg-urqmd_v3.3+.dat", "deltaf_tables", "bin_tables",
                  "OSCAR_header.txt", "z_exp_m_z.dat"):
        (tables / table).symlink_to(ROOT / "iSS_tables" / table)
    (tables / "chosen_particles_urqmd_v3.3+.dat").write_text("211\n")
    (work / "music_input").symlink_to(ROOT / "tests/music_input_14")
    write_surface(work, [CELL])
    return work


def write_surface(directory, rows):
    (directory / "surface.dat").write_text(
        "".join(" ".join(row) + "\n" for row in rows)
    )


def run(work, *options, error=None):
    result = subprocess.run(
        [str(EXE), str(ROOT / "iSS_parameters.dat"), ".", "surface.dat",
         "regulateEOS=0", "use_binary_format=0", "include_deltaf_shear=0",
         "randomSeed=1", "number_of_repeated_sampling=1", *options],
        cwd=work, capture_output=True, text=True, timeout=15,
    )
    output = result.stdout + result.stderr
    if error:
        assert result.returncode == 1 and error in output, output
        assert "Program executed normally" not in output, output
    else:
        assert result.returncode == 0, output
    return output


def get_yield(output):
    return float(re.search(r"Sampling using dN_dy=([^,]+)", output)[1])


with TemporaryDirectory() as workdir:
    directory = Path(workdir)

    work = setup(directory, "cap")
    for kind in (1, 20, 21):
        options = (f"bulk_deltaf_kind={kind}", "deltaf_max_ratio=0")
        assert get_yield(run(work, *options, "include_deltaf_bulk=1")) == get_yield(
            run(work, *options, "include_deltaf_bulk=0")
        )

    work = setup(directory, "pce")
    (work / "iSS_tables/EOS_particletable.dat").write_text("1\n1 211\n")
    (work / "decdat2.dat").write_text(
        "10 5000 0 0 0 0 0.3 0 0.15 0 0 0.05 0 0 0 0 0 0 0\n"
    )
    (work / "surface.dat").write_text("0 0 0 0\n")
    cases = []
    for mu in (0., 0.1):
        (work / "decdat_mu.dat").write_text(f"0 {mu}\n")
        yield_ = get_yield(run(
            work, "hydro_mode=0", "include_deltaf_bulk=0",
            "y_LB=-0.5", "y_RB=0.5", "use_OSCAR_format=1", "use_OSCAR2013=1"
        ))
        momenta = [row.split()[5:9]
                   for row in (work / "OSCAR.DAT").read_text().splitlines()
                   if not row.startswith("#")]
        assert len(momenta) > 100
        cases.append((yield_, momenta[:100]))
    assert cases[1][0] > cases[0][0] and cases[1][1] != cases[0][1], cases

    work = setup(directory, "signed")
    reverse = CELL.copy()
    reverse[4] = str(-float(CELL[4]))
    options = ("calculate_vn=1", "bulk_deltaf_kind=1", "include_deltaf_bulk=0")
    for mode in (0, 2, 3):
        sampling = (f"MC_sampling={mode}",
                    f"hydro_mode={1 if mode == 3 else 2}",
                    "y_LB=0", "y_RB=0.002")
        write_surface(work, [CELL, reverse])
        run(work, *options, *sampling, "use_pos_dN_only=0")
        flow = work / "thermal_211_integrated__vndata.dat"
        total = float(flow.read_text().splitlines()[0].split()[1])
        assert abs(total) < 1e-8, (mode, total)
        run(work, *options, *sampling, "use_pos_dN_only=1")
        assert float(flow.read_text().splitlines()[0].split()[1]) > 0.

    work = setup(directory, "volume")
    cells = []
    for normal in ("100", "300"):
        cell = CELL.copy()
        cell[4] = normal
        cells.append(cell)
    tensors = []
    for order in (cells, cells[::-1]):
        write_surface(work, order)
        run(work, "perform_checks=1", "include_deltaf_bulk=0")
        rows = (work / "checkReconstructedTmunu.dat").read_text().splitlines()
        tensors.append(tuple(map(float, rows[1].split()[:2])))
    assert all(math.isclose(a, b, rel_tol=1e-7)
               for a, b in zip(*tensors)), tensors

    work = setup(directory, "output")
    for filename, options in (
        ("particle_samples.bin", ("use_binary_format=1",)),
        ("particle_samples.gz", ("use_gzip_format=1",)),
        ("OSCAR.DAT", ("use_OSCAR_format=1", "use_OSCAR2013=1")),
    ):
        path = work / filename
        path.mkdir()
        (path / "keep").touch()
        run(work, *options, error=f"Sample output {filename}: cannot open file")
        assert (path / "keep").exists()

    if Path("/dev/full").exists():
        work = setup(directory, "failed_write")
        for filename, options in (
            ("particle_samples.bin", ("use_binary_format=1",)),
            ("particle_samples.gz", ("use_gzip_format=1",)),
            ("OSCAR.DAT", ("use_OSCAR_format=1", "use_OSCAR2013=1")),
        ):
            (work / filename).symlink_to("/dev/full")
            run(work, *options, error=f"Sample output {filename}: write failed")

    work = setup(directory, "hrg")
    tables = work / "iSS_tables/EOS_tables"
    tables.mkdir()
    table = tables / "HRGNEOS_BQS-urqmd_v3.3+.dat"
    table.touch()
    error = f"Invalid 2D HRG table iSS_tables/EOS_tables/{table.name}"
    run(work, "regulateEOS=1", error=error + ": expected at least two")
    table.write_text("# grid\n0.01 invalid\n")
    run(work, "regulateEOS=1", error=error + ": invalid numeric row")

print("FSSW, signed spectra, reconstruction, and output error checks passed.")
