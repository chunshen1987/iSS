"""Exercise 4D table validation with a small, self-contained grid."""

import struct
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
GRID = [-1., -1., -1., 0.1, 2., 2., 2., 0.2, 2., 2., 2., 2.]


def pack(values):
    return struct.pack(f"={len(values)}f", *values)


def run(directory, kind=1, *options):
    return subprocess.run(
        [str(EXE), str(ROOT / "iSS_parameters.dat"), ".", "surface.dat",
         "MC_sampling=4", f"bulk_deltaf_kind={kind}",
         "number_of_repeated_sampling=0", "use_binary_format=0", *options],
        cwd=directory, capture_output=True, text=True, timeout=10,
    )


def check_error(result, filename, reason):
    output = result.stdout + result.stderr
    assert result.returncode == 1, output
    assert f"{filename}: {reason}" in output, output
    assert "Regulate local temperature" not in output, output


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    tables = directory / "iSS_tables"
    eos = tables / "EOS_tables/HRG4D"
    coefficients = tables / "deltaf_tables/urqmd"
    eos.mkdir(parents=True)
    coefficients.mkdir(parents=True)
    particle_list = "pdg-urqmd_v3.3+.dat"
    (tables / particle_list).symlink_to(ROOT / "iSS_tables" / particle_list)
    (tables / "chosen_particles_urqmd_v3.3+.dat").write_text("211\n")
    music_input = (ROOT / "tests/music_input_14").read_text()
    (directory / "music_input").write_text(
        music_input.replace("EOS_to_use 14", "EOS_to_use 20")
    )
    cell = (ROOT / "tests/testViscousOneFluidCell3.dat").read_text().strip()
    (directory / "surface.dat").write_text(cell + " 0.04 0.0\n")

    for stem, value in (("p", 0.05), ("t", 0.15), ("mub", 0.), ("muq", 0.),
                        ("mus", 0.)):
        (eos / f"HRG_{stem}_b.dat").write_bytes(pack(GRID + [value]*16))
    pressure = eos / "HRG_p_b.dat"
    valid_pressure = pressure.read_bytes()
    for payload, reason in (
        (b"", "truncated grid header"),
        (valid_pressure[:-1], "expected 64 payload bytes, found 63"),
        (valid_pressure + b"\0", "expected 64 payload bytes, found 65"),
        (pack([-2.] + GRID[1:] + [0.05]*16),
         "grid header differs from the temperature table"),
    ):
        pressure.write_bytes(payload)
        check_error(run(directory), pressure.name, reason)
    pressure.write_bytes(valid_pressure)

    temperature = eos / "HRG_t_b.dat"
    valid_temperature = temperature.read_bytes()
    temperature.write_bytes(pack(GRID[:-1] + [1.] + [0.15]*16))
    check_error(run(directory), temperature.name, "invalid grid origin")
    temperature.write_bytes(valid_temperature)

    for regulate in (0, 1):
        result = run(directory, 20, f"regulateEOS={regulate}",
                     "include_deltaf_bulk=0", "include_deltaf_shear=0")
        assert result.returncode == 0, result.stdout + result.stderr
        assert "Read in deltaf Coeffs for EoS 4D" not in result.stdout

    for shear in (0, 1):
        result = run(directory, 20, f"include_deltaf_bulk={1 - shear}",
                     f"include_deltaf_shear={shear}")
        check_error(result, "NEoS4D_22mom_deltafCoeff.bin", "cannot open file")

    for kind, name, count in ((20, "22mom", 6), (21, "CE", 3)):
        table = coefficients / f"NEoS4D_{name}_deltafCoeff.bin"
        valid_coefficients = pack([1.] * (16*count))
        table.write_bytes(valid_coefficients[:-1])
        check_error(run(directory, kind), table.name, "expected ")
        table.write_bytes(valid_coefficients)
        result = run(directory, kind)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "Done reading deltaf Coeffs." in result.stdout

print("4D table validation checks passed.")
