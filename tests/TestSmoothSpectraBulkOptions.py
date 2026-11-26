"""Check bulk correction selection for smooth spectra and v_n."""

import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()


def run(directory, kind, *options):
    return subprocess.run(
        [str(EXE), str(ROOT / "iSS_parameters.dat"), ".",
         "testViscousOneFluidCell3.dat", "MC_sampling=0", "calculate_vn=1",
         f"bulk_deltaf_kind={kind}", *options],
        cwd=directory, capture_output=True, text=True, timeout=15,
    )


def check_error(result, expected):
    output = result.stdout + result.stderr
    assert result.returncode == 1, output
    assert expected in output, output
    assert "read in hyper-surface" not in output, output


with TemporaryDirectory() as workdir:
    directory = Path(workdir)
    for kind in (-1, 4.5, 5, 11, 20, 21):
        check_error(run(directory, kind), f"do not support bulk_deltaf_kind = {kind}")
    check_error(run(directory, 20, "MC_sampling=3", "calculate_vn=0"),
                "do not support bulk_deltaf_kind = 20")
    check_error(run(directory, 20, "MC_sampling=4"),
                "calculate_vn requires MC_sampling = 0, 1, 2, or 3")

    for name in ("music_input_14", "testViscousOneFluidCell3.dat"):
        destination = "music_input" if name == "music_input_14" else name
        (directory / destination).symlink_to(ROOT / "tests" / name)
    tables = directory / "iSS_tables"
    tables.mkdir()
    for name in ("EOS_tables", "deltaf_tables", "bin_tables",
                 "pdg-urqmd_v3.3+.dat"):
        (tables / name).symlink_to(ROOT / "iSS_tables" / name)
    (tables / "chosen_particles_urqmd_v3.3+.dat").write_text("211\n")

    spectra_file = directory / "dN_pTdpTdphidy.dat"
    result = run(directory, 20, "include_deltaf_bulk=0")
    assert result.returncode == 0, result.stdout + result.stderr
    uncorrected_spectra = spectra_file.read_bytes()
    assert (directory / "thermal_211_vndata.dat").stat().st_size > 0
    for kind in range(5):
        result = run(directory, kind)
        assert result.returncode == 0, result.stdout + result.stderr
        assert spectra_file.read_bytes() != uncorrected_spectra, kind

print("Smooth spectra reject unsupported options and apply kinds 0 through 4.")
