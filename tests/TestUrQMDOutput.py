"""Check direct UrQMD restart records against iSS's OSCAR samples."""

import math
import os
from pathlib import Path
from struct import pack, unpack
import subprocess
import sys
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "iSS.e").resolve()
SPECIES = {
    211: (101, 2, 1),
    2212: (1, 1, 1),
    -2212: (-1, -1, -1),
    213: (104, 2, 1),
    311: (106, -1, 0),
    -311: (-106, 1, 0),
}


def read_oscar(path):
    lines = iter(path.read_text().splitlines()[3:])
    for line in lines:
        event, count = map(int, line.split()[:2])
        particles = []
        for _ in range(count):
            _, pdg, *values = next(lines).split()
            particles.append((int(pdg), tuple(map(float, values))))
        yield event, particles


def check_file14(path, oscar, mode):
    lines = iter(path.read_text().splitlines())
    for event, particles in oscar:
        if not particles:
            continue
        assert next(lines).startswith("UQMD   version:")
        assert next(lines) == "projectile:  (mass, char)  197  79   target:  (mass, char)  197  79 "
        betas = next(lines)
        assert len(betas) == 69
        assert list(map(float, (betas[i:i+11] for i in (36, 47, 58)))) == [
            0., 0.9907499, -0.9907499]
        assert next(lines).startswith("impact_parameter")
        assert next(lines).startswith("equation_of_state")
        assert int(next(lines)[7:16]) == event
        for tag in ("op",)*3 + ("pa",)*4 + ("pvec:",):
            assert next(lines).startswith(tag)
        assert list(map(int, next(lines).split())) == [len(particles), 0]
        assert list(map(int, next(lines).split())) == [0]*8

        # The file-14 reader counts baryons at the front of the event.
        particles.sort(key=lambda particle: abs(SPECIES[particle[0]][0]) > 55)
        for pdg, coordinates in particles:
            # Legacy samplers use doubles for OSCAR and iSS_Hadron floats
            # for the in-memory UrQMD records.
            if mode in (1, 3):
                coordinates = unpack("9f", pack("9f", *coordinates))
            px, py, pz, energy, mass, x, y, z, time = coordinates
            record = next(lines)
            assert len(record) == 233
            assert tuple(map(int, record[144:185].split())) == (*SPECIES[pdg], 0, 0, 99)
            values = [float(record[i:i+16]) for i in range(0, 144, 16)]
            expected = [0., x-px/energy*time, y-py/energy*time, z-pz/energy*time,
                        energy, px, py, pz, mass]
            for value, target in zip(values, expected):
                assert math.isclose(value, target, rel_tol=1e-7, abs_tol=1e-10)
            decay, formation, cross_section_factor = (
                float(record[i:i+16]) for i in (185, 201, 217))
            assert math.isclose(formation, time, rel_tol=1e-7)
            assert cross_section_factor == 0.
            if pdg == 213:
                assert formation < decay < 1e34
            else:
                assert decay == 1e34
    assert next(lines, None) is None


with TemporaryDirectory() as workdir:
    work = Path(workdir)
    (work / "music_input").symlink_to(ROOT / "tests/music_input_14")
    tables = work / "iSS_tables"
    tables.mkdir()
    for name in ("bin_tables", "OSCAR_header.txt", "pdg-urqmd_v3.3+.dat", "z_exp_m_z.dat"):
        (tables / name).symlink_to(ROOT / "iSS_tables" / name)
    (tables / "chosen_particles_urqmd_v3.3+.dat").write_text(
        "\n".join(map(str, SPECIES)) + "\n")
    surface = (ROOT / "tests/testViscousOneFluidCell3.dat").read_text().split()
    surface[4] = "500"
    (work / "surface.dat").write_text(" ".join(surface) + "\n")

    for mode in (1, 2, 3, 4):
        args = [
            str(EXE), str(ROOT / "iSS_parameters.dat"), ".", "surface.dat",
            "hydro_mode=1", "regulateEOS=0", "bulk_deltaf_kind=1",
            "include_deltaf_bulk=0", "include_deltaf_shear=0", "randomSeed=1",
            "sample_upto_desired_particle_number=0", "number_of_repeated_sampling=3",
            "use_OSCAR_format=1", f"MC_sampling={mode}",
        ]
        for option in ("use_UrQMD_format=1", "use_UrQMD_format=0"):
            result = subprocess.run(
                [*args, option], cwd=work, env=dict(os.environ, OMP_NUM_THREADS="2"),
                capture_output=True, text=True, timeout=30,
            )
            assert result.returncode == 0, result.stdout + result.stderr
            if option.endswith("=1"):
                assert not (work / "OSCAR.DAT").exists()
                assert not (work / "particle_samples.bin").exists()

        oscar = list(read_oscar(work / "OSCAR.DAT"))
        assert len(oscar) == 3 and all(particles for _, particles in oscar)
        assert {pdg for _, particles in oscar for pdg, _ in particles} == set(SPECIES)
        check_file14(work / "fort.14", oscar, mode)
        (work / "OSCAR.DAT").unlink()
        (work / "fort.14").unlink()

print("Direct UrQMD output checks passed for all sampling modes.")
