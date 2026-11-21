"""Build an offline replay ZIP from the row evidence and its registered source dependencies."""

import argparse
from hashlib import sha256
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ROW_FOLDERS = (
    "row_04", "row_04_compact", "row_04_identity", "row_04_conjugates", "row_04_large",
    "row_30", "row_30_free", "row_30_branch", "row_30_memory", "row_30_arrows",
    "row_30_selective", "row_30_autonomous", "row_30_route",
)
SOURCES = (
    "LICENSE",
    "code/LICENSE",
    "code/defect_particle_classes/PROOF.md",
    "model/requirements.txt",
    "model/__init__.py",
    "model/executable_model/__init__.py",
    "model/executable_model/carrier.py",
    "model/executable_model/symmetry.py",
    "model/evaluation/particles.json",
    "model/candidate_laws/__init__.py",
    "model/candidate_laws/contrast_memory_spherical_v1/__init__.py",
    "model/candidate_laws/contrast_memory_spherical_v1/candidate.py",
    "model/candidate_laws/contrast_memory_spherical_v1/law.py",
    "model/candidate_laws/contrast_memory_spherical_v1/sectors.py",
    "model/candidate_laws/contrast_memory_spherical_v1/port_response.py",
    "model/candidate_laws/contrast_memory_spherical_v1/common.py",
    "model/candidate_laws/contrast_memory_spherical_v1/collar.py",
    "model/candidate_laws/contrast_memory_spherical_v1/PROOF.md",
    "model/candidate_laws/contrast_memory_spherical_v1/SPEC.md",
    "model/candidate_laws/interval_echo/__init__.py",
    "model/candidate_laws/interval_echo/candidate.py",
    "model/candidate_laws/interval_echo/clock.py",
    "model/candidate_laws/interval_echo/currents.py",
    "model/candidate_laws/interval_echo/CURRENTS.md",
    "model/candidate_laws/interval_echo/SPEC.md",
    "model/experiments/interval_echo_encounter/README.md",
    "model/experiments/interval_echo_gauge/__init__.py",
    "model/experiments/interval_echo_gauge/__main__.py",
    "model/experiments/interval_echo_gauge/README.md",
)


def members():
    paths = [ROOT / name for name in SOURCES]
    paths.extend((ROOT / "model/executable_model/carriers").glob("*.npz"))
    paths.extend(
        HERE / name
        for name in ("__init__.py", "README.md", "replay.py", "build_archive.py")
    )
    for subdir in ROW_FOLDERS:
        folder = HERE / subdir
        if (
            not (folder / ("receipt.json.gz" if subdir == "row_30_autonomous" else "results.json")).is_file()
            or not (folder / "README.md").is_file()
        ):
            raise FileNotFoundError(f"Incomplete experiment in {folder.name}")
        paths.extend(
            p
            for p in folder.iterdir()
            if p.is_file() and (p.suffix in {".md", ".py", ".json", ".npz", ".gz"} or p.name == "SHA256SUMS")
        )
    submission = HERE / "submission"
    paths.extend(p for p in submission.iterdir() if p.is_file())
    files = {p.relative_to(ROOT).as_posix(): p.read_bytes() for p in paths}
    files["README.md"] = (
        "# Interval Echo evidence archive\n\n"
        "The [research note and proposed table entries](model/candidate_laws/interval_echo/openwave/README.md) "
        "give the experiments, results, limits and reproduction commands.\n\n"
        "From this directory, install the pinned dependencies and replay with Python 3.12.\n\n"
        "```sh\npython -m pip install -r model/requirements.txt\n"
        "python -B -m model.candidate_laws.interval_echo.openwave.replay\n```\n\n"
        "SHA256SUMS lists every other source and evidence file. Licensing is described in LICENSE.\n"
    ).encode()
    files["SHA256SUMS"] = "".join(
        f"{sha256(data).hexdigest()}  {name}\n" for name, data in sorted(files.items())
    ).encode()
    return files


def build(path):
    files = members()
    with ZipFile(path, "w") as archive:
        for name, data in sorted(files.items()):
            info = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)
    return len(files)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, default=HERE / "interval_echo_evidence.zip")
    path = cli.parse_args().output
    print(f"Archived {build(path)} files in {path.name}")
