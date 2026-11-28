"""Extract the high-l Plik-lite arrays and its official likelihood self-test from ESA."""

import argparse
import hashlib
import io
import json
import tarfile
from pathlib import Path
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent
URL = "https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_Likelihood_Data-baseline_R3.00.tar.gz"
PREFIX = "baseline/plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik/clik/"
EXTERNAL = PREFIX + "lkl_0/_external/"
NAMES = (
    "blmin.dat",
    "blmax.dat",
    "bweight.dat",
    "cl_cmb_plik_v22.dat",
    "c_matrix_plik_v22.dat",
)


def extract(stream):
    want = {EXTERNAL + k: k for k in NAMES} | {
        PREFIX + "check_param": "check_param",
        PREFIX + "_mdb": "_mdb",
    }
    found = {}
    with tarfile.open(fileobj=stream, mode="r|gz") as archive:
        for entry in archive:
            if entry.name in want:
                found[want[entry.name]] = archive.extractfile(entry).read()
                if len(found) == len(want):
                    break
    if set(found) != set(want.values()):
        raise ValueError("Incomplete primary archive")
    return found


def pack(files):
    hashes = {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
    cfg = files["_mdb"].decode().splitlines()
    check = float(
        next(row.split()[-1] for row in cfg if row.startswith("check_value "))
    )
    # Narrow decoder for this official one-dimensional primary 64-bit FITS image.
    fits = files["check_param"]
    assert fits[80:160].startswith(b"BITPIX") and b"-64 " in fits[80:160]
    end = next(i for i in range(0, len(fits), 80) if fits[i : i + 8] == b"END     ")
    offset = 2880 * ((end + 80 + 2879) // 2880)
    fid = np.frombuffer(fits, ">f8", count=7528, offset=offset).astype(float)
    # One Fortran sequential record; scipy's FortranFile requires an actual fd.
    binary = files["c_matrix_plik_v22.dat"]
    assert (
        int.from_bytes(binary[:4], "little")
        == 613 * 613 * 8
        == int.from_bytes(binary[-4:], "little")
    )
    low = np.frombuffer(binary, "<f8", offset=4, count=613 * 613).reshape(613, 613)
    assert np.all(np.triu(low, 1) == 0)
    np.savez_compressed(
        HERE / "plik_lite.npz",
        lower_triangle=low[np.tril_indices(613)],
        observation=np.loadtxt(io.BytesIO(files["cl_cmb_plik_v22.dat"])),
        bin_start=np.loadtxt(io.BytesIO(files["blmin.dat"])).astype(int)[:215],
        bin_end=np.loadtxt(io.BytesIO(files["blmax.dat"])).astype(int)[:215],
        weight=np.loadtxt(io.BytesIO(files["bweight.dat"]))[:2479],
        official_check_spectra=fid[:-1].reshape(3, 2509),
        official_check_calibration=fid[-1],
        official_check_loglike=check,
    )
    (HERE / "provenance.json").write_text(
        json.dumps({"archive": URL, "directory": PREFIX, "sha256": hashes}, indent=2)
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path)
    args = parser.parse_args()
    with (
        args.archive.open("rb") if args.archive else urlopen(URL, timeout=90)
    ) as stream:
        pack(extract(stream))
