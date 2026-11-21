"""Terminal home-seam reads and the complete position census on fresh releases."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model.carrier import federation

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "confinement_protocol.json"
RECEIPT = HERE / "transport.json"

# Hash the LF representation of the registered files on every platform.
PROTOCOL_DIGESTS = {
    "confinement_protocol.json": "8cbb398148793923b726daf2a63d49678cefbc3fff4152a63133c4cb433f7a46",
    "sign_protocol.json": "8387837b36e99512a77e1e3405358a000690789d550a595ed3b5a3ee5264e218",
    "conjugation_protocol.json": "2294f3b91d95c3054118cbe49de8d480c2ba75c09ed7c47d1c94be5bf7991b57",
}


def read_protocol(path):
    raw = path.read_bytes().replace(b"\r\n", b"\n")
    if hashlib.sha256(raw).hexdigest() != PROTOCOL_DIGESTS[path.name]:
        raise ValueError(f"Registered protocol differs: {path.name}")
    return json.loads(raw)


def probability_vector(values, size):
    """Validate one complete 12-decimal published distribution before comparison."""
    if len(values) != size or any(type(p) not in (int, float) for p in values):
        raise ValueError("Expected a complete vector of numeric probabilities")
    vector = np.asarray(values, dtype=float)
    if not np.isfinite(vector).all() or np.any(vector < 0) or np.any(vector > 1):
        raise ValueError("Probability outside the finite unit interval")
    if abs(vector.sum() - 1) > size * (0.5e-12 + np.finfo(float).eps):
        raise ValueError("Probability vector does not sum to one at receipt precision")
    return vector


def write_receipt(path, result):
    """Keep each numeric spectrum or outcome vector together on one line."""
    content = json.dumps(result, indent=2, allow_nan=False)
    content = re.sub(r"\[\s*([-+0-9.eE,\s]+?)\s*\]",
                     lambda m: "[" + " ".join(m[1].split()) + "]", content)
    path.write_text(content + "\n", encoding="utf-8")


def release(sector, home):
    internal = [e for e in sector.car.seams if all(p // 12 == home for p in e)]
    gateways = {p for seam in sector.car.external for p in seam}
    source = next(seam for seam in internal if gateways.intersection(seam))
    return internal, source


def census(sector, home, internal, state):
    weights = abs(state) ** 2
    home_pair = (sector.aa // 12 == home) & (sector.bb // 12 == home)
    internal_pair = sector.aa // 12 == sector.bb // 12
    seam = sector.adjacent
    bins = [weights[sector.index(*e)] for e in internal]
    bins.extend(weights[mask].sum() for mask in (
        home_pair & ~seam, ~internal_pair & seam, internal_pair & ~home_pair & seam,
        ~home_pair & ~seam))
    if not abs(sum(bins) - 1) <= 1e-10:
        raise ArithmeticError("Terminal census lost probability")
    return [round(float(p), 12) for p in bins]


def release_metadata(sector, internal, source):
    return {"dimension": sector.P, "source": list(source), "home_seams": [list(e) for e in internal],
            "complement": ["home nonadjacent", "external seam", "other internal seam",
                           "nonadjacent outside home window"]}


def campaign():
    protocol = read_protocol(PROTOCOL)
    settings = protocol["release"]
    sector = PairSector(federation(settings["level"]))
    internal, source = release(sector, settings["cell"])
    if (KAPPA, source) != (settings["strengths"][0], (156, 157)):
        raise ValueError("Law or launch differs from the registered release")
    start = np.zeros(sector.P, dtype=complex)
    start[sector.index(*source)] = 1
    arms = {}
    for strength in settings["strengths"]:
        generator = -1j * sector.hamiltonian(strength)
        state, previous = start, 0
        readings = {}
        for exposure in settings["exposures"]:
            state = expm_multiply((exposure - previous) * generator, state)
            readings[str(exposure)] = census(sector, settings["cell"], internal, state)
            previous = exposure
        arms[str(strength)] = readings
    return {**release_metadata(sector, internal, source), "arms": arms}


def verify(result, archive):
    for key in ("dimension", "source", "home_seams", "complement"):
        if result[key] != archive[key]:
            raise ValueError(f"Receipt mismatch at {key}")
    if result["arms"].keys() != archive["arms"].keys():
        raise ValueError("Interaction arms differ")
    for strength, readings in result["arms"].items():
        retained = archive["arms"][strength]
        if readings.keys() != retained.keys():
            raise ValueError("Exposure menus differ")
        for exposure, bins in readings.items():
            other = probability_vector(retained[exposure], 34)
            if np.abs(bins - other).sum() > 1e-9:
                raise ValueError(f"Census mismatch at {strength}, {exposure}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = campaign()
    if args.write:
        write_receipt(RECEIPT, result)
    else:
        verify(result, json.loads(RECEIPT.read_text()))
    print("Verified", sum(map(len, result["arms"].values())), "complete terminal censuses")
