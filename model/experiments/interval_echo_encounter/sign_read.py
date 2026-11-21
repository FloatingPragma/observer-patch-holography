"""Phase-sensitive sign of the free contrast from a single home-seam read."""
import json

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.law import seam_term
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector, SeamEvents
from model.executable_model.carrier import federation
from .transport import HERE, release, write_receipt

PROTOCOL = HERE / "sign_protocol.json"
RECEIPT = HERE / "sign_read.json"


def contacts(car, x, y):
    def local(p):
        return p // 12 == x // 12
    z = min(p for p in car.neighbours[y] if local(p) and p != x and p not in car.neighbours[x])
    w = min(p for p in car.neighbours[z] if local(p) and p not in (x, y, z))
    return z, w


def pulse_check(strength):
    # Isolated blank-collar seam and the sector pulse agree on every binary input
    # at the registered duration, including when the occupied rail is a spectator.
    u = expm(-0.5j * np.pi * seam_term(strength))
    indices = [0, 1, 3, 4]
    a, b = (1-1j)/2, (1+1j)/2
    reduced = np.array([[1, 0, 0, 0], [0, a, b, 0], [0, b, a, 0], [0, 0, 0, 1]])
    if not np.allclose(u[:, indices], np.eye(9)[:, indices] @ reduced, rtol=0, atol=2e-14):
        raise ArithmeticError("The isolated native seam differs from the applied half contact")


def campaign():
    spec = json.loads(PROTOCOL.read_text())
    sector = PairSector(federation(spec["carrier_level"]))
    internal, (x, y) = release(sector, spec["cell"])
    z, w = contacts(sector.car, x, y)
    screen = [sector.index(*s) for s in internal]
    address = sector.index(x, y)
    yz, zw = (sector.car.seams.index(tuple(sorted(e))) for e in [(y,z), (z,w)])
    arms = {}
    for strength in spec["strengths"]:
        pulse_check(strength)
        events = SeamEvents(sector, strength)
        def contact(v, seam):
            events.apply(v, seam, np.pi / 2)
        ready = events.start((x,y))
        contact(ready, yz)
        contact(ready, zw)
        paths = [sector.index(x, p) for p in (y, z, w)]
        expected = np.zeros_like(ready)
        expected[paths] = [(1-1j)/2, 0.5, 0.5j]
        if not np.allclose(ready, expected, rtol=0, atol=2e-15):
            raise ArithmeticError("Preparation differs from the registered three amplitudes")
        # A coherent column plus three individually propagated dephased branches.
        starts = np.zeros((sector.P, 4), complex)
        starts[:, 0] = ready
        starts[paths, [1, 2, 3]] = ready[paths]
        generator = -1j * sector.hamiltonian(strength)
        derivative = generator @ starts
        contact(derivative, yz)
        initial = starts.copy()
        contact(initial, yz)
        slope = 2 * (initial[address].conj() * derivative[address]).real
        readings = {}
        for t in spec["exposures"]:
            evolved = expm_multiply(t * generator, starts) if t else starts.copy()
            contact(evolved, yz)
            prob = abs(evolved)**2
            coherent, mixture = prob[:, 0], prob[:, 1:].sum(axis=1)
            records = []
            for weights in (coherent, mixture):
                if abs(weights.sum() - 1) > 1e-12:
                    raise ArithmeticError("Preparation or terminal branches lost probability")
                quiet = np.ones(sector.P, bool)
                quiet[screen] = False
                records.append(np.r_[weights[screen], weights[quiet].sum()].round(12).tolist())
            readings[str(t)] = records
        arms[str(strength)] = {"slopes": [round(float(slope[0]),12), round(float(sum(slope[1:])),12)],
                               "readings": readings}
    return {"ports": [x,y,z,w], "home_seams": [list(e) for e in internal],
            "order": ["coherent", "dephased"], "arms": arms}


def verify(result, archive):
    if any(result[key] != archive[key] for key in ("ports", "home_seams", "order")):
        raise ValueError("Sign-read source or outcome order differs")
    if result["arms"].keys() != archive["arms"].keys():
        raise ValueError("Sign-read arms differ")
    for strength, arm in result["arms"].items():
        stored = archive["arms"][strength]
        if arm["slopes"] != stored["slopes"] or arm["readings"].keys() != stored["readings"].keys():
            raise ValueError("Initial slopes or exposures differ")
        for t, rows in arm["readings"].items():
            target = np.asarray(stored["readings"][t])
            if target.shape != (2,31) or not np.allclose(rows, target, atol=2e-11, rtol=0):
                raise ValueError(f"Sign-read probabilities differ at {strength}, {t}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = campaign()
    if args.write:
        write_receipt(RECEIPT, result)
    else:
        verify(result, json.loads(RECEIPT.read_text()))
    print("Nominated-seam probabilities at exposure 1/16",
          {k: [a[0] for a in v["readings"]["0.0625"]] for k, v in result["arms"].items()})
