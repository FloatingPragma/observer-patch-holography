"""Compare the record-generated walks with every deposited cell-seven terminal record."""
import base64
import hashlib
import json
import zlib
from pathlib import Path

import numpy as np
from scipy.special import xlogy

from model.executable_model import carrier
from .record_release import pair, source


DEPOSIT = Path(__file__).resolve().parents[3] / "openwave/interval_echo/row_04/registered_release.json"
DIGEST = "1b5a1262baa76e6774f90f612033bd4a58bc5a6d90aefd76f999c83779311377"
ENROLLMENT = "bb697927ff37446130d4c6e8c68126b6a3eac00d"
TIMES = (3, 11, 29)


def deviance(counts, p):
    expected = pair.SHOTS*p
    return 2*np.sum(xlogy(counts, counts)-xlogy(counts, expected), axis=-1)


def read_deposit():
    raw = DEPOSIT.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DIGEST:
        raise ValueError("The deposited observation stream changed")
    return json.loads(raw)


def compare(arm):
    archived = read_deposit()
    edges, index, _, _ = pair.screen(carrier.HERE / "carriers/FEDERATION_LEVEL0.npz")
    key = hashlib.sha256(b"M10.3:row4:" + bytes.fromhex(ENROLLMENT)).digest()
    eligible = [v for v in range(20) if v not in (0, 13)]
    home = eligible[int.from_bytes(key[:4], "big") % len(eligible)]
    patch = edges[np.all(edges//12 == home, axis=1)]
    gateways = np.unique(edges[edges[:, 0]//12 != edges[:, 1]//12])
    candidates = patch[np.isin(patch, gateways).any(axis=1)]
    launch = candidates[int.from_bytes(key[4:8], "big") % len(candidates)]
    seed = int.from_bytes(key[8:16], "big")
    if (archived["registration_commit"] != ENROLLMENT or archived["cell"] != home
            or archived["launch"] != launch.tolist() or archived["outcomes"] != patch.tolist()+["quiet"]
            or archived["shots"] != pair.SHOTS or arm not in ("primitive", "assembled", "off")):
        raise ValueError("Unexpected registered experiment")
    h, _, _, _ = source(arm)
    states = pair.evolve(h, index[tuple(launch)], TIMES)
    arm_index = int(arm == "off")
    rows = []
    for j, (t, state, reading) in enumerate(zip(TIMES, states, archived["arms"][arm_index]["records"], strict=True)):
        if reading["exposure"] != t:
            raise ValueError("Wrong appointment")
        outcomes = zlib.decompress(base64.b64decode(reading["draws_zlib_base64"], validate=True))
        observed = np.frombuffer(outcomes, np.uint8)
        counts = np.bincount(observed, minlength=len(patch)+1)
        if len(observed) != pair.SHOTS or counts.tolist() != reading["counts"]:
            raise ValueError("Incomplete observations")
        p = pair.local_law(state, index[patch[:, 0], patch[:, 1]])
        p /= p.sum()
        if not np.isfinite(p).all() or np.any(p <= 0):
            raise ArithmeticError("The registered detector has a null or invalid predicted bin")
        D = float(deviance(counts, p))
        simulations = np.random.default_rng([45731, j]).multinomial(pair.SHOTS, p, size=9999)
        greater = int(np.count_nonzero(deviance(simulations,p) >= D))
        if arm != "assembled":
            rng = np.random.default_rng([seed, arm_index, j])
            if rng.choice(len(p), size=pair.SHOTS, p=p).astype(np.uint8).tobytes() != outcomes:
                raise ValueError("The same-stream control changed a deposited record")
        rows.append(dict(exposure=t, probabilities=p.tolist(), deviance=D, upper_rank=greater+1))
    return dict(arm=arm, launch=launch.tolist(), home=home, empirical_stream_sha256=DIGEST,
                settings=rows, family_rejected=any(r["upper_rank"]/10000 <= .01/3 for r in rows))


if __name__ == "__main__":
    print(json.dumps([compare(arm) for arm in ("primitive", "assembled", "off")], allow_nan=False))
