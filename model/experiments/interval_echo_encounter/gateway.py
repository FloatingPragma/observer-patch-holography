"""Replay the registered gateway experiment, including enrollment refusals and quiet reads."""
import argparse
import hashlib
import json
from math import isclose, log

import numpy as np

from .gateway_kernel import kernel
from .gateway_source import HERE, forecast, require, settings
from model.candidate_laws.contrast_memory_spherical_v1.presence import load_json

RECEIPT = HERE / "gateway_counts.json"
RESULT = HERE / "gateway_result.json"
# Acquisition at 093c582cb. The digest covers every original sampling weight
# and count. Replay source fingerprints can change without moving this anchor.
ACQUISITION = "4ac3881772103bff00f0280958728d26116b0e50d0aeb1e307b221456bcaafc4"
MODEL = HERE.parents[1]
CONTRAST = "candidate_laws/contrast_memory_spherical_v1/"
SOURCES = [HERE / name for name in ("gateway_protocol.json", "gateway_kernel_protocol.json",
           "gateway_source.py", "gateway_kernel.py", "gateway.py")]
SOURCES += [MODEL / path for path in ("executable_model/carrier.py", "candidate_laws/interval_echo/candidate.py",
            *(CONTRAST + name for name in ("candidate.py", "law.py", "reader.py", "sectors.py", "presence.py",
                                           "common.py")))]
SOURCES += [MODEL / "executable_model" / "carriers" / f"FEDERATION_LEVEL{level}.npz" for level in range(3)]
SOURCES.append(MODEL / "requirements.txt")


def pins():
    result = {}
    for path in SOURCES:
        data = path.read_bytes()
        if path.suffix != ".npz":
            data = data.replace(b"\r\n", b"\n")
        result[path.relative_to(MODEL).as_posix()] = hashlib.sha256(data).hexdigest()
    return result


def registered_receipt(receipt):
    require(type(receipt) is dict and set(receipt) == {"sources", "weights", "counts"},
            "incomplete gateway receipt")
    require(receipt["sources"] == pins(), "measurement code or dependencies changed")
    payload = json.dumps([receipt["weights"], receipt["counts"]], sort_keys=True,
                         separators=(",", ":"), allow_nan=False).encode("ascii")
    require(hashlib.sha256(payload).hexdigest() == ACQUISITION, "registered acquisition differs")
    return receipt


def counts_law(p):
    # Freeze the finite sampling law to integer quanta. Its L1 discrepancy is
    # at most 31e-13. The quiet/refused split and every exact zero survive.
    scale = 10**13
    bins = np.rint(p * scale).astype(np.int64)
    bins[np.argmax(bins)] += scale - int(bins.sum())
    require(min(bins) >= 0 and np.all((p == 0) <= (bins == 0)), "invalid sampler law")
    return bins


def sample(weights, shots, rng):
    counts = []
    weight, left = int(sum(weights)), shots
    for item in weights[:-1]:
        item = int(item)
        draw = int(rng.binomial(left, item / weight)) if item and weight else 0
        left -= draw
        weight -= item
        counts.append(draw)
    return counts + [left]


def read(row, shots, p, L):
    p = np.asarray(p)
    require(p.shape == (32,) and np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1))
            and abs(p.sum() - 1) < 1e-13 and p[-1] < 1, "invalid joint probability law")
    require(type(row) is list and len(row) == 32 and all(type(x) is int and 0 <= x <= shots for x in row)
            and sum(row) == shots, "incomplete original-start receipt")
    data = np.array(row, dtype=np.int64)
    f = data / shots
    require(np.all(data[p == 0] == 0), "click in an impossible channel")
    require(np.all(abs(f - p) <= np.sqrt(2*p*(1-p)*L/shots) + 2*L/(3*shots) + 3.2e-12),
            "joint record law rejected")
    admitted = int(data[:-1].sum())
    require(admitted > 0, "no enrolled starts available for the conditional witness")
    hits = data[:30]
    f = np.r_[hits, hits.sum()] / admitted
    c = 2*L/admitted
    # Invert |f-p| <= sqrt(c*p) + c/3. Include the certified sampler
    # discrepancy after conditioning, for comparison with the unrounded law.
    upper = f + c/3 + c/2 + np.sqrt(c*(f+c/3) + c*c/4)
    rad = np.sqrt(c*upper) + c/3 + 6.4e-12/(1-p[-1])
    return np.maximum(0, f-rad), np.minimum(1, f+rad)


def summarize(cfg, laws, budget, counts):
    require(set(counts) == set(laws) and len(laws) == 81, "cohort menu mismatch")
    L = log(2*(32+31)*len(laws)/cfg["alpha"])
    intervals = {key: read(counts[key], cfg["shots"], p, L) for key,p in laws.items()}
    rows = {}
    for law in cfg["laws"]:
        rows[law] = {}
        for word in cfg["appointments"]:
            a,A = intervals[f"{law}_plus_{word}"]
            b,B = intervals[f"{law}_minus_{word}"]
            phase = np.maximum(0, np.maximum(a-B,b-A))
            row = {"cell_phase_interval": [float(a[-1]-B[-1]), float(A[-1]-b[-1])],
                   "largest_single_seam_phase_lower": float(max(phase[:-1]))}
            if law == "native":
                low = np.column_stack((np.sqrt(a),np.sqrt(b)))
                high = np.column_stack((np.sqrt(A),np.sqrt(B)))
                motion = np.maximum(0,low-2*budget)
                phase_band = np.maximum(0,np.maximum(low[:,0]-high[:,1], low[:,1]-high[:,0])-budget.sum(axis=1))
                row.update(cell_motion_lower=motion[-1].tolist(),
                           largest_seam_motion_lower=motion[:-1].max(axis=0).tolist(),
                           band_phase_cell_lower=float(phase_band[-1]),
                           largest_seam_band_phase_lower=float(max(phase_band[:-1])))
            rows[law][str(word)] = row
    return rows


def close(actual, saved, *, absolute=1e-9):
    if isinstance(actual, dict):
        require(type(saved) is dict and actual.keys() == saved.keys(), "result fields changed")
        for key in actual:
            close(actual[key], saved[key], absolute=absolute)
    elif isinstance(actual, (tuple,list)):
        require(type(saved) is list and len(actual) == len(saved), "result list changed")
        for left,right in zip(actual,saved,strict=True):
            close(left,right, absolute=absolute)
    elif type(actual) is float:
        require(type(saved) in (float,int) and isclose(saved, actual, abs_tol=absolute, rel_tol=1e-12),
                "numeric replay disagrees")
    else:
        require(saved == actual and type(saved) is type(actual), "result identity changed")


def save(path, obj):
    with path.open("x", encoding="utf8", newline="\n") as out:
        json.dump(obj, out, allow_nan=False, separators=(",",":"))
        out.write("\n")


def replay(record=False):
    cfg = settings()
    if record:
        require(not RECEIPT.exists() and not RESULT.exists(), "acquisition already exists")
    else:
        receipt = registered_receipt(load_json(RECEIPT))
    probes = kernel()
    laws, budget, atlas = forecast(cfg)
    rng = np.random.default_rng(cfg["seed"])
    if record:
        receipt = {"sources":pins(), "weights":{}, "counts":{}}
        for name,p in laws.items():
            weights = counts_law(p)
            receipt["weights"][name] = weights.tolist()
            receipt["counts"][name] = sample(weights, cfg["shots"], rng)
        save(RECEIPT, registered_receipt(receipt))
    else:
        require(set(receipt) == {"sources","weights","counts"} and set(receipt["weights"]) == set(laws),
                "sampling menu changed")
        for name,p in laws.items():
            weights = receipt["weights"][name]
            require(type(weights) is list and len(weights) == 32
                    and all(type(x) is int and 0 <= x <= 10**13 for x in weights) and sum(weights)==10**13
                    and np.all(np.asarray(weights)[p == 0] == 0)
                    and np.sum(abs(np.asarray(weights)/10**13-p))<3.2e-12, "archived sampling law changed")
            require(sample(weights,cfg["shots"],rng)==receipt["counts"].get(name), "seeded counts disagree")
    row = {"receiver":atlas, "amplitude_contamination":budget.tolist(), "kernel":probes,
           "bookings":summarize(cfg,laws,budget,receipt["counts"])}
    if record:
        save(RESULT,row)
    else:
        close(row, load_json(RESULT))
    print("Verified 81 cohorts, 32 outcomes per start and all gateway kernel probes.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    replay(parser.parse_args().record)
