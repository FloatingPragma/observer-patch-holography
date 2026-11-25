"""Independently read the public stream using the exact two-rail instrument contract."""

import gzip
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
PHASE = ROOT / "code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json"
CENSUS = ROOT / "code/born_context_phase_lift/BORN_CONTEXT_WEB_PAYLOAD.v1.json"


# The verifier is the independent checkout anchor for this fixed run. The
# native acquisition separately reconstructs every branch's pulse-word hash.
ARCHIVED_RUN = "87be13ecb9c40295e593dd3ceb28bfd5499888182ee1a280084dec77a0d50afb"


def check(condition, message):
    if not condition:
        raise ValueError(message)


def exact_replay(shots, seed):
    # This reader never imports the acquisition, candidate, route or branch compiler.
    tables = json.loads(PHASE.read_bytes())["contexts"]
    counts = json.loads(CENSUS.read_bytes())["context_web"]["diagonal_context"]["realized_outcome_counts"]
    rng = np.random.Generator(np.random.PCG64(seed))
    result = np.full((shots, 6), 255, dtype=np.uint8)

    def fair_tickets(size, width):
        weights = 2 ** np.arange(width - 1, -1, -1)
        return weights @ (rng.random((width, size)) < 0.5)

    tickets = fair_tickets(shots, 8)
    result[:, 0] = tickets
    keep = np.flatnonzero(tickets < sum(counts))
    rail = (tickets[keep] >= counts[0]).astype(int)
    result[keep, 1] = (rng.random(len(keep)) < rail)
    context = fair_tickets(len(keep), 3)
    result[keep, 2] = context
    x_basis = np.array([[1, 1], [1, -1]]) / np.sqrt(2)

    def number(terms):
        a, b, c, d = map(Fraction, terms)
        return float(a) + float(b) * np.sqrt(3) + 1j * (float(c) + float(d) * np.sqrt(3))

    for slot, entry in enumerate(tables):
        yes = np.array([[number(value) for value in row] for row in entry["effect"]])
        for bit in (0, 1):
            rows = keep[(context == slot) & (rail == bit)]
            joint = np.zeros((17, 17))
            for label, branch in ((0, yes), (2, np.eye(2) - yes)):
                joint[label, :2] = abs(x_basis @ branch[:, bit])**2
            cumulative = np.cumsum(joint.flatten() / joint.sum())
            sampled = np.searchsorted(cumulative, rng.random(len(rows)), side="right")
            result[rows, 3], result[rows, 4] = np.divmod(sampled, 17)
    result[keep, 5] = rng.random(len(keep)) < rail
    return result


def verify(stream, receipt):
    descriptor = {key: value for key, value in receipt.items() if key != "run_sha256"}
    bound = json.dumps(descriptor, separators=(",", ":"), sort_keys=True).encode()
    check(hashlib.sha256(bound).hexdigest() == receipt["run_sha256"], "run descriptor differs")
    required = {str(PHASE.relative_to(ROOT)), str(CENSUS.relative_to(ROOT)),
                str((HERE / "acquire.py").relative_to(ROOT)),
                "model/candidate_laws/contrast_memory_spherical_v1/law.py",
                "model/candidate_laws/interval_echo/candidate.py",
                "model/executable_model/carriers/FEDERATION_LEVEL0.npz"}
    check(required <= receipt["source_sha256"].keys(), "missing source or calibration")
    check(list(receipt["word_sha256_by_context_including_all_raw_branches"]) == receipt["contexts_in_bit_order"],
          "incomplete operation bindings")
    for name, digest in receipt["source_sha256"].items():
        check(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, f"source differs: {name}")
    check((receipt["launches"], receipt["seed_pcg64"]) == (65536, 17911168), "wrong acquisition")
    expected = exact_replay(receipt["launches"], receipt["seed_pcg64"])
    check(expected.tobytes() == stream, "trace differs from independent source/operation replay")
    check(hashlib.sha256(stream).hexdigest() == receipt["uncompressed_trace_sha256"], "stream digest differs")
    refused = expected[:, 1] == 255
    check(int(np.count_nonzero(refused)) == receipt["refused_launches"], "refusal count differs")
    check(receipt["run_sha256"] == ARCHIVED_RUN, "descriptor differs from archived native run")
    return len(expected)


if __name__ == "__main__":
    public = gzip.decompress((HERE / "records.bin.gz").read_bytes())
    manifest = json.loads((HERE / "results.json").read_bytes())
    print(f"{verify(public, manifest)} complete trials verified against the archived projector table")
