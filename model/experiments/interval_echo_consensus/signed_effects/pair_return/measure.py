"""Prepare and reread the coordinate and uniform states in all six address pairs."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .. import fixed_pointer
from ..full_address import measure as address
from ..phase_table import acquire as apparatus
from ..phase_table import measure as pulses
from ..sequential import compare, require

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "results.json"
DIM = 12
LABELS = 17
SHOTS = 8192
SEED = 1310612


def read_branches(calls, prior):
    """Each column comes from the pointer and restored data on a 240-port word."""
    result = np.zeros((LABELS, DIM, DIM), dtype=complex)
    for column in range(DIM):
        word, phase = address.execute_basis(calls, 1, column, prior)
        raw = address.read_word(word)
        positions = np.flatnonzero(word[:DIM] != 1)
        require(len(positions) == 1 and word[positions[0]] == 2, "data left its code")
        row = int(positions[0])
        expected = address.initial_word(1, row, raw)
        require(raw < 16 and np.array_equal(word, expected) and word.sum() == 240,
                "incorrect return of a pointer, host, load or reference")
        result[raw, row, column] = phase
    return result


def pair_effect(branches, start):
    return sum(k.conj().T @ k for k in branches[start:start + 2])


def probe_words(start, name):
    """Place the probe in this pair and bring its output to the terminal analyzer."""
    route = [] if start == 0 else (
        pulses.transposition(0, start) + pulses.transposition(1, start + 1)
    )
    prep = apparatus.preparation_word(name) + route
    analyzer = fixed_pointer.inverse(route)
    if name == "plus":
        analyzer += apparatus.x_analyzer_word()
    return prep, analyzer


def read_probe(first, terminal, start, name):
    prep, analyzer = probe_words(start, name)
    state = pulses.source_unitary(tuple(prep), 1)[:, 0]
    if name == "plus":
        intended = (np.eye(DIM)[start] + np.eye(DIM)[start + 1]) / np.sqrt(2)
    else:
        intended = np.eye(DIM)[start + (name == "port1")]
    require(abs(abs(np.vdot(intended, state)) - 1) < 2e-12, "incorrect material probe")
    rotation = pulses.source_unitary(tuple(analyzer), 1)
    joint = np.zeros((LABELS, LABELS))
    for raw, branch in enumerate(first):
        if np.any(branch):
            # The second native address XOR starts from the measured pointer.
            paths = terminal[raw] @ (rotation @ branch @ state)
            joint[raw] = np.sum(abs(paths)**2, axis=1)
    require(abs(joint.sum() - 1) < 3e-12, "incomplete two-read law")
    target = int(name == "port1")
    returned = sum(joint[r, r ^ target] for r in (start, start + 1))
    return joint, float(returned)


def report():
    local = address.local_check()
    calls = address.schedule()
    coarse_calls = [call for call in calls if call["bit"] != 0]
    heads = {"pair": read_branches(coarse_calls, 0), "address_then_forget": read_branches(calls, 0)}
    terminal = {raw: read_branches(calls, raw) for raw in range(DIM)}
    rng = np.random.Generator(np.random.PCG64(SEED))
    banks = []
    for start in range(0, DIM, 2):
        expected = np.diag(np.isin(np.arange(DIM), [start, start + 1]))
        for kind, head in heads.items():
            require(np.allclose(pair_effect(head, start), expected, atol=2e-13, rtol=0),
                    "incorrect grouped effect")
            for name in ("port0", "port1", "plus"):
                joint, returned = read_probe(head, terminal, start, name)
                prediction = 0.5 if name == "plus" and kind == "address_then_forget" else 1
                require(abs(returned - prediction) < 2e-12, "incorrect return probability")
                counts = rng.multinomial(SHOTS, joint.ravel() / joint.sum()).reshape(joint.shape)
                target = int(name == "port1")
                hits = sum(int(counts[r, r ^ target]) for r in (start, start + 1))
                require(abs(hits / SHOTS - prediction) <= np.sqrt(np.log(2 * 36 / .01) / (2 * SHOTS)),
                        "return frequency outside simultaneous 99 percent interval")
                banks.append({"pair_start": start, "instrument": kind, "probe": name,
                              "return_probability": round(returned, 12), "returned": hits,
                              "nonzero_joint_counts": [[int(i), int(j), int(counts[i, j])]
                                                       for i, j in zip(*np.nonzero(counts))]})
    files = [Path(__file__), Path(apparatus.__file__), Path(pulses.__file__),
             address.ROOT / "Lean/EventAlgebra/FiniteProbeInstrumentSelection.lean"]
    hashes = address.sources()
    for file in files:
        hashes[str(file.relative_to(address.ROOT))] = hashlib.sha256(file.read_bytes()).hexdigest()
    return {"source_sha256": hashes, "native_checks": local, "initial_material_load": 240,
            "excess_data_ports": list(range(DIM)), "balancing_vacancy": address.REFERENCE,
            "public_group": "floor(first_pointer/2) for valid words; 16 is a malformed raw pointer",
            "joint_raw_labels": list(range(LABELS)), "zero_bins_omitted": True,
            "terminal_address": "first_pointer XOR second_pointer on valid words",
            "head_copy_calls": {"pair": len(coarse_calls), "address_then_forget": len(calls)},
            "terminal_copy_calls": len(calls), "shots_per_setting": SHOTS, "seed_pcg64": SEED,
            "original_starts": SHOTS * len(banks), "banks": banks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    actual = report()
    if args.write:
        with RECEIPT.open("x") as output:
            json.dump(actual, output, indent=2, allow_nan=False)
            output.write("\n")
    else:
        compare(json.loads(RECEIPT.read_text()), actual)
    print(f"Verified {len(actual['banks'])} return tests on six rank-two outcomes.")


if __name__ == "__main__":
    main()
