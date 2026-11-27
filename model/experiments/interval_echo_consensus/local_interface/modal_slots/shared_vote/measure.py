"""Measure two full modal labels after an edit on common or distinct packets."""

from itertools import product
from pathlib import Path
import json

import numpy as np
from scipy.linalg import expm
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .. import measure as modal

HERE = Path(__file__).resolve().parent
SHOTS, SEED = 2**18, 320177
LAWS = HERE / "acquisition_laws.json"
# expm on the 9x9 seam terms reproduces the rail probabilities to about
# 1e-14; a larger difference is a changed law rather than roundoff.
LAW_TOLERANCE = 1e-12


def archived_laws():
    return json.loads(LAWS.read_text())["branch_probabilities"]


def label_pair(labels, observer):
    for hi, high in enumerate(labels):
        for lo, low in enumerate(labels):
            changed = [i for i in range(5) if high[i] != low[i]]
            if len(changed) != 1:
                continue
            field = changed[0]
            if observer in (86, 88) and field != 2:
                continue
            if low[field] != 0 or high[field].bit_count() != 1:
                continue
            (bit,) = np.flatnonzero(
                np.array(modal.packet_bits(low)) != modal.packet_bits(high)
            )
            return hi, lo, field, int(bit)
    raise ValueError("no two literal labels obey the frozen choice")


def physical_packets(programs, choices):
    packets = {}
    for o, prog in programs.items():
        encoded = modal.encode(prog, choices[o][0])
        for node, rails in zip(prog["supports"], encoded, strict=True):
            if node in packets:
                np.testing.assert_array_equal(packets[node], rails)
            else:
                packets[node] = rails.copy()
    return packets


def modal_labels(packets, programs):
    return {
        o: modal.mode(modal.from_ports([packets[node] for node in prog["supports"]]))
        for o, prog in programs.items()
    }


def statistics(progs, choices):
    # A table lookup of the measured literal on each addressed code column.
    reads = {
        o: [modal.read(prog, i) for i in range(len(prog["labels"]))]
        for o, prog in progs.items()
    }
    assert all(list(progs[o]["labels"]) == seq for o, seq in reads.items())
    triples = set(product(*(reads[o] for o in modal.IDS)))
    return {
        "sharp_triples": len(triples),
        "minimum_private_edits": {
            str(o): len(progs[o]["private"]) - 48 for o in choices
        },
    }


def outcome_table(progs, choices, pair, deficit, board, mapping):
    packets = physical_packets(progs, choices)
    for o in pair:
        _, _, _, bit = choices[o]
        # These distinct, full-pi flips act on a private physical rail pair.
        private = progs[o]["private"]
        take = 0 if deficit is None else len(private) - 48 - deficit
        for node in private[:take]:
            assert tuple(packets[node][bit]) == (0, 2)
            packets[node][bit] = [2, 0]
    primary, partner = pair
    shared = {o: set(progs[o]["supports"]) - set(progs[o]["private"]) for o in pair}
    if shared[primary] & shared[partner]:
        node = min(shared[primary] & shared[partner])
        assert choices[primary][3] == choices[partner][3]
        probes = [(node, choices[primary][3])]
    else:
        probes = [(min(shared[o]), choices[o][3]) for o in pair]
    for node, bit in probes:
        first = 12 * mapping[node] + 2 * bit
        assert first + 1 in board.neighbours[first] and tuple(packets[node][bit]) == (
            2,
            0,
        )
    branches = []
    for flags in product((0, 1), repeat=len(probes)):
        for (node, bit), flag in zip(probes, flags, strict=True):
            packets[node][bit] = [2 - 2 * flag, 2 * flag]
        labels = modal_labels(packets, progs)
        indices = {o: progs[o]["labels"].index(labels[o]) for o in progs}
        for o in set(progs) - set(pair):
            assert indices[o] == choices[o][0]
        outputs = tuple(int(indices[o] == choices[o][0]) for o in pair)
        assert all(indices[o] in choices[o][:2] for o in pair)
        branches.append(
            {
                "rail_bits": list(flags),
                "label_indices": [indices[o] for o in modal.IDS],
                "joint_bin": 2 * outputs[0] + outputs[1],
            }
        )
    return probes, branches


def acquire(archived=None):
    """Return the counts and the branch laws drawn from.

    A replay passes the archived laws. Fresh probabilities must agree with
    them within LAW_TOLERANCE, and the seeded draws use the archived values.
    """
    progs, _, mapping = modal.setup()
    choices = {o: label_pair(p["labels"], o) for o, p in progs.items()}
    summary = statistics(progs, choices)
    board = federation(2)
    contrast = candidate.KAPPA * law._load_ops()[1]
    native = law.seam_term(candidate.KAPPA)
    addresses = {port for program in progs.values() for port in program["address"]}
    assert all(
        12 * owner + j not in addresses for owner in mapping.values() for j in range(12)
    )
    seed = np.random.default_rng(SEED)
    banks, laws = [], []
    for pair in [(86, 88), (86, 247)]:
        for name, deficit in [("code", None), ("one_short", 1), ("boundary", 0)]:
            probes, branches = outcome_table(
                progs, choices, pair, deficit, board, mapping
            )
            for pulse in (0, 1, 2, 3):
                for arm, h in [("native", native), ("parked", contrast)]:
                    amp = expm(-1j * h * (pulse * np.pi / 3))[2, 6]  # 20 -> 02
                    p = float(np.clip(abs(amp) ** 2, 0, 1))
                    probabilities = np.array(
                        [
                            np.prod([p if bit else 1 - p for bit in b["rail_bits"]])
                            for b in branches
                        ]
                    )
                    branch_law = probabilities / probabilities.sum()
                    if archived is not None:
                        stored = np.asarray(archived[len(banks)], dtype=float)
                        if stored.shape != branch_law.shape or not np.all(
                            abs(branch_law - stored) <= LAW_TOLERANCE
                        ):
                            raise ArithmeticError(
                                "branch probabilities differ from the acquisition law"
                            )
                        branch_law = stored
                    laws.append(branch_law.tolist())
                    counts = seed.multinomial(SHOTS, branch_law)
                    joint = [0] * 4
                    for branch, count in zip(branches, counts, strict=True):
                        joint[branch["joint_bin"]] += int(count)
                    banks.append(
                        {
                            "observers": list(pair),
                            "private_edits": name,
                            "pulse_pi_over_three": pulse,
                            "flight": arm,
                            "probe_nodes_bits": probes,
                            "branches": branches,
                            "rail_counts": counts.tolist(),
                            "joint_counts": joint,
                        }
                    )
    data = {
        "starts_per_setting": SHOTS,
        "seed": SEED,
        "labels": choices,
        "summary": summary,
        "banks": banks,
    }
    return data, laws


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data, laws = acquire(None if args.write else archived_laws())
    text = json.dumps(data, indent=1, allow_nan=False) + "\n"
    if args.write:
        (HERE / "counts.json").write_text(text)
        LAWS.write_text(
            json.dumps({"branch_probabilities": laws}, indent=1, allow_nan=False)
            + "\n"
        )
    else:
        assert text == (HERE / "counts.json").read_text()
    print(data["summary"])
    for b in data["banks"]:
        if (
            b["private_edits"] == "boundary"
            and b["pulse_pi_over_three"] == 1
            and b["flight"] == "native"
        ):
            print(b["observers"], b["joint_counts"])
