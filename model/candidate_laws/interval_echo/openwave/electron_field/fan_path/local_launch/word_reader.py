"""Read the fresh reference using only the receiving port and a single face."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from runpy import run_path

import numpy as np

from ..fixture import FanPath, archived_fixture
from .measure import HERE, pauli_record, run_case


def word_pairs(path, words):
    """Construct the exhaustive face bank and its predeclared disjoint matching."""
    used = set()
    pairs = []
    for face, (mask, clockwise) in enumerate(
        zip(path.field.masks, path.field.clockwise, strict=True)
    ):
        for a in words:
            b = a ^ mask
            if a >= b or a & mask not in (clockwise, clockwise ^ mask):
                continue
            selected = [
                w
                for w in words
                if w & ~mask == a & ~mask and w & mask in (clockwise, clockwise ^ mask)
            ]
            assert selected == [a, b]
            exclusive = not ({a, b} & used)
            if exclusive:
                used.update((a, b))
            pairs.append((face, words.index(a), words.index(b), exclusive))
    return pairs


def local_read(path, wave, written, admitted, wave_error):
    end = path.size - 1
    pad = (2 * np.sqrt(written) + wave_error) * wave_error + 1e-11
    reads = []
    port_totals = []
    selected = 0.0
    entangled = 0.0
    first = None
    for x, address in sorted(enumerate(path.cells[end]), key=lambda item: item[1]):
        block = path.blocks[end, x]
        words = block.words[end]
        # External configurations are identical for all words at this port.
        local = wave[block.rows].reshape(-1, len(words), wave.shape[1] // 2, 2)
        local = np.moveaxis(local, 2, 1)
        port_rate = np.vdot(local, local).real / 2
        port_totals.append(port_rate)
        matched_at_port = 0.0
        for face, i, j, exclusive in word_pairs(path, words):
            amplitudes = np.take(local, [i, j], axis=2).reshape(-1, 2, 2)
            rec = pauli_record(amplitudes, pad)
            probability = rec["reference_joint_success"]
            lo, hi = rec["partial_transpose_minimum_interval"]
            verdict = "entangled" if hi < 0 else "separable" if lo > 0 else "unresolved"
            if port_rate == 0:
                verdict = "zero support"
            label = [address, face, words[i], words[j]]
            if exclusive:
                matched_at_port += probability
                if verdict == "entangled":
                    entangled += probability
                    if first is None:
                        first = label
            reads.append(
                {
                    "port_face_and_words": label,
                    "in_exclusive_bank": exclusive,
                    "joint_record": verdict,
                    **rec,
                }
            )
        assert matched_at_port <= port_rate + 2e-14
        selected += matched_at_port
    cell = float(sum(port_totals))
    np.testing.assert_allclose(
        cell, np.linalg.norm(wave[path.owner == end]) ** 2 / 2, atol=1e-14
    )
    outcomes = [
        1 - admitted,
        admitted - written,
        written - cell,
        cell - selected,
        selected - entangled,
        entangled,
    ]
    assert min(outcomes) >= -2e-14
    np.testing.assert_allclose(sum(outcomes), 1, rtol=0, atol=1e-14)
    return {
        "observable_absolute_error_upper_per_pair": float(pad),
        "exclusive_bank_outcomes": np.maximum(outcomes, 0).tolist(),
        "first_entangled_exclusive_label": first,
        "pair_readings": reads,
    }


def run():
    with archived_fixture() as archived:
        cases = [
            run_case(
                FanPath(archived, n, 1, True),
                FanPath(archived, n, -1, True),
                local_read,
            )
            for n in (3, 4)
        ]
    return {
        "receiver_design_sha256": sha256(
            (HERE / "RECEIVER.md").read_bytes()
        ).hexdigest(),
        "original_design_sha256": sha256((HERE / "DESIGN.md").read_bytes()).hexdigest(),
        "exclusive_bank_outcome_order": [
            "original_clock_refused",
            "source_write_refused",
            "charge_outside_receiver",
            "charge_in_receiver_unmatched_word",
            "matched_label_without_certified_entanglement",
            "matched_label_with_certified_entanglement",
        ],
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run()
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(HERE.parents[1] / "replay.py"))["compare"]
        compare(
            json.loads(json.dumps(record, allow_nan=False)),
            json.loads((HERE / "word_reads.json").read_text()),
            "word_reader",
        )
        print(
            "Reproduced every local signed material and triangle read, with both rates and source arms."
        )


if __name__ == "__main__":
    main()
