"""Exact local guide dynamics, original charge-bank transfer and fresh remote calibration."""

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.contrast_memory_spherical_v1 import law
from ..native_source_probe import circuit
from .geometry import fixtures

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "native_source_probe/regional_census/results.json"
SHOTS = 65536
BITS = tuple(product((0, 1), repeat=3))
TARGETS = {"vacancy": (0,), "host": (1,), "excess": (2,), "internal_pair": (2, 0)}
CONTROLS = (
    "native",
    "cut_bridge",
    "half_bridges",
    "exchange_only_flight",
    "exchange_only_seams",
)


@lru_cache(maxsize=None)
def refocused_gate(word, kappa):
    """Evaluate the full three-port word with a changed blank contact and fixed private clock."""
    state = np.eye(27, dtype=complex).reshape(3, 3, 3, 27)
    for kind, a, b, ratio in word:
        dwell = np.pi * float(ratio)
        if kind == "reference":
            axes = [a, *[j for j in range(3) if j != a], 3]
            value = state.transpose(axes)
            phase = np.exp(-1j * 8 * dwell * np.arange(3) ** 2)
            updated = value * phase[:, None, None, None]
        else:
            axes = [a, b, 3 - a - b, 3]
            value = state.transpose(axes).reshape(9, -1)
            gate = expm(-1j * dwell * law.seam_term(kappa))
            updated = (gate @ value).reshape(3, 3, 3, 27)
        state = updated.transpose(np.argsort(axes))
    return state.reshape(27, 27)


@lru_cache(maxsize=None)
def write_blocks(setting, kappa=8):
    ratio = Fraction(2 if setting == "parity" else 1, 65)
    word = circuit.pointer.diagonal_word(0, 2, ratio)
    gate = refocused_gate(tuple(word), kappa)
    if kappa == 8:
        np.testing.assert_allclose(
            gate, circuit.apply(word, np.eye(27)), atol=2e-13, rtol=0
        )
    blocks = []
    for ids in circuit.CODE:
        block = gate[np.ix_(ids, ids)]
        np.testing.assert_allclose(
            block.conj().T @ block, np.eye(2), atol=2e-13, rtol=0
        )
        blocks.append(block)
    return blocks


def unitary_layers(control):
    kappa = (
        0
        if control in ("exchange_only_flight", "exchange_only_seams")
        else candidate.KAPPA
    )
    h = law.seam_term(kappa)
    ix = [law.SRC.index((2, 0)), law.SRC.index((0, 2))]
    local = h[np.ix_(ix, ix)]
    active = expm(-1j * np.pi / 2 * local)
    hop = np.kron(np.kron(active, active), active)
    # Compression onto the one-double-per-dimer guards cancels exchange
    # across different dimers. Keep all diagonal terms, including vacuum.
    links = []
    for bits in BITS:
        energy = 0
        for i in (0, 1):
            if i == 0 and control == "cut_bridge":
                continue
            source = law.SRC.index((2 * bits[i], 2 * (1 - bits[i + 1])))
            energy += h[source, source]
        links.append(energy)
    dwell = np.pi / 65 * (0.5 if control == "half_bridges" else 1.0)
    return hop @ np.diag(np.exp(-1j * dwell * np.array(links)))


@lru_cache(maxsize=None)
def starts():
    local = law.seam_term(candidate.KAPPA)
    ix = [law.SRC.index((2, 0)), law.SRC.index((0, 2))]
    active = expm(-1j * np.pi / 2 * local[np.ix_(ix, ix)])
    # Each separately prepared dimer is in its lower free mode.
    return active, np.ones(2, complex) / np.sqrt(2)


def calibration(targets, setting, control):
    head_hop, plus = starts()
    head = plus.copy()
    for n in targets:
        head = (
            write_blocks(setting, 0 if control == "exchange_only_seams" else 8)[n]
            @ head
        )
    if setting == "parity":
        # Private quarter phase on the head's first port when its load is two.
        head[0] *= np.exp(-1j * 4 * candidate.KAPPA * np.pi / 64)
    state = np.kron(np.kron(head_hop @ head, plus), plus)
    tick = unitary_layers(control)
    rows = []
    for t in range(3):
        probs = np.sum(np.abs(state.reshape(4, 2)) ** 2, axis=0)
        np.testing.assert_allclose(probs.sum(), 1, atol=2e-13, rtol=0)
        rows.append(probs)
        state = tick @ state
    return np.array(rows)


def band_transfer(census):
    surveys = []
    decoder = census["pointer_decoder"]
    for launch in census["launches"]:
        panels = {}
        for mode, cohort in launch["panels"].items():
            bank = []
            for i, entry in enumerate(cohort):
                signed = [0, 0, 0]
                for code, count in entry["counts"]:
                    if code == census["total_output_bins"] - 1:
                        continue
                    if (
                        code >= 73 * len(decoder)
                        or decoder[code % len(decoder)] is None
                    ):
                        raise ValueError("Unexpected off-code source outcome")
                    signed[1 + decoder[code % len(decoder)]] += count
                if (
                    sum(signed) != entry["accepted"]
                    or signed[0] + signed[2] != entry["charged"]
                ):
                    raise ValueError("Source census does not reconcile")
                bank.append(
                    dict(
                        profile=i,
                        signed_source_counts=signed,
                        initial_filter_refusals=census["shots"] - sum(signed),
                        parity_dark_probability_interval=entry["simultaneous_interval"],
                    )
                )
            panels[mode] = bank
        surveys.append(
            dict(
                launch=launch["owner"],
                native_per_event_error=launch["native_per_event_error"],
                panels=panels,
            )
        )
    return surveys


def run():
    source = SOURCE.read_bytes()
    census = json.loads(source)
    rng = np.random.default_rng(190325)
    panels = []
    for control, (name, targets), setting in product(
        CONTROLS, TARGETS.items(), ("parity", "sign")
    ):
        probs = calibration(targets, setting, control)
        # Snap roundoff at exact probability endpoints before generating counts.
        probs = np.round(probs, 13)
        outcomes = []
        for t, (p0, _) in enumerate(probs):
            zero = int(rng.binomial(SHOTS, p0))
            outcomes.append(
                dict(
                    tick=t,
                    probability_word20=float(p0),
                    terminal_counts=[zero, SHOTS - zero],
                )
            )
        panels.append(
            dict(guide=control, source=name, quadrature=setting, observations=outcomes)
        )
    return dict(
        source_sha256=sha256(source).hexdigest(),
        calibration_shots_per_setting=SHOTS,
        initial_pair_shots_per_setting=census["shots"],
        calibrations=panels,
        fixtures=fixtures(census["profiles"]),
        equivalent_source_banks=band_transfer(census),
        entire_band_remote_odd_limit=census["universal_distinct_owner_bound"],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(
            json.dumps(result, separators=(",", ":"), sort_keys=True) + "\n"
        )
    elif json.loads(path.read_text()) != result:
        raise ValueError("Remote detector receipt differs")
    for s in result["calibrations"]:
        if s["source"] in ("excess", "host") and s["quadrature"] == "parity":
            print(
                s["guide"],
                s["source"],
                [r["terminal_counts"][0] for r in s["observations"]],
            )
    print(
        "Universal single-remote-read upper limit:",
        result["entire_band_remote_odd_limit"]["upper"],
    )


if __name__ == "__main__":
    main()
