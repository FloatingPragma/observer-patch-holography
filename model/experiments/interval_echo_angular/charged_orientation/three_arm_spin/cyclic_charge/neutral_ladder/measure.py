"""Native regional charge conversion on the transported three-character winding register."""

import json
from functools import cache
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_gauge.arrival_ladder import measure as ladder
from model.candidate_laws.interval_echo.openwave.row_11_spectrum.conjugate_flight.measure import (
    charge_reflection,
)
from ...apparatus import native_graph
from ...measure import ONE
from ..instruments import BASES, material, extract
from ..measure import SHOTS, detection_probabilities, sample, witness
from .apparatus import stock

HERE = Path(__file__).resolve().parent
WORDS = ladder.WORDS


@cache
def gate(angle, guard, k=8):
    # Use the directly echoed (x,guard) column for both guard stocks. The
    # material conjugate of the ordinary positive guard is the zero guard.
    energy = ladder.stored_energies(k)
    e = ladder.LAW.echo()
    raw = np.array(ladder.LAW.SRC)
    phase = np.zeros((3, 3), complex)
    t = angle / (8 * ladder.K)
    inv = (-angle / 2 % (2 * np.pi)) / (4 * ladder.K)
    for x in range(3):
        for a in range(3):
            index = ladder.LAW.SRC.index((x, guard))
            dwell = np.exp(-1j * t * (energy[raw[:, 0], a] + energy[raw[:, 1], a]))
            registered = e @ (e[:, index] * np.tile(dwell, 11))
            registered *= np.exp(-1j * inv * energy[x, a])
            assert np.linalg.norm(np.delete(registered, index)) < 4e-13
            phase[x, a] = registered[index]
    opening = expm(-0.5j * np.pi * ladder.free_seam(k))
    closing = expm(-1.5j * np.pi * ladder.free_seam(k))
    u = (
        np.kron(np.eye(3), closing)
        @ np.diag(np.repeat(phase.ravel(), 3))
        @ np.kron(np.eye(3), opening)
    )
    np.testing.assert_allclose(u.conj().T @ u, np.eye(27), atol=4e-13)
    if guard == 2:
        np.testing.assert_allclose(u, ladder.interaction(angle, k), atol=2e-13)
    return u


def word_trials(u, word):
    col = u[:, WORDS.index(word)]
    return abs(col) ** 2


def run():
    apparatus = [stock(sign) for sign in (-1, 1)]
    conjugation = charge_reflection()
    h = native_graph(-1)
    assert (h != native_graph(1)).nnz == 0
    m = material("hub")
    d = material("tip")
    x = np.zeros((h.shape[0], 3), complex)
    x[ONE] = np.eye(3)
    free = expm_multiply(-32j * h, x, traceA=-32j * h.diagonal().sum())
    remote = free[64:128].reshape(8, 8, 3)
    encoding = extract("tip")[0]
    rng = np.random.default_rng(122981)
    banks = []
    ctrl = []
    coherent = []
    coherence_rng = np.random.default_rng(122982)
    for sign in (-1, 1):
        charged = (1 + sign, 1, 1 - sign)
        neutral = (1 + sign, 1 - sign, 1)
        guard = 1 - sign
        for key, angle in [("off", 0.0), ("half", np.pi / 2), ("full", np.pi)]:
            u = gate(angle, guard)
            z = u[WORDS.index(neutral), WORDS.index(charged)]
            prob = abs(z) ** 2
            np.testing.assert_allclose(
                prob, {"off": 0.0, "half": 0.5, "full": 1.0}[key], atol=2e-12
            )
            # Lift to the actual entire 8-head register and three physical
            # output rails, preserving every relative winding phase.
            after = np.kron(u, np.eye(8))
            srcslice = slice(8 * WORDS.index(charged), 8 * (1 + WORDS.index(charged)))
            outslice = slice(8 * WORDS.index(neutral), 8 * (1 + WORDS.index(neutral)))
            operation = after[outslice, srcslice]
            np.testing.assert_allclose(operation, z * np.eye(8), atol=1e-13)
            np.testing.assert_allclose(encoding @ operation, z * encoding, atol=1e-13)
            rows = []
            loads = []
            for basis in BASES:
                source, p_remote, p_read, _ = detection_probabilities(
                    remote, m[basis], d[basis]
                )
                for i in range(3):
                    # Independent signed ternary hub/target/donor acquisition.
                    word_p = word_trials(u, charged)
                    all_words = sample(
                        [
                            *(p_remote[i] * word_p),
                            source[i] - p_remote[i],
                            1 - source[i],
                        ],
                        rng,
                    )
                    loads.append([basis, i, all_words])
                    for j in range(3):
                        phase = p_read[j, i]
                        probs = [
                            prob * phase,
                            prob * (p_remote[i] - phase),
                            (1 - prob) * phase,
                            (1 - prob) * (p_remote[i] - phase),
                            source[i] - p_remote[i],
                            1 - source[i],
                        ]
                        rows.append([basis, i, j, sample(probs, rng)])
            bounds = witness(rows) if key != "off" else None
            if bounds:
                assert bounds["entanglement_fidelity_lower_99"] > 2 / 3
            banks.append(
                dict(
                    sign=sign,
                    setting=key,
                    charged_word=charged,
                    neutral_word=neutral,
                    conversion_probability=float(prob),
                    conversion_amplitude=[float(z.real), float(z.imag)],
                    counts=rows,
                    load_banks=loads,
                    neutral_qutrit_witness=bounds,
                )
            )
        # A donor read halfway through distinguishes the two routes to the
        # final pair and visibly removes their interference.
        half = gate(np.pi / 2, guard)
        j = WORDS.index(charged)
        k = WORDS.index(neutral)
        uninterrupted = abs((half @ half)[k, j]) ** 2
        mid_measured = float((abs(half) ** 2 @ abs(half) ** 2)[k, j])
        np.testing.assert_allclose([uninterrupted, mid_measured], [1, 0.5], atol=2e-12)
        src, at, _, _ = detection_probabilities(remote, m["character"], d["character"])
        for i in range(3):
            for tag, prob in [("held", uninterrupted), ("donor_read", mid_measured)]:
                branches = [
                    at[i] * prob,
                    at[i] * (1 - prob),
                    src[i] - at[i],
                    1 - src[i],
                ]
                coherent.append([sign, i, tag, sample(branches, coherence_rng)])
        for name, word, g, k, target in [
            (
                "opposite_incident",
                (1 - sign, 1, 1 - sign),
                guard,
                8,
                (1 - sign, 1 - sign, 1),
            ),
            ("wrong_guard", charged, 2 - guard, 8, neutral),
            ("depleted_donor", (1 + sign, 1, 1), guard, 8, neutral),
            ("all_filled_hub", (1, 1, 1 - sign), guard, 8, (1, 1 - sign, 1)),
            ("prepared_inverse", neutral, guard, 8, charged),
            ("half_contrast_frozen", charged, guard, 4, neutral),
        ]:
            u = gate(np.pi, g, k)
            p = word_trials(u, word)
            counts = sample(p, rng)
            ctrl.append(
                dict(
                    setting=name,
                    sign=sign,
                    local_input=word,
                    target=target,
                    target_probability=float(p[WORDS.index(target)]),
                    ternary_counts=counts,
                )
            )
    return dict(
        source_and_receiver=[4, 160],
        native_conjugation=conjugation,
        apparatus=apparatus,
        native_wait=32,
        starts_per_bank=SHOTS,
        seed=122981,
        ternary_word_order=WORDS,
        load_count_tail=["no arrival", "source refused"],
        phase_counts=[
            "neutral_and_phase",
            "other_neutral",
            "charged_and_phase",
            "other_charged",
            "no_arrival",
            "source_refused",
        ],
        two_half_pulses=dict(
            seed=122982,
            bins=["neutral", "charged", "no_arrival", "source_refused"],
            counts=coherent,
        ),
        experiments=banks,
        local_controls=ctrl,
    )


if __name__ == "__main__":
    out = run()
    records = out.pop("experiments")
    controls = out.pop("local_controls")
    body = json.dumps(out, indent=2)[:-2]
    for key, rows in [("experiments", records), ("local_controls", controls)]:
        body += (
            ",\n  "
            + json.dumps(key)
            + ": [\n    "
            + ",\n    ".join(json.dumps(r) for r in rows)
            + "\n  ]"
        )
    (HERE / "results.json").write_text(body + "\n}\n")
