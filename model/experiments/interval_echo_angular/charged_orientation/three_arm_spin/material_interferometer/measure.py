"""Cross the charge's own interferometer with an independent winding pointer."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ..apparatus import native_graph
from ..measure import BASE, CY, ONE, Y, build
from ..recoil_search.measure import clock
from .apparatus import Terminal

SHOTS = 2**38
SENSES = (1, -1)


def interval(k, shots=SHOTS):
    # Fewer than 128 marginal rates including training and both holdouts.
    tail = 0.01 / (2 * 128)
    return [
        float(beta.ppf(tail, k, shots - k + 1)) if k else 0,
        float(beta.isf(tail, k + 1, shots - k)) if k < shots else 1,
    ]


def sample(p_src, p_remote, p_match, rng):
    probs = np.array([p_match, p_remote - p_match, p_src - p_remote, 1 - p_src])
    assert min(probs) > -1e-14
    counts = rng.multinomial(SHOTS, np.maximum(probs, 0) / np.maximum(probs, 0).sum())
    return dict(
        counts=counts.tolist(),
        match_probability=float(p_match),
        remote_probability=float(p_remote),
        match_interval99=interval(int(counts[0])),
        remote_interval99=interval(int(counts[:2].sum())),
    )


def forecast(rows, narrow=False):
    previous = Path(__file__).resolve().parents[1]
    if narrow:
        bank = json.loads((previous / "results.json").read_text())["filtered"]
        train_n = bank["starts_per_axis"]
        # Choi counts compare identical analyzer axes on both reference legs.
        # A Y eigenpreparation is the conjugate of the retained reference Y.
        assert bank["counts"][1]["axis"] == "Y"
        same_ref, opposite_ref = bank["counts"][1]["counts"][:2]
        recorded = {
            (x, y): opposite_ref if x == y else same_ref for x in SENSES for y in SENSES
        }
    else:
        bank = json.loads((previous / "reciprocal.json").read_text())["forward"]
        train_n = bank["shots_per_setting"]
        recorded = {
            (-row["input_sign"], -sgn): n
            for row in bank["observations"]
            if row["input_axis"] == 1
            for sgn, n in zip(SENSES, row["settings"][1][:2], strict=True)
        }
    for row in rows:
        expected = interval(recorded[row["emitted"], row["analyzed"]], train_n)
        seen = row["match_interval99"]
        row["pointer_forecast99"] = expected
        row["forecast_overlap"] = bool(
            max(seen[0], expected[0]) <= min(seen[1], expected[1])
        )
    return sum(row["forecast_overlap"] for row in rows)


def narrow_banks(sources, detectors):
    previous = Path(__file__).resolve().parents[1]
    spec = json.loads((previous / "results.json").read_text())["clock"]
    ev, _, weights, q, v, _, perms = build()
    wave = v @ (weights * clock(ev, spec))
    first = q @ wave
    cyc = perms[CY]
    second = (first[cyc] - first[cyc[cyc]]) / np.sqrt(3)
    check = BASE[:, 1] @ ((q[cyc[ONE]] - q[cyc[cyc[ONE]]]) @ v) / np.sqrt(3)
    if check @ weights < 0:
        second *= -1
    response = np.column_stack([first, second])
    projection = np.zeros((8, 2))
    projection[ONE] = BASE
    shots = 2**42
    rng = np.random.default_rng(632476)
    rows = []
    for src, apertures in enumerate(sources):
        for i, x in enumerate(SENSES):
            prep = apertures[i]
            coordinates = projection.T @ prep
            # The numerical material pulse has tiny leakage outside the E row.
            # Every subsequent clock, projection and decoder is contractive.
            err = np.linalg.norm(prep - projection @ coordinates)
            assert err < 1e-13
            filtered = response @ coordinates
            admitted = np.linalg.norm(filtered) ** 2
            slab = filtered[64:128].reshape(8, 8)
            remote = np.linalg.norm(slab) ** 2
            p_src = np.linalg.norm(prep) ** 2
            for detector, decoded in enumerate(detectors):
                for j, y in enumerate(SENSES):
                    matched = np.linalg.norm(decoded[j].conj() @ slab) ** 2
                    probs = [
                        matched,
                        remote - matched,
                        admitted - remote,
                        p_src - admitted,
                        1 - p_src,
                    ]
                    counts = rng.multinomial(shots, probs)
                    rows.append(
                        dict(
                            source=src,
                            detector=detector,
                            emitted=x,
                            analyzed=y,
                            counts=counts.tolist(),
                            match_probability=float(matched),
                            remote_probability=float(remote),
                            clock_probability=float(admitted),
                            projected_input_norm_error=float(err),
                            clock_interval99=interval(int(counts[:3].sum()), shots),
                            match_interval99=interval(int(counts[0]), shots),
                            remote_interval99=interval(int(counts[:2].sum()), shots),
                        )
                    )
    return dict(
        clock=spec,
        virtual_starts=shots,
        count_order=[
            "matched phase",
            "other receiving charge",
            "no receiving charge after clock",
            "clock refusal",
            "source refusal",
        ],
        forecasts_overlapping=forecast(rows, narrow=True),
        experiments=rows,
    )


def run():
    chart = [Terminal("hub"), Terminal("tip")]
    local = [t.ports_and_maps() for t in chart]
    pointer = np.zeros((2, 8), complex)
    for i, sign in enumerate(SENSES):
        pointer[i, ONE] = np.exp(2j * np.pi * sign * np.arange(3) / 3) / 3
        coord = BASE.T @ pointer[i, ONE]
        np.testing.assert_allclose(Y @ coord, -sign * coord, atol=2e-15)
    # Compare native encoders and decoders as quantum maps up to one common
    # mounting phase (which drops out of the local projective reads).
    errors = []
    for e in local:
        phases = np.einsum("ic,ic->i", pointer.conj(), e) * 3
        errors.append(float(np.max(abs(e - phases[:, None] * pointer))))
        assert max(errors) < 2e-13 and np.max(abs(abs(phases) - 1)) < 2e-13
    sources = [pointer, local[0]]
    detectors = [pointer, local[1]]
    # One joint sparse flight for four independently prepared source columns
    # and three sharp-word source controls. All columns include admission.
    h = native_graph()
    states = np.zeros((h.shape[0], 7), complex)
    states[:8, :4] = np.concatenate(sources).T
    states[ONE, range(4, 7)] = 1 / 3
    gen = -32j * h
    result = expm_multiply(gen, states, traceA=gen.diagonal().sum())
    remote = result[64:128].reshape(8, 8, 7)
    p_src = (abs(states) ** 2).sum(axis=0)
    p_remote = (abs(remote) ** 2).sum(axis=(0, 1))
    rng = np.random.default_rng(632475)
    rows, mixed = [], []
    for d, analysis in enumerate(detectors):
        amps = np.einsum("qr,rhc->qhc", analysis.conj(), remote)
        matches = (abs(amps) ** 2).sum(axis=1)
        for s in range(2):
            for i, x in enumerate(SENSES):
                for j, y in enumerate(SENSES):
                    k = 2 * s + i
                    rows.append(
                        dict(
                            source=s,
                            detector=d,
                            emitted=x,
                            analyzed=y,
                            **sample(p_src[k], p_remote[k], matches[j, k], rng),
                        )
                    )
        for j, sign in enumerate(SENSES):
            mixed.append(
                dict(
                    detector=d,
                    analyzed=sign,
                    **sample(
                        p_src[4:].sum(), p_remote[4:].sum(), matches[j, 4:].sum(), rng
                    ),
                )
            )
    return dict(
        duration=32,
        virtual_starts=SHOTS,
        banks=["three-path pointer", "charge-path contacts"],
        counts_order=[
            "matched local phase",
            "other receiving charge",
            "no receiving charge",
            "source refusal",
        ],
        terminals=[t.ports for t in chart],
        max_aperture_errors=errors,
        forecasts_overlapping=forecast(rows),
        experiments=rows,
        dephased=mixed,
        filtered=narrow_banks(sources, detectors),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(run(), indent=2) + "\n"
    if args.write:
        Path(__file__).with_name("results.json").write_text(out)
    else:
        print(out, end="")
