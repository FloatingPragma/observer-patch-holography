"""Held transfer from intrinsic seam work and its change during a native flight."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv

from ..common_front import certificate as front
from ..dressed_charge.measure import wave_at
from ..full_algebra import native_precision as native
from ..null_operator.measure import parked_terms
from ..shared_radar.measure import preparations


HERE = Path(__file__).resolve().parent


def receiving_effect(charge):
    """Zero-energy projector of the contacted right seam, on its two endpoints."""
    contacts = native.iso.ranging.parity_blocks(charge, 8, 1)
    np.testing.assert_array_equal(native.iso.ranging.parity_blocks(charge, 0, 1), 0)
    parked = contacts[0, 0]
    np.testing.assert_array_equal(parked, 8 * np.eye(2))
    # Port, left coin, right coin.
    matrix = np.block(
        [[np.kron(np.eye(2), contacts[i, j]) for j in (1, 2)] for i in (1, 2)]
    )
    projector = np.eye(8) - matrix / 8
    np.testing.assert_array_equal(projector @ projector, projector)
    return projector


def source_jets(effects, source):
    """Interval bounds on a fixed linear effect using the registered fixed-volume directions."""
    mid, targets, axis, rots, radius = source
    full = front.full
    entries = [
        e if isinstance(e, sp.MatrixBase) else front.PRIOR.sym_matrix(e)
        for e in effects
    ]
    norm = iv.sqrt(full.ri(sum(full.inner(e, e) for e in entries)))

    def mean(blocks):
        return full.ri(sum(full.inner(b, e) for b, e in zip(blocks, entries)))

    ordinary = mean(targets["ordinary"]) + iv.mpf([-1, 1]) * norm * 2 * iv.sqrt(
        8
    ) * full.ri(radius)
    rotation = sum(
        a * (mean(blocks) + iv.mpf([-1, 1]) * norm * full.ri(radius))
        for a, blocks in zip(axis, rots)
    )
    pairing, pair_rotation = front.variable_coefficients(entries, source)
    angles = json.loads((HERE.parent / "null_optics/certificate.json").read_text())[
        "matched_pairing"
    ]
    for key, lo, hi in (
        ("angle_interval", -297, -291),
        ("ordinary_angle_interval", -997, -983),
    ):
        assert lo < angles[key][0] <= angles[key][1] < hi
    return {
        "ordinary": ordinary + iv.mpf([-997, -983]) * rotation,
        "pairing": pairing + iv.mpf([-297, -291]) * pair_rotation,
    }


def run():
    iv.dps = 30
    source = front.inputs()
    states = preparations()
    arms = json.loads((HERE.parent / "signed_transfers.json").read_text())["tests"]
    angles = {
        ("pairing" if row["source"] == "collar" else "ordinary"): row["angle_per_unit"]
        for row in arms
    }
    paired = native.iso.pairing.formation
    board = paired.federation(1)
    route = paired.region(board)
    np.testing.assert_array_equal(route, [0, 1, 38, 39])
    times = (1 / 32, 1 / 64)
    rows = {}
    each = {}
    rng = np.random.default_rng(192627)
    shots = 2**60
    log_tail = np.log(2 * 2 * len(times) * len(states) / 0.01)
    for q in (-1, 1):
        proj = receiving_effect(q)
        block_effect = []
        exchange_only = []
        shot_effect = [[] for _ in times]
        for odd in (False, True):
            h = native.generator(board, route, q, odd)
            start = np.arange(4)
            end = [4 * p + i for p in route[2:] for i in range(4)]
            assert not h[end, :][:, start].nnz
            second = (h[end, :] @ h[:, start]).toarray()

            # Direct unique 0→1→38 path, also evaluated on the changed generator
            # with all squared-contrast couplings turned off and the reader frozen.
            def path_at_strength(strength):
                contacts = native.iso.ranging.parity_blocks(
                    q, strength, -1 if odd else 1
                )
                from_left = -2 * np.kron(contacts[2, 1], np.eye(2))
                return np.vstack([from_left, np.zeros((4, 4))])

            np.testing.assert_array_equal(second, path_at_strength(8))
            changed = path_at_strength(0)
            exchange_only.append(changed.T @ proj @ changed / 64)
            # H=h/2 up to scalar, evolution starts with -t^2 H^2/2.
            lead = second.T @ proj @ second / 64
            block_effect.append(lead)
            for n, t in enumerate(times):
                wave = wave_at(h, route[0], t)[list(route[2:])].reshape(8, 4)
                shot_effect[n].append(wave.conj().T @ proj @ wave)
        each[q] = np.array(block_effect)
        finite = []
        for time, effect in zip(times, shot_effect):
            readings = {
                name: float(np.einsum("bij,bji->", effect, state).real)
                for name, state in states.items()
            }
            turn = (readings["plus"] - readings["minus"]) / 2
            matched = {
                name: (readings[name] - readings["equilibrium"] + angle * turn)
                / time**4
                for name, angle in angles.items()
            }
            counts = {
                name: int(rng.binomial(shots, probability))
                for name, probability in readings.items()
            }
            estimates = {name: count / shots for name, count in counts.items()}
            errors = {
                name: np.sqrt(2 * probability * log_tail / shots) + 2 * log_tail / shots
                for name, probability in estimates.items()
            }
            samples = {}
            for name, angle in angles.items():
                measured = (
                    estimates[name]
                    - estimates["equilibrium"]
                    + angle * (estimates["plus"] - estimates["minus"]) / 2
                ) / time**4
                bound = (
                    errors[name]
                    + errors["equilibrium"]
                    + abs(angle) * (errors["plus"] + errors["minus"]) / 2
                ) / time**4
                samples[name] = {"estimate": measured, "half_width": bound}
            finite.append(
                dict(
                    time=time,
                    zero_energy_probability=readings,
                    matched_per_time_four=matched,
                    counts=counts,
                    count_slopes=samples,
                )
            )
        rows[str(q)] = dict(
            finite=finite,
            exchange_only_reference_coefficient=float(
                np.einsum("bij,bji->", exchange_only, states["equilibrium"]).real
            ),
            active_reference_coefficient=float(
                np.einsum("bij,bji->", block_effect, states["equilibrium"]).real
            ),
        )
    each[0] = (each[-1] + each[1]) / 2
    slopes = {}
    for q, e in each.items():
        jets = source_jets(e, source)
        entropy_ordinary = iv.mpf([".4094266", ".4094802"])
        entropy_pairing = iv.mpf([".2310794", ".2311330"])
        gain = jets["pairing"] * entropy_ordinary / (jets["ordinary"] * entropy_pairing)
        assert (
            jets["ordinary"].b < 0
            and jets["pairing"].b < 0
            and gain.a > 0
            and gain.b < 1
        )
        slopes[str(q)] = dict(
            leading_effect=e.tolist(),
            ordinary=str(jets["ordinary"]),
            pairing=str(jets["pairing"]),
            held_transfer=str(gain),
        )
    return dict(
        shots_per_setting=shots,
        seed=192627,
        parked_native_terms=parked_terms(),
        launch=route[0],
        receiving_seam=route[2:],
        coefficient_slopes=slopes,
        exposures=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(report)
    else:
        print(report, end="")
