"""Independently acquired material Bell-coincidence records on the two native cameras."""

import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

from . import apparatus as gates
from .. import measure as screen

HERE = Path(__file__).resolve().parent
STARTS = 2**46
OUTCOMES = ("parity_mismatch", "both", "left_only", "right_only", "neither")
ALPHA = 0.01 / (2 * (2 * 3 * 2 * 25 + 3))


def store_map(odd, shuttle):
    """Store the arrived coin on a fresh auxiliary even mark via two SWAPs."""
    first = gates.extractor(gates.contact(odd), odd)
    put = gates.expm(-1j * gates.pi * gates.contact(0) / 8).reshape(2, 2, 2, 2)
    # Incoming rail arbitrary; auxiliary coin starts in |0>. Final rail is 0.
    np.testing.assert_allclose(put[1, :, :, 0], 0, atol=2e-13)
    result = put[0, :, :, 0] @ shuttle @ first
    np.testing.assert_allclose(result, np.eye(2), atol=2e-13)
    return result


def pair_effects():
    identity = np.eye(4)
    shuttles = [gates.shuttle(m) for m in gates.geometry()]
    right = gates.bell_effect(0, 0)
    table = []
    for odd in (0, 1):
        # The saved target's new mark is even, whatever its original parity.
        left = gates.bell_effect(odd, 0)
        local = np.kron(
            store_map(odd, shuttles[0]) @ gates.parity_read(odd)[odd],
            store_map(0, shuttles[1]),
        )
        lift = np.kron(
            np.kron(shuttles[0] @ gates.parity_read(odd)[odd], shuttles[1]), local
        )
        table.append(
            [
                lift.conj().T @ gates.paired_effect(a, b) @ lift
                for a, b in (
                    (left, right),
                    (left, identity - right),
                    (identity - left, right),
                    (identity - left, identity - right),
                )
            ]
        )
    return np.asarray(table)


def instrument(waves, source, ref, effects):
    rho = waves @ source[:, None, None, :, :] @ waves.conj().swapaxes(-1, -2)
    hits = np.trace(rho, axis1=-2, axis2=-1).real
    probability = np.trace(ref, axis1=-2, axis2=-1).real
    mismatch = np.einsum("pyx,p->yx", hits, probability[::-1])
    # Retain joint off-diagonal coin amplitudes when contracting each of
    # the four binary Bell outcomes. Neither head is separately averaged.
    coincidences = np.einsum(
        "poabcd,pca,pyxdb->yxo",
        effects.reshape(2, 4, 4, 4, 4, 4),
        ref,
        rho,
        optimize=True,
    ).real
    combined = np.concatenate([mismatch[:, :, None], coincidences], axis=-1)
    np.testing.assert_allclose(combined.sum(-1), hits.sum(0), rtol=0, atol=2e-13)
    np.testing.assert_allclose(
        combined[:, :, 1],
        screen.transmission(waves, source, ref.transpose(0, 2, 1)) / 4,
        rtol=0,
        atol=2e-14,
    )
    assert np.min(combined) > -2e-15
    return np.maximum(combined, 0)


def acquire(p, rng):
    records = []
    for source in range(p.shape[1]):
        outcomes = p[:, source].ravel()
        all_events = np.r_[outcomes, 1 - outcomes.sum()]
        counts = rng.multinomial(STARTS, all_events)
        records.append(counts.tolist())
    hits = np.array(records)[:, :-1].reshape(-1, 5, len(OUTCOMES))[:, :, 1].T
    lo = beta.ppf(ALPHA, hits, STARTS - hits + 1)
    hi = beta.ppf(1 - ALPHA, hits + 1, STARTS - hits)
    return records, np.nan_to_num(lo), np.nan_to_num(hi, nan=1.0)


def controls(effects):
    eye = np.eye(4, dtype=complex)[None, None, None]
    unit = np.broadcast_to(eye, (2, 1, 1, 4, 4))

    def wrap(vector):
        vector = np.asarray(vector)
        return np.stack([np.outer(vector, vector.conj()), np.zeros((4, 4))])

    zero = wrap([1, 0, 0, 0])
    one = wrap([0, 0, 0, 1])
    mix = (zero + one) / 2
    independent = instrument(unit, mix, mix, effects)[0, 0, 1]
    copied = (
        instrument(unit, zero, zero, effects)[0, 0, 1]
        + instrument(unit, one, one, effects)[0, 0, 1]
    ) / 2
    dark = wrap(np.array([1, 0, 1j, 0]) / np.sqrt(2))
    correlated = (np.eye(4) + np.kron(gates.X, gates.X)) / 4
    uncorrelated = np.eye(4) / 4
    corr = np.stack([correlated, np.zeros((4, 4))])
    separate = np.stack([uncorrelated, np.zeros((4, 4))])
    return dict(
        off_diagonal_xx_live=float(instrument(unit, corr, corr, effects)[0, 0, 1]),
        off_diagonal_xx_marginals_only=float(
            instrument(unit, corr, separate, effects)[0, 0, 1]
        ),
        independent_shared_marginals=float(independent),
        copied_lottery_label=float(copied),
        pure_left_y_zero_flight=float(instrument(unit, dark, dark, effects)[0, 0, 1]),
    )


def run():
    receipt = gates.report()
    effects = pair_effects()
    x = screen.preparation()
    turn = np.kron(np.diag([1, 1j]), np.eye(2))
    held = turn @ x @ turn.conj().T
    rng = np.random.default_rng(6261126)
    receipt.update(
        controls=controls(effects),
        starts_per_launch=STARTS,
        detector_outcomes=list(OUTCOMES) + ["outside_all_five"],
        dwell=screen.DWELL,
    )
    stocks = (("equilibrium", x, x), ("held_live", held, held), ("held_stale", held, x))
    calibration = {}
    identity = np.broadcast_to(np.eye(4), (2, 1, 1, 4, 4))
    for key, source, ref in stocks:
        u = float(instrument(identity, source, ref, effects)[0, 0, 1])
        count = int(rng.binomial(STARTS, u))
        interval = [
            float(beta.ppf(ALPHA, count, STARTS - count + 1)),
            float(beta.ppf(1 - ALPHA, count + 1, STARTS - count)),
        ]
        assert interval[0] < u < interval[1]
        calibration[key] = dict(
            population=u, accepted=count, starts=STARTS, sampled_99=interval
        )
    receipt["zero_flight_standard"] = calibration
    board = screen.camera.front.optics.pairing.formation.federation(1)
    route = screen.camera.front.optics.pairing.formation.region(board)
    rings = screen.camera.screens(board, route)
    charges = {}
    for q in (-1, 1):
        w = screen.propagate(board, route, rings, q)
        reports = {}
        for key, source, ref in stocks:
            p = instrument(w, source, ref, effects)
            rows = {}
            bounds = []
            for label, table in (("forward", p[5:, :5]), ("reverse", p[:5, 5:])):
                counts, lo, hi = acquire(table, rng)
                hits = table[:, :, 1]
                assert np.all(lo <= hits) and np.all(hits <= hi)
                j, _ = screen.camera.image(hits)
                rows[label] = dict(
                    counts=counts,
                    population_area_gain=float(np.linalg.det(j)),
                    count_area_gain_interval=screen.gain_bounds(lo, hi),
                )
                bounds.append((lo, hi))
            a, b = bounds[0]
            c, d = bounds[1]
            rows["count_directional_asymmetry_lower"] = float(
                max(0.0, (a - d.T).max(), (c.T - b).max())
            )
            rows["population_maximum_asymmetry"] = float(
                np.max(abs(p[5:, :5, 1] - p[:5, 5:, 1].T))
            )
            reports[key] = rows
        charges[str(q)] = reports
    receipt["independent_flight_and_reference_pairs"] = (10 * 2 * 3 + 3) * STARTS
    receipt["charges"] = charges
    return receipt


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
