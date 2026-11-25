"""Forward and reverse material exposures using an independently drawn local analyzer."""

import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ....beam_response import measure as camera
from ....full_algebra import station
from ..geometry import HERE as BALLS

HERE = Path(__file__).resolve().parent
SHOTS = 2**42
DWELL = 4.0
RELATIVE_STATE_ERROR = 0.000012


def preparation():
    receipt = json.loads(
        (BALLS.parents[1] / "full_algebra" / "results.json").read_text()
    )
    x = np.asarray(receipt["x_real"]) + 1j * np.asarray(receipt["x_imag"])
    prep = receipt["preparation"]
    pauli = [
        (np.eye(2) + sign * a) / 2 for a in station.phase.AXES[1:] for sign in (-1, 1)
    ]
    atoms = np.asarray([np.kron(left, right) for left in pauli for right in pauli])
    kept = np.asarray(prep["kept"]) @ atoms.reshape(36, 16)
    dephased = np.asarray(prep["dephased"]) @ atoms.reshape(36, 16)
    # The sum of the 36 local product Pauli projectors is 9I. Thus even the
    # smaller, dephased half supplies a lower spectral bound of .017 per block.
    assert min(prep["kept"]) > 0 and 4.5 * min(prep["dephased"]) > 0.017
    assert 2.01e-7 / 0.017 < RELATIVE_STATE_ERROR
    source = np.stack([kept + dephased / 2, dephased / 2]).reshape(2, 4, 4)
    np.testing.assert_allclose(source, x, rtol=0, atol=1e-15)
    return source


def propagate(board, route, rings, charge):
    """Evolve launches from each direction as separate columns."""
    addresses = camera.indices(np.r_[rings[0], rings[1]])
    panels = []
    for odd in (False, True):
        h = camera.front.native.generator(board, route, charge, odd)
        assert (h - h.T).count_nonzero() == 0 and np.isrealobj(h.data)
        init = np.zeros((h.shape[0], 40))
        init[addresses] = np.eye(40)
        ht = -(0.5j * DWELL) * h
        out = expm_multiply(ht, init, traceA=ht.diagonal().sum())[addresses]
        panels.append(out.reshape(10, 4, 10, 4).transpose(0, 2, 1, 3))
    return np.asarray(panels)


def transmission(waves, source, analyzer):
    return np.einsum(
        "pai,pyxij,pjk,pyxak->yx", analyzer, waves, source, waves.conj(), optimize=True
    ).real


def count(image, rng):
    all_counts = []
    for col in image.T:
        all_counts.append(rng.multinomial(SHOTS, np.r_[col, 1 - col.sum()]).tolist())
    obs = np.asarray(all_counts)[:, :5].T
    # Marginal exact binomial intervals with a union bound across 200 raw
    # read probabilities and 1 independent colocated standard.
    alpha = 0.01 / (2 * 201)
    low = beta.ppf(alpha, obs, SHOTS - obs + 1)
    high = beta.ppf(1 - alpha, obs + 1, SHOTS - obs)
    return np.nan_to_num(low), np.nan_to_num(high, nan=1.0), all_counts


def gain_bounds(lo, hi):
    # Three homogeneous angular coordinates: det(cross moment)/W³ is the
    # determinant of the centered 2x2 angular covariance (Schur complement).
    from mpmath import iv
    from ... import curved

    iv.dps = 35
    axes = np.asarray(
        [
            [iv.mpf(1), *[iv.mpf(int(x)) / camera.UNIT for x in row]]
            for row in camera.PIXELS
        ],
        dtype=object,
    )
    # Translate both angular coordinates to their held-data centroids to
    # reduce interval overestimation without changing the determinant.
    pivot = (lo + hi) / 2
    w = pivot.sum()
    src = axes.copy()
    dst = axes.copy()
    for j in range(2):
        src[:, j + 1] -= iv.mpf(float(pivot.sum(0) @ camera.BEARINGS[:, j] / w))
        dst[:, j + 1] -= iv.mpf(float(pivot.sum(1) @ camera.BEARINGS[:, j] / w))
    table = np.array(
        [
            [iv.mpf([float(x), float(y)]) for x, y in zip(a, b, strict=True)]
            for a, b in zip(lo, hi, strict=True)
        ],
        dtype=object,
    )
    moments = dst.T @ table @ src
    value = 4 * curved.det_interval(moments) / moments[0, 0] ** 3
    return [float(value.a), float(value.b)]


def run():
    x = preparation()
    board = camera.front.optics.pairing.formation.federation(1)
    route = camera.front.optics.pairing.formation.region(board)
    rings = camera.screens(board, route)
    rng = np.random.default_rng(28112026)
    matched = x.transpose(0, 2, 1)
    identity = np.broadcast_to(np.eye(4), (2, 4, 4))
    unit = float(np.einsum("pij,pji", matched, x).real)
    calibration_count = int(rng.binomial(SHOTS, unit))
    alpha = 0.01 / (2 * 201)
    calibration_interval = [
        float(beta.ppf(alpha, calibration_count, SHOTS - calibration_count + 1)),
        float(beta.ppf(1 - alpha, calibration_count + 1, SHOTS - calibration_count)),
    ]
    assert calibration_interval[0] <= unit <= calibration_interval[1]
    # Quarter turn about left Z permutes all six Pauli test projectors.
    rotate = np.kron(np.diag([1, 1j]), np.eye(2))
    held = rotate @ x @ rotate.conj().T
    held_rng = np.random.default_rng(30281126)
    data = {}
    for q in (-1, 1):
        w = propagate(board, route, rings, q)
        charge = {}
        for key, source, effect, stream in (
            ("plain", x, identity, rng),
            ("matched", x, matched, rng),
            ("held_fixed", held, matched, held_rng),
            ("held_recalibrated", held, held.transpose(0, 2, 1), held_rng),
        ):
            p = transmission(w, source, effect)
            forward, reverse = p[5:, :5], p[:5, 5:]
            row = {}
            enclosures = []
            exact_source = []
            multiplicity = 1 if key == "plain" else 2
            for name, table in (("forward", forward), ("reverse", reverse)):
                lo, hi, counts = count(table, stream)
                np.testing.assert_array_less(lo, table)
                np.testing.assert_array_less(table, hi)
                enclosures.append((lo, hi))
                lower = lo * (1 - RELATIVE_STATE_ERROR) ** multiplicity
                upper = hi * (1 + RELATIVE_STATE_ERROR) ** multiplicity
                exact_source.append((lower, upper))
                j, _ = camera.image(table)
                row[name] = dict(
                    counts=counts,
                    jacobian=j.tolist(),
                    determinant=float(np.linalg.det(j)),
                    sampled_area_gain_interval=gain_bounds(lo, hi),
                    equilibrium_area_gain_interval=gain_bounds(lower, upper),
                )
            # Use only columns 0,1,2 of left ring, every pixel of the right.
            row["unequal_patch_forward"] = float(forward[:, :3].sum() / 3)
            row["unequal_patch_reverse"] = float(reverse[:3, :].sum() / 5)
            delta = forward - reverse.T
            lo, hi = enclosures[0]
            revlo, revhi = enclosures[1]
            row["maximum_pair_asymmetry"] = float(np.max(abs(delta)))
            row["sampled_pair_asymmetry_lower"] = float(
                max(0.0, np.max(lo - revhi.T), np.max(revlo.T - hi))
            )
            lo, hi = exact_source[0]
            revlo, revhi = exact_source[1]
            row["equilibrium_pair_asymmetry_lower"] = float(
                max(0.0, np.max(lo - revhi.T), np.max(revlo.T - hi))
            )
            if key in ("matched", "held_recalibrated"):
                np.testing.assert_allclose(forward, reverse.T, atol=4e-14, rtol=0)
            charge[key] = row
        # Variant needing only local product Pauli analyzers and the publicly
        # selected kept/dephased preparation arm, with no final parity detector.
        kept = x[0] - x[1]
        blank = np.zeros_like(kept)
        kept_normalized = kept.T / np.trace(kept).real
        dep_normalized = x[1].T / np.trace(x[1]).real
        parity_blind = transmission(
            w, np.stack([kept, blank]), np.stack([kept_normalized, kept_normalized])
        )
        parity_blind += transmission(
            w, np.stack([x[1], x[1]]), np.stack([dep_normalized, dep_normalized])
        )
        np.testing.assert_allclose(
            parity_blind[5:, :5], parity_blind[:5, 5:].T, rtol=0, atol=4e-14
        )
        charge["arm_conditioned_area_gain"] = float(
            np.linalg.det(camera.image(parity_blind[5:, :5])[0])
        )
        data[str(q)] = charge
    return dict(
        cells=board.cells,
        left=list(map(int, rings[0])),
        right=list(map(int, rings[1])),
        dwell=DWELL,
        starts_per_launch=SHOTS,
        zero_flight_hit=unit,
        zero_flight_accepted=calibration_count,
        zero_flight_sample_interval=calibration_interval,
        charges=data,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
