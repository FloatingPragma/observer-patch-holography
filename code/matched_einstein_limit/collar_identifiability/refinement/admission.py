"""Sparse admission on the complete return register of a formed native collar."""

import json
from math import atanh, expm1, isqrt, log1p, pi

import numpy as np

from model.experiments.interval_echo_einstein_station.collar_pairing import measure as source


def information_bounds(reference, delta):
    """Symmetric-logarithmic-derivative maximum for one classical record and the full quantum log."""
    weights, axes = np.linalg.eigh(reference)
    d = axes.conj().T @ delta @ axes
    L = axes @ (2*d / (weights[:, None] + weights[None, :])) @ axes.conj().T
    np.testing.assert_allclose((reference @ L + L @ reference) / 2, delta, rtol=0, atol=3e-14)
    max_classical = np.trace(reference @ L @ L).real
    # The ideal projective read of L is a mathematical upper control; it is not a source instrument.
    _, vectors = np.linalg.eigh(L)
    p = np.diag(vectors.conj().T @ reference @ vectors).real
    shift = np.diag(vectors.conj().T @ delta @ vectors).real
    np.testing.assert_allclose((shift**2/p).sum(), max_classical, atol=3e-13, rtol=2e-11)
    full, _ = source.relative_curvature(reference, delta)
    return float(full), float(max_classical)


def material_panel(level, sign, strength):
    board = source.formation.federation(level)
    route = source.formation.region(board)
    joint = source.prepared_states(board, route, sign)[0]
    joint = joint / np.trace(joint).real
    reference = np.kron(*source.formation.marginals(joint))
    full, ceiling = information_bounds(reference, joint-reference)
    h = source.formation.hamiltonian(board, route, 1, strength)
    arm, = source.material_read(h, route[1:3], {sign: (joint, reference)})
    records = []
    for paired, fresh in zip(arm["formed"]["appointments"], arm["independent"]["appointments"], strict=True):
        outcomes = ("relay_miss", "return_miss", "complete")
        p, q = (np.array([row[key] for key in outcomes]) for row in (paired, fresh))
        d = p - q
        slope = float(-d @ np.log(q))
        coefficient = float((d**2/q).sum())
        assert coefficient <= ceiling + 1e-13
        records.append(dict(time=paired["time"], q=q.tolist(), d=d.tolist(),
                            linear_modular=slope, quartic_information=coefficient))
    return dict(cells=board.cells, source_sign=sign, probe_strength=strength,
                quantum_log=full, one_copy_classical_ceiling=ceiling, records=records)


def weaken(q, d, hit):
    """Signed, product-background-subtracted modular slope in the weakening direction."""
    return float(hit * (d @ np.log1p(hit * d/q)))


def compare_returns(source_return, reference_return):
    """Ties use one fair material bit; outcomes 0 and 1 in each input are exhaustive."""
    ties = source_return*reference_return + (1-source_return)*(1-reference_return)
    return source_return*(1-reference_return) + ties/2


def balanced_charge(q, d, hit):
    shift = hit * d[-1]  # Last record is the completed return.
    np.testing.assert_allclose(compare_returns(q[-1] + shift, q[-1]), (1 + shift)/2, atol=2e-15)
    return float(shift * atanh(shift))


def sphere_count(radius):
    """All cubic lattice points in a closed Euclidean ball, counted by vertical chords."""
    n = 0
    for a in range(-radius, radius + 1):
        bmax = isqrt(radius**2 - a**2)
        for b in range(-bmax, bmax + 1):
            n += 1 + 2*isqrt(radius**2 - a**2 - b**2)
    return n


def any_hit(probability, trials):
    return -expm1(trials * log1p(-probability))


def flat_panel(record):
    q, d = (np.array(record[key]) for key in ("q", "d"))
    rows = []
    for n in (8, 16, 32, 64, 128):
        ell = 1/n
        volume_fraction = sphere_count(n) / n**6
        regular = any_hit(volume_fraction, n)
        rows.append(dict(inverse_radius=n, per_mark_hit=volume_fraction, any_hit=regular,
                         one_mark_density=weaken(q, d, volume_fraction)/ell**4,
                         growing_bank_density=weaken(q, d, regular)/ell**4,
                         balanced_one_mark=balanced_charge(q, d, volume_fraction)/ell**4,
                         balanced_growing_bank=balanced_charge(q, d, regular)/ell**4))
    return dict(limiting_growing_bank=(4*pi/3)**2 * record["quartic_information"],
                limiting_balanced_growing_bank=float((4*pi*d[-1]/3)**2), rows=rows)


def run():
    measurements = [material_panel(level, sign, 8) for level in (0, 1) for sign in (-1, 1)]
    measurements.append(material_panel(1, 1, 0))
    return dict(material=measurements, flat_admission=flat_panel(measurements[3]["records"][-1]))


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
