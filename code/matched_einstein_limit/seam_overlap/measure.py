"""Independently meter a native modular contrast and test its curvature in an ADM completion."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from model.experiments.interval_echo_einstein_station import backreaction as seam
from .. import observer_field as radar

HERE = Path(__file__).resolve().parent
V0 = (2 * math.pi)**3
BETA = 1/16
SHOTS = 2**22
SIGNS = (-0.5, 0.5)
CYCLE = np.array([[1., -1.], [-1., 1.]])
APPOINTMENTS = (0.35, 0.6, 0.85, 1.1)
SPACINGS = (0.04, 0.02, 0.01)


def native_energies(phase):
    exchange, agreement = seam.law.local_term_parts()
    h = exchange + seam.KAPPA * agreement
    loads = np.eye(seam.law.N_SRC)
    basis = [seam.collar_vector(magnitude, 1, phase) for magnitude in (0, 2)]
    energies = np.empty((2, 2))
    for i, word in enumerate(basis):
        for j, occupancy in enumerate([(2, 0), (1, 1)]):
            state = np.kron(word, loads[seam.law.SRC.index(occupancy)])
            energy = state @ h @ state
            np.testing.assert_allclose(h @ state, energy * state, atol=2e-13, rtol=0)
            energies[i, j] = energy
    vacuum = loads[seam.law.SRC.index((0, 0))]
    symmetric = (loads[seam.law.SRC.index((1, 0))] + loads[seam.law.SRC.index((0, 1))])/np.sqrt(2)
    ordinary = []
    for load in (vacuum, symmetric):
        state = np.kron(np.eye(11)[0], load)
        energy = state @ h @ state
        np.testing.assert_allclose(h @ state, energy * state, atol=2e-13, rtol=0)
        ordinary.append(energy)
    return np.asarray(ordinary), energies


def canonical_source(energy):
    # q = log(a), p is the sum of the three diagonal canonical momenta.
    kappa, vacuum = radar.COUPLING, radar.VACUUM
    p0 = -math.sqrt(12*V0*(vacuum*V0/kappa + energy)/kappa)

    def rates(_, state):
        q, p, _ = state
        volume = V0*math.exp(3*q)
        return [-kappa*p/(6*volume), -kappa*p*p/(4*volume)-3*vacuum*volume/kappa, math.exp(-q)]

    sol = solve_ivp(rates, (0, 1.5), (0., p0, 0.), dense_output=True,
                    method='DOP853', rtol=2e-13, atol=1e-13)
    if not sol.success:
        raise RuntimeError(sol.message)

    def optical_state(time):
        q, _, distance = sol.sol(time)
        return np.r_[np.full(3, math.exp(q)), np.zeros(5), np.full(3, distance)]
    return optical_state


def tapes(energy):
    source = canonical_source(energy)
    return [[radar.records(source, t, step) for t in APPOINTMENTS] for step in SPACINGS]


def geometry(tape):
    _, lengths, *_ = radar.local_state(tape)
    temporal, spatial, *_ = radar.reconstruct(tape)
    return V0 * np.prod(lengths), temporal, spatial


def gibbs(energies):
    weights = np.exp(-BETA*energies)
    return weights/weights.sum()


def draw(rng, probs):
    return rng.multinomial(SHOTS, probs.ravel()).reshape(probs.shape)


def read_overlap(reference, high, low, beta):
    sigma = reference/reference.sum()
    difference = high/high.sum() - low/low.sum()
    local = -np.log(np.outer(sigma.sum(axis=1), sigma.sum(axis=0)))
    overlap = -np.log(sigma) - local
    return float(np.sum(overlap * difference)/beta), float(np.sum(local * difference)/beta)


def counting_error(reference, high, low, thermometer, gap, beta):
    """Delta-method standard error for independently sampled multinomial banks."""
    sigma, a, b = [bank / bank.sum() for bank in (reference, high, low)]
    d = a - b
    rows, cols = sigma.sum(axis=1), sigma.sum(axis=0)
    influence = (-d/sigma + (d.sum(axis=1)/rows)[:, None] + d.sum(axis=0)/cols) / beta
    charge = np.log(np.outer(rows, cols)/sigma) / beta

    def variance(p, x, count):
        return np.sum(p*(x - np.sum(p*x))**2)/count

    term = variance(sigma, influence, reference.sum())
    term += variance(a, charge, high.sum()) + variance(b, charge, low.sum())
    beta_variance = (1/thermometer).sum()/gap**2
    return float(np.sqrt(term + (np.sum(d*charge)/beta)**2 * beta_variance))


def compare(ordinary, calibration, high, low, matter):
    rows = []
    for index, spacing in enumerate(SPACINGS):
        v0, g0, _ = geometry(calibration[0][index][0])
        v1, g1, _ = geometry(calibration[1][index][0])
        kappa = (g1 - g0)/(ordinary[1]/v1 - ordinary[0]/v0)
        vacuum = g0 - kappa*ordinary[0]/v0
        holds = []
        for t, a, b in zip(APPOINTMENTS[1:], high[index][1:], low[index][1:], strict=True):
            va, ga, sa = geometry(a)
            vb, gb, sb = geometry(b)
            temporal = va*(ga-vacuum) - vb*(gb-vacuum)
            null = va*(ga+sa) - vb*(gb+sb)
            holds.append(dict(time=t, integrated_temporal=float(temporal),
                              integrated_axial_null=null.tolist(), forecast=float(kappa*matter)))
        rows.append(dict(spacing=spacing, ordinary_coupling=float(kappa),
                         vacuum=float(vacuum), held=holds))
    return rows


def run():
    rng = np.random.default_rng(142634)
    ordinary, _ = native_energies(1)
    calibration = [tapes(e) for e in ordinary]
    thermometer = draw(rng, gibbs(ordinary))
    beta = math.log(thermometer[0]/thermometer[1])/(ordinary[1]-ordinary[0])
    probs = [(1+s*CYCLE)/4 for s in SIGNS]
    arms = {}
    for label, phase in [('interacting', 1), ('protected', -1)]:
        _, h = native_energies(phase)
        reference = draw(rng, gibbs(h))
        high, low = [draw(rng, p) for p in probs]
        energy = [float(np.sum(p*h)) for p in probs]
        acquisitions = [tapes(e) for e in energy]
        overlap, local = read_overlap(reference, high, low, beta)
        exact, _ = read_overlap(gibbs(h), *probs, BETA)
        arms[label] = dict(spectrum=h.tolist(), reference_counts=reference.tolist(),
                           high_counts=high.tolist(), low_counts=low.tolist(), source_energies=energy,
                           overlap_meter=overlap, marginal_meter=local, exact_overlap_energy=exact,
                           overlap_counting_se=counting_error(reference, high, low, thermometer,
                                                             ordinary[1]-ordinary[0], beta),
                           high_tapes=acquisitions[0], low_tapes=acquisitions[1],
                           comparisons=compare(ordinary, calibration, *acquisitions, overlap))
    # A different gravitational charge, constant on the four interacting words.
    # Spectroscopy, the Gibbs reference and calibration on the blank sector are unchanged.
    common_energy = float(native_energies(1)[1].mean())
    common_tapes = tapes(common_energy)
    control = dict(gravitating_energy=common_energy, shared_tapes=common_tapes,
                   comparisons=compare(ordinary, calibration, common_tapes, common_tapes,
                                       arms['interacting']['overlap_meter']))
    return dict(selective_gravity_control=control, ordinary_energies=ordinary.tolist(), thermometer_counts=thermometer.tolist(),
                beta_from_ordinary=float(beta), shots_per_setting=SHOTS,
                calibration_tapes=calibration, arms=arms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    data = run()
    path = HERE / 'results.json'
    if parser.parse_args().write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    else:
        from cosmology.rar_deep_regime.lensing_profile import check
        check(json.loads(path.read_text()), data)
    for name, arm in data['arms'].items():
        print(name, arm['overlap_meter'], arm['comparisons'][-1])
