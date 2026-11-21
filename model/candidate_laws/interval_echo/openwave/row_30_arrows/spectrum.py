"""Spatially projected free density records, with a separate held temporal window."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.candidate_laws.interval_echo import candidate
from .experiment import local_check, sector


HERE = Path(__file__).resolve().parent
TARGETS = (.12, .32, 2.9, 5, 7)
TIMES = np.arange(769) / 24
TRAIN = TIMES <= 16
SHOTS = 65536


def harmonic_ports(car):
    values, vectors = np.linalg.eigh(car.laplacian().toarray())
    probes = []
    for target in TARGETS:
        eig = values[np.argmin(abs(values - target))]
        space = vectors[:, abs(values - eig) < 1e-9]
        projection = space @ space[0]
        norm = np.linalg.norm(projection)
        weights = projection / norm if norm > 1e-10 else np.zeros(car.n)
        assert np.linalg.norm(car.laplacian() @ weights - eig * weights) < 1e-10
        probes.append((target, float(eig), weights))
    assert len(set(round(x[1], 7) for x in probes)) == len(probes)
    return probes


def check_head_detector():
    """Verify that the native generator preserves local occupation minus incident heads."""
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    h = h0 + candidate.KAPPA * h1
    labels = [(coin, pair) for coin in (0, 1) for pair in ((0, 0), (1, 0), (0, 1))]
    embedded = np.zeros((law.DIM, len(labels)))
    for column, (coin, pair) in enumerate(labels):
        for d in (-1, 1):
            native = (1 + law.LABELS.index((coin, d))) * law.N_SRC + law.SRC.index(pair)
            embedded[native, column] = 1 / np.sqrt(2)
    block = embedded.T @ h @ embedded
    assert np.max(abs(h @ embedded - embedded @ block)) < 1e-11
    for endpoint in (0, 1):
        occupation = np.array([pair[endpoint] for _, pair in labels])
        head = np.array([int(coin == 1-endpoint) for coin, _ in labels])
        difference = occupation - head
        assert np.max(abs(block * (difference[:, None] - difference))) < 1e-11
        # Both terms change on a real hop; commutation of their difference has content.
        assert np.max(abs(block * (occupation[:, None] - occupation))) > 1


def density_records(h, sites, seed):
    active, indices = np.unique(sites, return_inverse=True)
    initial = np.zeros(len(sites))
    initial[0] = 1
    states = expm_multiply(-1j * h, initial, start=0, stop=32, num=len(TIMES))
    p = np.array([np.bincount(indices, weights=abs(row)**2, minlength=len(active)) for row in states])
    if max(abs(p.sum(axis=1)-1)) > 1e-10:
        raise ArithmeticError("density detector lost a preparation")
    rng = np.random.default_rng(seed)
    counts = np.array([rng.multinomial(SHOTS, row / row.sum()) for row in p])
    return active, p, counts


def single_frequency(values, max_se):
    train_y = values[TRAIN]
    train_t = TIMES[TRAIN]
    grid = np.linspace(np.pi / 4, 56, 1850)
    phases = np.outer(grid, train_t)
    cosine, sine = np.cos(phases), np.sin(phases)
    c = cosine - cosine.mean(axis=1, keepdims=True)
    s = sine - sine.mean(axis=1, keepdims=True)
    yc = train_y - train_y.mean()
    cy, sy = c @ yc, s @ yc
    cc, ss, cs = (c*c).sum(axis=1), (s*s).sum(axis=1), (s*c).sum(axis=1)
    improvement = (ss*cy**2 - 2*cs*cy*sy + cc*sy**2) / (cc*ss - cs**2)
    j = int(np.argmax(improvement))

    def solve(omega):
        basis = np.column_stack([np.ones(len(TIMES)), np.cos(omega*TIMES), np.sin(omega*TIMES)])
        params = np.linalg.lstsq(basis[TRAIN], train_y, rcond=None)[0]
        predicted = basis @ params
        return float(np.sum((predicted[TRAIN] - train_y)**2)), predicted, params

    interval = (grid[max(j-1, 0)], grid[min(j+1, len(grid)-1)])
    best = minimize_scalar(lambda x: solve(x)[0], bounds=interval, method='bounded')
    omega = min([float(best.x), float(grid[j])], key=lambda x: solve(x)[0])
    _, predicted, params = solve(omega)
    # Use the best constant *on each half* for the variance denominator.
    r2 = []
    for m in (TRAIN, ~TRAIN):
        ref = float(np.sum((values[m] - values[m].mean())**2))
        r2.append(1 - float(np.sum((values[m] - predicted[m])**2)) / ref if ref > 1e-20 else None)
    amp = float(np.hypot(params[1], params[2]))
    resolved = all(v is not None and v >= .9 for v in r2) and amp > 5*max_se
    return dict(omega=round(omega, 7), amplitude=round(amp, 9), constant=round(float(params[0]), 9),
                train_r2=None if r2[0] is None else round(r2[0], 7),
                holdout_r2=None if r2[1] is None else round(r2[1], 7),
                max_standard_error=round(max_se, 9), resolved=resolved)


def shell(rows):
    passing = [r for r in rows if r['count_fit']['resolved'] and r['exact_fit']['resolved']]
    if len(passing) < 3:
        return dict(passing_spatial_modes=len(passing), status='insufficient resolved frequencies')
    x = np.array([r['lambda'] for r in passing])
    frequency = np.array([r['count_fit']['omega'] for r in passing])
    slope = (frequency[1]**2-frequency[0]**2) / (x[1]-x[0])
    intercept = frequency[0]**2 - slope*x[0]
    if slope <= 0 or intercept <= 0 or np.min(np.diff(np.sort(frequency))) < 2*np.pi/16:
        return dict(passing_spatial_modes=len(passing), status='uncalibrated or unresolved mass shell')
    errors = abs(np.sqrt(intercept+slope*x[2:])/frequency[2:]-1)
    return dict(passing_spatial_modes=len(passing), rest_squared=round(float(intercept), 8),
                speed_squared=round(float(slope), 8), held_relative_error=np.round(errors, 8).tolist(),
                status='density envelope passed' if max(errors) <= .05 else 'density envelope failed')


def run():
    car = federation(0)
    local_check()
    check_head_detector()
    probes = harmonic_ports(car)
    active_source = 0
    neighbour = next(v for v in car.neighbours[0] if car.degree[v] == car.degree[0])
    ports = list(range(car.n))
    controls = [x for x in ports if x not in (0, neighbour)] + [0, neighbour]
    backgrounds = [('ascending', ports),
                   ('reversed_cell_labels', sorted(ports, key=lambda p: (p//12, (-p) % 12))),
                   ('native_two_state_control', controls)]
    arms = []
    for index, (name, order) in enumerate(backgrounds):
        h, sites, sector_info = sector(car, order)
        assert sector_info['omitted_states'] == 0
        active, p, counts = density_records(h, sites, 300522 + index)
        used = probes if index < 2 else [(None, None, np.eye(1, car.n, k=active_source).ravel())]
        rows = []
        for target, eig, weights in used:
            w = weights[active]
            exact = p @ w
            measured = counts @ w / SHOTS
            variance = p @ w**2 - exact**2
            se = float(np.sqrt(max(variance.max(), 0.)/ SHOTS))
            rows.append(dict(target=target, **{'lambda': None if eig is None else round(eig, 8)},
                             source_projection_visible=bool(np.linalg.norm(weights) > 0),
                             exact_fit=single_frequency(exact, se), count_fit=single_frequency(measured, se)))
        evidence = dict(orientation=name, reachable_states=len(sites), occupied_port_labels=active.tolist(),
                        shots_per_exposure=SHOTS, histogram=counts.tolist(), modes=rows)
        if index < 2:
            evidence['dispersion_test'] = shell(rows)
        else:
            assert len(sites) == 2
            expected = 2*abs(h[0, 1])
            assert all(rows[0][k]['resolved'] and abs(rows[0][k]['omega']/expected - 1) < .01
                       for k in ('exact_fit', 'count_fit'))
            evidence.update(control_frequency=float(expected), source_port=0, destination_port=neighbour)
        arms.append(evidence)
    return dict(exposure_spacing='1/24', exposures=769, last_exposure=32, final_training_exposure=16,
                probability_dictionary='load occupation from incident arrowhead change', boards=arms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write the complete histogram to a fresh receipt')
    output = parser.parse_args().output
    data = run()
    if output is not None:
        with output.open('x') as f:
            json.dump(data, f, allow_nan=False, separators=(',', ':'))
            f.write('\n')
    elif data != json.loads((HERE/'spectrum.json').read_text()):
        raise SystemExit('Density spectroscopy differs from the saved receipt')
    for arm in data['boards']:
        print(arm['orientation'], arm.get('dispersion_test', arm.get('control_frequency')),
              [(row['count_fit']['omega'], row['count_fit']['holdout_r2']) for row in arm['modes']])
