"""Read a summed kinetic response without selecting one spectral line."""

import argparse
import json
from pathlib import Path

from scipy.linalg import expm

import networkx as nx
import numpy as np

from .. import dispersion, measure as flight

HERE = Path(__file__).resolve().parent
SHOTS = 2**58
TAUS = (4096, 8192)
SEED = 719193
ALPHA = .01
# 3 population banks, 3*6*2 native, 3*6*2 half-rate,
# 2*6*2 scalar, and 1 which-arm bank; each has <=20 bounded scores.
LOG_PENALTY = np.log(2 * 200 * 20 / ALPHA)


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def sample(raw, rng):
    raw = np.asarray(raw, float)
    check(np.isfinite(raw).all() and np.min(raw) >= 0, 'invalid raw event law')
    correction = abs(float(raw.sum()) - 1)
    check(correction < 3e-12, 'event law not normalized')
    probabilities = raw / raw.sum()
    return dict(raw_probabilities=raw.tolist(), normalization_correction=correction,
                counts=rng.multinomial(SHOTS, probabilities).tolist())


def spatial(board):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    edges = np.array(sorted(tuple(sorted(e)) for e in graph.edges()), int)
    profiles = []
    for shell, (eigenvalue, projector) in enumerate(dispersion.geometry(board), 1):
        for center in (0, 10):
            f = projector[:, center] / max(abs(projector[:, center]))
            profiles.append(dict(shell=shell, center=center, eigenvalue=eigenvalue, values=f.tolist()))
    return edges, profiles


def population(psi, x, y, cells):
    return (np.bincount(x // 12, weights=abs(psi)**2, minlength=cells)
            + np.bincount(y // 12, weights=abs(psi)**2, minlength=cells)) / 2


def population_record(p, admission, amplitude_errors, rng):
    record = sample(np.r_[admission * p, 1 - admission], rng)
    n = sum(record['counts'][:-1])
    check(n > 0, 'no admitted population records')
    seen = np.array(record['counts'][:-1]) / n
    sampling = float(np.sqrt(LOG_PENALTY / (2 * n)))
    lower = np.maximum(0, np.sqrt(np.maximum(0, seen - sampling)) - amplitude_errors)**2
    upper = np.minimum(1, (np.sqrt(np.minimum(1, seen + sampling)) + amplitude_errors)**2)
    record.update(conditional_owner_population=seen.tolist(), owner_probability_bounds=np.c_[lower, upper].tolist(),
                  amplitude_systematic_bounds=amplitude_errors.tolist(), sampling_halfwidth=sampling,
                  reader='Read both native port addresses and choose excess/vacancy with an independent fair classical bit')
    return record


def dirichlet_interval(record, edges, f):
    bounds = np.array(record['owner_probability_bounds'])
    squared = (f[edges[:, 0]] - f[edges[:, 1]])**2
    lo, hi = [float(np.sum(np.sqrt(bounds[edges[:, 0], i] * bounds[edges[:, 1], i]) * squared))
              for i in (0, 1)]
    check(lo > 0, 'population does not calibrate a kinetic metric')
    return lo, hi


def fringe(energy, compressed, rest_g, tau):
    gap = energy - energy[0]
    k = compressed
    a = np.exp(-1j * gap * tau) * (k @ rest_g)
    b = k @ (np.exp(-1j * gap[:len(rest_g)] * tau) * rest_g)
    # Coherent reference recombination: +/- are actual orthogonal Y exits.
    plus = float(np.linalg.norm((a - 1j * b) / 2)**2)
    minus = float(np.linalg.norm((a + 1j * b) / 2)**2)
    c = k.conj().T @ (np.exp(1j * gap * tau)[:, None] * k)
    c = c * np.exp(-1j * gap[:len(rest_g)] * tau)[None, :]
    operator = (c - c.conj().T) / (2j * tau)
    return plus, minus, operator


def response_record(energy, k, g, full_grating, outside, admission, tau, population_data, edges, f,
                    coefficient_error, operator_error, rng, which_arm=False):
    plus, minus, operator = fringe(energy, k, g, tau)
    gap = energy - energy[0]
    moved = np.exp(-1j * gap[:len(g)] * tau) * g
    grating_success = float((np.vdot(g, full_grating @ g).real
                              + np.vdot(moved, full_grating @ moved).real) / 2)
    in_band = plus + minus
    band_miss = float((np.linalg.norm(outside @ g)**2 + np.linalg.norm(outside @ moved)**2) / 2)
    check(abs(grating_success - in_band - band_miss) < 3e-12, 'final energy-filter decomposition fails')
    # Both retained branches are squared norms; no signed probability is clipped.
    grating_success = in_band + band_miss
    if which_arm:
        plus = minus = in_band / 2
    record = sample(np.r_[admission * np.array([plus, minus, grating_success - in_band,
                                                1 - grating_success]), 1 - admission], rng)
    n = sum(record['counts'][:4])
    observed = (record['counts'][0] - record['counts'][1]) / n / tau
    sampling = float(np.sqrt(2 * LOG_PENALTY / n) / tau)
    spread = float(np.ptp(np.linalg.eigvalsh(operator)))
    rest_width = float(gap[len(g) - 1])
    systematic = operator_error + spread * coefficient_error + rest_width
    bandwidth = float(gap[-1])
    sinc = float(np.sinc(bandwidth * tau / np.pi))
    check(bandwidth * tau < np.pi / 2 and sinc > 0, 'sine meter outside bounded bandwidth')
    moment = [max(0., observed - sampling - systematic),
              (observed + sampling + systematic) / sinc]
    d = dirichlet_interval(population_data, edges, f)
    record.update(tau=tau, sine_moment=observed, sampling_halfwidth=sampling,
                  numerical_systematic_bound=systematic,
                  operator_approximation_bound=operator_error, coefficient_error=coefficient_error,
                  response_operator_width=spread, rest_energy_width=rest_width,
                  band_energy_width=bandwidth, sinc_lower=sinc, first_moment_interval=moment,
                  population_dirichlet_interval=list(d),
                  kinetic_coefficient_interval=[moment[0] / d[1], moment[1] / d[0]],
                  category_order='Y plus; Y minus; grating passed but final band failed; grating failed; initial rest-filter failed')
    return record


def build_source(energies, basis, psi, g, admission, x, y, edges, profiles, pop_errors,
                 coefficient_error, operator_errors, rng):
    p = population(psi, x, y, 20)
    pop = population_record(p, admission, pop_errors, rng)
    outputs = {}
    for law, rate in (('native', 1), ('half_rate', .5)):
        rows = []
        for profile in profiles:
            f = np.array(profile['values'])
            value = (f[x // 12] + f[y // 12]) / 2
            k = basis.T @ (value[:, None] * basis[:, :len(g)])
            full = basis[:, :len(g)].T @ (value[:, None]**2 * basis[:, :len(g)])
            outside = value[:, None] * basis[:, :len(g)] - basis @ k
            records = [response_record(rate * energies, k, g, full, outside, admission, t, pop, edges, f,
                                       coefficient_error, rate * operator_errors, rng) for t in TAUS]
            rows.append(dict(shell=profile['shell'], center=profile['center'], records=records))
        outputs[law] = rows
    # No source/propagator changes; decohere the timing arms before recombination.
    f = np.array(profiles[0]['values'])
    value = (f[x // 12] + f[y // 12]) / 2
    k = basis.T @ (value[:, None] * basis[:, :len(g)])
    full = basis[:, :len(g)].T @ (value[:, None]**2 * basis[:, :len(g)])
    outside = value[:, None] * basis[:, :len(g)] - basis @ k
    which = response_record(energies, k, g, full, outside, admission, TAUS[-1], pop, edges, f,
                            coefficient_error, operator_errors, rng, which_arm=True)
    return dict(admission=admission, population=pop, response=outputs, which_arm_control=which)


def intersect(records):
    values = [r['kinetic_coefficient_interval'] for r in records]
    lo, hi = max(v[0] for v in values), min(v[1] for v in values)
    check(lo <= hi, 'two clocks disagree on the kinetic interval')
    return [lo, hi]


def decision(calibration, held):
    # A declared 5% effective-law tolerance is added beyond all interval errors.
    lo, hi = calibration[0] * .95, calibration[1] * 1.05
    return 'reject' if held[1] < lo or held[0] > hi else 'compatible'


def inference(sources, calibration=None):
    if calibration is None:
        calibration = intersect(sources[0]['response']['native'][0]['records'])
    rows = []
    for source in sources:
        for law, profiles in source['response'].items():
            for i, profile in enumerate(profiles):
                estimate = intersect(profile['records'])
                rows.append(dict(source=source['label'], law=law, shell=profile['shell'], center=profile['center'],
                                 kinetic_coefficient_interval=estimate,
                                 outcome='calibration' if source is sources[0] and law == 'native' and i == 0
                                 else decision(calibration, estimate)))
    return dict(calibration=calibration, operational_mass_interval=[1 / (2 * calibration[1]), 1 / (2 * calibration[0])],
                relative_law_tolerance=.05, held=rows)


def clock_standard(rng):
    # Actual sharp10/01 native seam, including its full contrast diagonal.
    h = np.array([[8.5, -.5], [-.5, 8.5]])
    rows = []
    for time in (.25, .5, .75):
        state = expm(-1j * time * h) @ np.array([1., 0.])
        record = sample(abs(state)**2, rng)
        seen = record['counts'][1] / SHOTS
        err = float(np.sqrt(LOG_PENALTY / (2 * SHOTS)))
        phase = [float(2 * np.arcsin(np.sqrt(max(0., seen - err)))),
                 float(2 * np.arcsin(np.sqrt(min(1., seen + err))))]
        record.update(commanded_appointment=time, native_gap_interval=[a / time for a in phase],
                      sharp_words=['10', '01'], compensated_input=[1, 0, 2], sampling_halfwidth=err)
        rows.append(record)
    calibration = rows[0]['native_gap_interval']
    check(calibration[0] <= 1 <= calibration[1], 'independent dimer gap calibration fails')
    check(all(max(calibration[0], r['native_gap_interval'][0]) <= min(calibration[1], r['native_gap_interval'][1])
              for r in rows[1:]), 'held dimer clock calibration fails')
    return dict(records=rows, gap_calibration=calibration,
                unit='one radian of the independently calibrated native odd-load dimer phase',
                scope='Native reference response measured; protected scheduling and counted long-time repetitions supplied')


def run():
    board = flight.federation(0)
    x, y, _, _, h, cost = flight.defect.generators(board)
    original = json.loads((flight.HERE / 'results.json').read_text())
    reduced = flight.PairReduction(h, cost)
    e, basis, residual, _ = reduced.lowest(h, 24)
    cut = (e[1] + e[2]) / 2
    check(reduced.inertia(cut)['upper_lines_below'] == 2 and reduced.inertia(flight.CUT)['upper_lines_below'] == 24,
          'moving/rest band count differs')
    band_error = max(original['cluster_projector_error'], residual / (flight.CUT - e[-1])) + 1e-12
    rest_error = residual / (min(cut - e[1], flight.CUT - e[-1]) - residual) + 1e-12
    operator_error = float(8 * np.ptp(e) * (rest_error + band_error) + 4 * residual + 240 * band_error**2)
    edges, profiles = spatial(board)
    rng = np.random.default_rng(SEED)
    native = []
    for launch in original['flights']:
        if launch['admission'] < .001:
            continue
        index = flight.defect.ordered_index(*launch['ports'], board.n)
        g = basis[index, :2]
        admission = float(g @ g)
        g = g / np.sqrt(admission)
        psi = basis[:, :2] @ g
        coefficient_error = float(min(2., 2 * rest_error / np.sqrt(admission)))
        pop_errors = []
        for owner in range(board.cells):
            c = ((x // 12 == owner).astype(float) + (y // 12 == owner).astype(float)) / 2
            local = basis[:, :2].T @ (c[:, None] * basis[:, :2])
            pop_errors.append(2 * rest_error + np.sqrt(max(np.linalg.eigvalsh(local))) * coefficient_error)
        result = build_source(e, basis, psi, g, admission, x, y, edges, profiles,
                              np.array(pop_errors), coefficient_error, operator_error, rng)
        result.update(label=f"owner_{launch['owner']}", initial_ports=launch['ports'])
        native.append(result)
    # Scalar references permit inhomogeneity and use the identical clocks and f readers.
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    lap = nx.laplacian_matrix(graph, nodelist=range(20)).toarray()
    scalar = []
    scalar_A = 2.**-19
    for name, envelope in (('uniform', np.ones(20)), ('inhomogeneous', 1 + .2 * np.array(profiles[0]['values']))):
        q = envelope / np.linalg.norm(envelope)
        potential = -scalar_A * (lap @ q) / q
        hs = scalar_A * lap + np.diag(potential)
        es, vs = np.linalg.eigh(hs)
        check(abs(es[0]) < 1e-15 and abs(abs(vs[:, 0] @ q) - 1) < 1e-12, 'scalar ground construction fails')
        # Pair-shaped address interface with x=y is only this declared reader control.
        result = build_source(es, vs, vs[:, 0], np.ones(1), 1., np.arange(20)*12, np.arange(20)*12,
                              edges, profiles, np.zeros(20), 0., 2e-18, rng)
        result.update(label=name, declared_kinetic_coefficient=scalar_A, supplied_potential=potential.tolist())
        scalar.append(result)
    native_inference = inference(native)
    scalar_inference = inference(scalar)
    check(any(r['outcome'] == 'reject' and r['law'] == 'native' for r in native_inference['held']),
          'claimed scalar-inertia obstruction not observed')
    check(all(r['outcome'] != 'reject' for r in scalar_inference['held'] if r['law'] == 'native'),
          'scalar positive reference rejected')
    clock = clock_standard(rng)
    return dict(shots=SHOTS, seed=SEED, taus=list(TAUS), clock_standard=clock, population_count_coin='fair excess/vacancy selection',
                owner_edges=edges.tolist(), profiles=profiles, residual=float(residual),
                band_projector_error=band_error, rest_projector_error=rest_error,
                native_operator_error=operator_error, native_sources=native, scalar_controls=scalar,
                native_inference=native_inference, scalar_inference=scalar_inference,
                scope='Rejects this single positive scalar ground-envelope kinetic coefficient; not a no-go for a matrix, nonlocal or momentum-resolved composite mass')


def compare(a, b, path='sum_rule_inertia'):
    if isinstance(b, dict):
        check(isinstance(a, dict) and a.keys() == b.keys(), path)
        for k in b:
            compare(a[k], b[k], path + '.' + k)
    elif isinstance(b, list):
        check(isinstance(a, list) and len(a) == len(b), path)
        for i, (c, d) in enumerate(zip(a, b, strict=True)):
            compare(c, d, path + f'[{i}]')
    elif isinstance(b, float):
        check(np.isfinite(a) and np.isclose(a, b, rtol=2e-8, atol=2e-12), path)
    else:
        check(type(a) is type(b) and a == b, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    encoded = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        with args.output.open('x') as out:
            out.write(encoded)
    else:
        compare(json.loads(encoded), json.loads((HERE / 'results.json').read_text()))
        print('Verified population-corrected summed inertia test and scalar/changed-law controls.')


if __name__ == '__main__':
    main()
