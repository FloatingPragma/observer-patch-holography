"""Held joint-phase spectroscopy with conjugate native terminal contacts."""

import argparse
import base64
import json
from pathlib import Path
import zlib

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .experiment import HERE, local_check, sector
from .phase import terminal_gates
from .spectrum import TIMES, check_head_detector, harmonic_ports, shell, single_frequency


SHOTS = 16384
SETTINGS = ('sham', 'plus', 'minus', 'dephased_plus', 'dephased_minus')


def address_phases():
    law = candidate.native_law
    _, contrast = law.local_term_parts()
    inputs = np.eye(law.DIM)[:, [law.SRC.index((0, 0)), law.SRC.index((1, 0))]]
    phases = {}
    for angle in (np.pi/2, 3*np.pi/2):
        contact = expm(-1j * angle * contrast) @ inputs
        extracted = inputs.T @ contact
        if np.linalg.norm(contact - inputs @ np.diag([1, np.exp(-1j*angle)])) > 1e-11:
            raise ArithmeticError('degree-reference phase is not its declared native operation')
        phases[angle] = extracted[1, 1]
    return phases


def sample_probabilities(rng, p):
    if p.min() < -1e-10 or p.max() > 1 + 1e-10:
        raise ArithmeticError('terminal instrument gave an invalid probability')
    return np.asarray(rng.binomial(SHOTS, np.clip(p, 0, 1)), dtype='<u2')


def counts_payload(counts):
    raw = np.ascontiguousarray(counts, dtype='<u2').tobytes()
    packed = zlib.compress(raw, level=6)
    assert zlib.decompress(packed) == raw
    return dict(shape=list(counts.shape), dtype='<u2', compression='base64(zlib(array bytes in C order))',
                settings=list(SETTINGS), clicks=base64.b64encode(packed).decode('ascii'))


def link_channels(car, h, sites, state):
    """Sum only coherences with the same spectator coins and the opposite active coin."""
    edges = np.array(car.seams)
    ids = {tuple(e): i for i, e in enumerate(car.seams)}
    hc = h.tocoo()
    forward = sites[hc.row] > sites[hc.col]
    x, y = hc.row[forward], hc.col[forward]
    assert np.allclose(hc.data[forward], -candidate.KAPPA/2)
    edge_ids = [ids[int(sites[j]), int(sites[i])] for i, j in zip(x, y, strict=True)]
    # No one local state has two partners on the same physical seam.
    assert len(set(zip(x, edge_ids))) == len(x) == len(set(zip(y, edge_ids)))
    aggregate = sparse.csr_matrix((np.ones(len(x)), (edge_ids, np.arange(len(x)))),
                                  shape=(len(edges), len(x)))
    to_ports = sparse.csr_matrix((np.ones(len(sites)), (sites, np.arange(len(sites)))),
                                 shape=(car.n, len(sites)))
    density = abs(state.T)**2
    prob = (to_ports @ density).T
    left_right = (aggregate @ (state[:, x].conj()*state[:, y]).T).T
    px = (aggregate @ density[x]).T
    py = (aggregate @ density[y]).T
    derivative = 2 * (to_ports @ (state.T.conj() * (h @ state.T)).imag).T
    return dict(ordinary=prob[:, edges[:, 1]], x=px, y=py, product=left_right,
                density_derivative=derivative)


def interrogate(gates, address, weights, phase, sign, target, channels, seed):
    all_readings = channels['ordinary'][:, address]
    px, py = (channels[key][:, address] for key in ('x', 'y'))
    spectator = all_readings - px
    if spectator.min() < -1e-10:
        raise ArithmeticError('the active pair probabilities exceed the endpoint probability')
    z = channels['product'][:, address] * phase
    rng = np.random.default_rng(seed)
    counts = [sample_probabilities(rng, all_readings)]
    contrasts = []
    variances = []
    for dephased in (False, True):
        responses = []
        for setting in ('plus', 'minus'):
            a, b = gates[setting][0]
            p = spectator + abs(a)**2 * px + abs(b)**2 * py
            if not dephased:
                p = p + 2 * (a.conjugate() * b * z).real
            responses.append(p)
            counts.append(sample_probabilities(rng, p))
        plus, minus = responses
        actual = (plus-minus)*sign
        predicted = np.zeros_like(z.real) if dephased else target[:, address]
        if np.max(abs(actual - predicted)) > 1e-10:
            raise ArithmeticError('native plus/minus terminal response disagrees with the joint quadrature')
        contrasts.append(actual)
        variances.append((plus*(1-plus) + minus*(1-minus))/SHOTS)
    matrix = np.stack(counts, axis=-1)
    signed = sign * (matrix[:, :, 1].astype(float)-matrix[:, :, 2])/SHOTS
    dep = sign * (matrix[:, :, 3].astype(float)-matrix[:, :, 4])/SHOTS
    rows = []
    for spatial_target, eig, w in weights:
        norm = np.linalg.norm(w)
        if norm <= 1e-12:
            raise ValueError('blind modulation at this momentum')
        v = w/norm
        se = float(np.sqrt(max(float(np.max(variances[0] @ v**2)), 0.)))
        dephased_se = float(np.sqrt(max(float(np.max(variances[1] @ v**2)), 0.)))
        rows.append(dict(target=spatial_target, **{'lambda': None if eig is None else round(eig, 8)},
                         exact_fit=single_frequency(contrasts[0] @ v, se),
                         count_fit=single_frequency(signed @ v, se),
                         dephased_count_fit=single_frequency(dep @ v, dephased_se)))
    if any(r['dephased_count_fit']['resolved'] for r in rows):
        raise ArithmeticError('the dephased control reports a resolved frequency')
    return dict(read_link_ids=address.tolist(), records=counts_payload(matrix), modes=rows)


def quadrature_control(neighbour, h, sites, gates, phases):
    # Evolve to one quarter of a population oscillation from a sharply occupied source.
    separation = 2*abs(h[0, 1])
    src = np.array([1., 0.])
    state = expm_multiply(-1j * (np.pi / (2*separation)) * h, src)
    results = []
    for i, angle in enumerate((np.pi/2, 3*np.pi/2)):
        prepared = state.copy()
        prepared[sites == 0] *= phases[angle]
        # The real reader first applies its quarter-phase at a (the source here).
        operated = prepared.copy()
        operated[sites == 0] *= phases[np.pi/2]
        p = {}
        rng = np.random.default_rng(301843+i)
        for name in ('plus', 'minus'):
            a, b = gates[name][0]
            x, y = operated[sites == neighbour][0], operated[sites == 0][0]
            p[name] = float(abs(a*x+b*y)**2)
            p['dephased_'+name] = float(abs(a*x)**2 + abs(b*y)**2)
        clicks = {name: int(sample_probabilities(rng, np.array(prob))) for name,prob in p.items()}
        results.append(dict(preparation_angle=float(angle), probabilities={k: round(v, 10) for k,v in p.items()},
                            clicks=clicks, expected_contrast=round(p['plus']-p['minus'], 10),
                            observed_contrast=round((clicks['plus']-clicks['minus'])/SHOTS, 8)))
    assert min(abs(r['observed_contrast']) for r in results) > .9
    assert results[0]['observed_contrast'] * results[1]['observed_contrast'] < 0
    assert max(abs(r['probabilities']['dephased_plus'] - r['probabilities']['dephased_minus']) for r in results) < 1e-10
    return results


def run():
    car = federation(0)
    edges = np.array(car.seams)
    local_check()
    check_head_detector()
    gates, phases = terminal_gates(), address_phases()
    harmonics = harmonic_ports(car)
    a, b = edges.T
    read_i = np.arange(len(edges))
    read_r = np.flatnonzero((car.degree[a] == 5) | (car.degree[b] == 5))
    sign = np.where(car.degree[a[read_r]] == 5, 1, -1)
    phase = np.where(sign == 1, phases[np.pi/2], phases[np.pi/2].conjugate())
    native_neighbour = next(v for v in car.neighbours[0] if car.degree[v] == car.degree[0])
    ports = list(range(car.n))
    control = [p for p in ports if p not in (0, native_neighbour)] + [0, native_neighbour]
    orientations = [('ascending', ports),
                    ('reversed_cell_labels', sorted(ports, key=lambda p: (p//12, (-p) % 12))),
                    ('native_two_state_control', control)]
    arms = []
    for i, (name, order) in enumerate(orientations):
        h, sites, info = sector(car, order)
        if info['omitted_states']:
            raise ArithmeticError('the configuration sector is incomplete')
        source = np.zeros(len(sites))
        source[0] = 1
        states = expm_multiply(-1j*h, source, start=0, stop=32, num=len(TIMES))
        if np.max(abs(np.sum(abs(states)**2, axis=1)-1)) > 1e-10:
            raise ArithmeticError('free density matrix lost normalization')
        data = link_channels(car, h, sites, states)
        imag, real = -2*data['product'].imag, 2*data['product'].real
        arm = dict(orientation=name, states=len(sites))
        if i < 2:
            longitudinal = [(tar, eig, v[b]-v[a]) for tar, eig, v in harmonics]
            modulated = [(tar, eig, .5*(v[b[read_r]]+v[a[read_r]])) for tar, eig, v in harmonics]
            for (_, _, v), (_, _, w) in zip(harmonics, longitudinal, strict=True):
                if np.max(abs(data['density_derivative'] @ v - candidate.KAPPA/2 * (imag @ w))) > 1e-10:
                    raise ArithmeticError('the current violates the native continuity equation')
            arm['current'] = interrogate(gates, read_i, longitudinal, 1., 1, imag, data, 301810+2*i)
            arm['real_bond'] = interrogate(gates, read_r, modulated, phase, sign, real, data, 301811+2*i)
            for key in ('current', 'real_bond'):
                grade = shell(arm[key]['modes'])
                # The spatial weights of the real response modulate measured bonds.
                arm[key]['frequency_test'] = {k: (v.replace('density envelope', 'response envelope')
                                                   if isinstance(v, str) else v) for k,v in grade.items()}
        else:
            assert len(sites) == 2
            single_edge = np.array([car.seams.index((0, native_neighbour))])
            arm['current'] = interrogate(gates, single_edge, [(None,None,np.ones(1))],
                                         1., 1, imag, data, 301821)
            fit = arm['current']['modes'][0]
            assert all(fit[k]['resolved'] and abs(fit[k]['omega']/8-1) < .01 for k in ('exact_fit', 'count_fit'))
            arm['real_preparation_control'] = quadrature_control(native_neighbour, h, sites, gates, phases)
        arms.append(arm)
    return dict(boards=arms, exposures=len(TIMES), spacing='1/24', last_training_exposure=16,
                shots_per_link_contact=SHOTS, quiet='shots_per_link_contact minus clicks for each setting')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write to a fresh receipt')
    output = parser.parse_args().output
    result = run()
    if output:
        encoded = json.dumps(result, allow_nan=False, separators=(',', ':'))
        with output.open('x') as f:
            f.write(encoded + '\n')
    elif result != json.loads((HERE/'phase_spectrum.json').read_text()):
        raise SystemExit('Phase spectrum differs from the archived receipt')
    for arm in result['boards']:
        print(arm['orientation'])
        for key in ('current', 'real_bond'):
            if key in arm:
                print(key, [(r['count_fit']['omega'], r['count_fit']['holdout_r2']) for r in arm[key]['modes']])
