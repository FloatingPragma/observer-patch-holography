"""Local seam records of separated colour sources and domain walls in a prepared medium."""

import argparse
from itertools import permutations
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law

from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier

from .spectrum import physical_frame, spectrum

HERE = Path(__file__).resolve().parent
SEED = 220126
COLOUR_SHOTS = 2**27
PHASE_SHOTS = 2**20
KAPPAS = (int(candidate.KAPPA), 0, -int(candidate.KAPPA))


def skeleton(board):
    graph = nx.Graph()
    ports = {}
    for i, j in board.external:
        a, b = i // 12, j // 12
        graph.add_edge(a, b)
        ports[a, b], ports[b, a] = i, j
    return graph, ports


def aligned(band, port):
    axis = band[port] / np.linalg.norm(band[port])
    middle = np.eye(3)[np.argmin(abs(axis))]
    middle -= axis * (middle @ axis)
    middle /= np.linalg.norm(middle)
    return band @ np.column_stack((axis, middle, np.cross(axis, middle)))


def colour_states():
    epsilon = np.zeros((3, 3, 3))
    for p in permutations(range(3)):
        inversions = sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
        epsilon[p] = (-1)**inversions / np.sqrt(6)
    return {'singlet': epsilon, 'octet_word': np.roll(epsilon, 1, axis=0)}


def joint_record(rng, odd_probability, both_probability, shots):
    probs = np.array([odd_probability, both_probability,
                      1 - odd_probability - both_probability]).clip(0, 1).round(14)
    return list(map(int, rng.multinomial(shots, probs / probs.sum())[:2]))


def term_stats(counts, kappa, shots):
    odd, both = np.asarray(counts) / shots
    mean = odd - 2 * kappa * both
    var = (odd + (2 * kappa)**2 * both - mean**2) / shots
    return mean, max(0, var)


def colour_trial(board, graph, ports, band, rng):
    b, c = 0, min(graph[0])
    distances = nx.multi_source_dijkstra_path_length(graph, {b, c}, weight=None)
    choices = [min((a for a in graph if distances[a] == d),
                   key=lambda a: (d == 1 and not graph.has_edge(a, b), a))
               for d in range(1, max(distances.values()) + 1)]
    touching = choices[0]
    assert graph.has_edge(touching, b) and not graph.has_edge(touching, c)
    fa = aligned(band, ports[touching, b] % 12)
    fb = aligned(band, ports[b, touching] % 12)
    states = colour_states()
    result = []
    for d, a in enumerate(choices, 1):
        cut = sorted(graph[a])
        row = dict(distance=d, owner=a, seams=[[ports[a, n], ports[n, a]] for n in cut],
                   preparations={})
        for label, state in states.items():
            # These are independent copies read in the prescribed colour frame.
            # Source flags are local total occupancies, internal mode populations,
            # pair antisymmetry and the joint determinant projector.
            source_record = {}
            for j, name in enumerate(('A', 'B', 'C')):
                prob = np.moveaxis(state, j, 0).reshape(3, -1)
                source_record[name] = list(map(int, rng.multinomial(
                    COLOUR_SHOTS, (prob**2).sum(axis=1))))
            antisym = np.linalg.norm((state - state.swapaxes(1, 2)) / 2)**2
            singlet = abs(np.vdot(states['singlet'], state))**2
            source_record['BC_antisymmetric'] = int(rng.binomial(COLOUR_SHOTS, np.clip(antisym, 0, 1)))
            source_record['joint_determinant'] = int(rng.binomial(COLOUR_SHOTS, np.clip(singlet, 0, 1)))
            measurements, exact_joint = [], 0.
            for neighbour in cut:
                u = fa[ports[a, neighbour] % 12]
                p = float(u @ u / 3)
                q, both = 0., 0.
                if neighbour in (b, c):
                    frame = fb if neighbour == b else band
                    v = frame[ports[neighbour, a] % 12]
                    q = float(v @ v / 3)
                    indices = 'i,j,ijk->k' if neighbour == b else 'i,k,ijk->j'
                    both = float(np.linalg.norm(np.einsum(indices, u, v, state))**2)
                exact_joint += both
                arms = []
                for left, right, overlap in ((p, q, both), (p, 0., 0.), (0., q, 0.), (0., 0., 0.)):
                    arms.append(joint_record(rng, (left+right)/2-overlap, overlap, COLOUR_SHOTS))
                measurements.append(arms)
            energies = {}
            for k in KAPPAS:
                mean, variance = 0., 0.
                for four_arms in measurements:
                    for arm, sign in zip(four_arms, (1, -1, -1, 1), strict=True):
                        val, err = term_stats(arm, k, COLOUR_SHOTS)
                        mean += sign * val
                        variance += err
                expected = -(2*k+1)*exact_joint
                assert abs(mean - expected) < 6*np.sqrt(variance)
                energies[str(k)] = dict(exact=expected, measured=mean, se=np.sqrt(variance))
            row['preparations'][label] = dict(reader=source_record,
                                              odd_both_by_seam_and_arm=measurements,
                                              connected_energy=energies)
        row['fixed_source_spectrum'] = spectrum(graph, ports, (a, b, c), (fa, fb, band),
                                                states['singlet'], KAPPAS[0])
        result.append(row)
    forecasts = {}
    for label in states:
        values = [r['preparations'][label]['connected_energy'][str(KAPPAS[0])] for r in result]
        first, second = values[:2]
        slope = second['measured'] - first['measured']
        held = []
        for distance, item in enumerate(values[2:], 3):
            # The intercept and slope are correlated because both use the same two points.
            prediction = (2-distance)*first['measured'] + (distance-1)*second['measured']
            error = np.sqrt(((2-distance)*first['se'])**2 + ((distance-1)*second['se'])**2
                            + item['se']**2)
            held.append(dict(distance=distance, forecast=prediction, observed=item['measured'],
                             residual_z=(item['measured']-prediction)/error))
        forecasts[label] = dict(train_lengths=[1, 2], slope=slope,
                                slope_se=np.hypot(first['se'], second['se']), held=held)
    return dict(antitriplet_owners=[b, c], rows=result, contact_line_forecasts=forecasts)


def make_guide(graph, ports):
    paths = nx.single_source_shortest_path(graph, 0)
    path = paths[min(k for k in paths if len(paths[k]) == max(map(len, paths.values())))]
    aligned_ports, forward, spectators = [], [], []
    for i, owner in enumerate(path):
        off = min(set(graph[owner]) - set(path))
        aligned_ports.append(ports[owner, off])
        spectators.append(ports[off, owner])
        if 0 < i < len(path)-1:
            forward.append(ports[owner, path[i+1]])
    return path, aligned_ports, forward, spectators


def phase_read(occupied, aligned_ports, forward):
    phase = [0]
    assert {aligned_ports[0], aligned_ports[-1]} <= occupied
    for normal, flipped in zip(aligned_ports[1:-1], forward, strict=True):
        if (normal in occupied) == (flipped in occupied):
            raise ValueError('A guide owner has a missing or ambiguous marker')
        phase.append(int(flipped in occupied))
    phase.append(0)
    return [i for i, (left, right) in enumerate(zip(phase[:-1], phase[1:], strict=True))
            if left != right]


def sharp_read(occupied, edges, rng):
    # On a sharp unequal pair, swap-odd has Bernoulli(1/2); equal pairs
    # have definite parity. Independent term settings can be pooled exactly.
    mixed = sum((p in occupied) != (q in occupied) for p, q in edges)
    both = sum(p in occupied and q in occupied for p, q in edges)
    return dict(mixed_seams=mixed, both_seams=both,
                odd_count=int(rng.binomial(mixed * PHASE_SHOTS, .5)),
                both_count=int(both * PHASE_SHOTS))


def phase_energy(arm, reference, kappa):
    db = arm['both_seams'] - reference['both_seams']
    dm = arm['mixed_seams'] - reference['mixed_seams']
    measured = ((arm['odd_count']-reference['odd_count'])
                - 2*kappa*(arm['both_count']-reference['both_count'])) / PHASE_SHOTS
    return dict(exact=dm/2-2*kappa*db, measured=measured,
                se=np.sqrt((arm['mixed_seams']+reference['mixed_seams'])/(4*PHASE_SHOTS)))


def medium_trial(board, graph, ports, rng):
    path, normal, forward, spectators = make_guide(graph, ports)
    sensitive = set(normal+forward)
    edges = [e for e in board.seams if any(p in sensitive for p in e)]
    confs = {}
    for mode in ('connected', 'removed_spectators', 'outer_records_only'):
        base = set(normal) | (set(spectators) if mode != 'removed_spectators' else set())
        rows = []
        for m in range(1, min(12, len(forward))+1):
            changed = {0, m-1} if mode == 'outer_records_only' else set(range(m))
            state = (base - {normal[i+1] for i in changed}) | {forward[i] for i in changed}
            assert len(state) == len(base)
            record, reference = sharp_read(state, edges, rng), sharp_read(base, edges, rng)
            walls = phase_read(state, normal, forward)
            assert walls[0] == 0 and walls[-1] == m
            energies = {str(k): phase_energy(record, reference, k) for k in KAPPAS}
            for item in energies.values():
                assert abs(item['measured'] - item['exact']) < 6 * item['se']
            rows.append(dict(length=m, phase_boundaries=walls, trial=record, reference=reference,
                             energy=energies))
        confs[mode] = rows
    predictions = {}
    for mode in ('connected', 'removed_spectators'):
        rows = confs[mode]
        for k in KAPPAS:
            # Work relative to the unperturbed length zero is identically zero.
            values = [r['energy'][str(k)] for r in rows]
            first = values[:2]
            weights = np.array([1/v['se']**2 for v in first])
            x = np.array([1., 2.])
            scale = float(sum(weights*x*x))
            slope = float(sum(weights*x*[v['measured'] for v in first])/scale)
            held_z = [abs(v['measured']-m*slope)/np.sqrt(v['se']**2 + m*m/scale)
                      for m, v in enumerate(values[2:], 3)]
            assert max(held_z) < 6
            predictions[f'{mode}_{k}'] = dict(train_lengths=[1, 2], slope=slope, se=1/np.sqrt(scale),
                                              max_held_z=max(held_z), held_count=len(held_z))
    return dict(route=path, aligned_ports=normal, forward_ports=forward, spectators=spectators,
                sensitive_seam_count=len(edges), preparations=confs, predictions=predictions)


def pulse_controls(rng):
    exchange, contrast = law.local_term_parts()
    one, other, both = (np.eye(len(exchange))[:, law.SRC.index(pair)]
                        for pair in ((0, 1), (1, 0), (1, 1)))
    result = {}
    for k in KAPPAS:
        unitary = expm(-1j * np.pi/32 * (exchange + k*contrast))
        arms = {}
        for name, state in (('even', (one+other)/np.sqrt(2)), ('odd', (one-other)/np.sqrt(2)),
                            ('both', both)):
            out = unitary @ state
            overlap = np.vdot(state, out)
            prob = np.clip([(1+overlap.real)/2, (1+overlap.imag)/2,
                            np.linalg.norm(out[:len(law.SRC)])**2], 0, 1).round(14)
            arms[name] = dict(x_y_blank=list(map(int, rng.binomial(PHASE_SHOTS, prob))),
                              probabilities=list(prob))
        result[str(k)] = arms
    return result


def acquire():
    rng = np.random.default_rng(SEED)
    band = physical_frame(currents._service(carrier.cell(), False, 1.)[4][0])
    result = {}
    for level in range(3):
        board = carrier.federation(level)
        graph, ports = skeleton(board)
        result[str(level)] = dict(owners=board.cells,
                                  colour=colour_trial(board, graph, ports, band, rng),
                                  medium=medium_trial(board, graph, ports, rng))
    return currents._portable(dict(seed=SEED, colour_shots=COLOUR_SHOTS, phase_shots=PHASE_SHOTS,
                                  levels=result, pulse_controls=pulse_controls(rng)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='regenerate the frozen receipt')
    args = parser.parse_args()
    path = HERE / 'records.json'
    records = acquire()
    if args.write:
        path.write_text(json.dumps(records, indent=1, allow_nan=False) + '\n')
    elif records != json.loads(path.read_text()):
        raise ValueError('OW-22 record replay differs')
    print('OW-22 records verified on 20, 80, 320 owners' if not args.write else 'Wrote OW-22 records')
    from . import dressing
    dressing.replay(write=args.write)


if __name__ == '__main__':
    main()
