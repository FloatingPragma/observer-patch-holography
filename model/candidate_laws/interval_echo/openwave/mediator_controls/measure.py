"""Release a scalar load prepared on one pentagon and read the opposite pentagon."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import cell
from model.executable_model.symmetry import rotations
from model.experiments.interval_echo_message import headline as message_headline
from model.experiments.interval_echo_space.shell_source import local_source, prepare

HERE = Path(__file__).resolve().parent
EXPOSURES = (0, .25, .5, 1, 2)
SHOTS = 8192
TOL = 1e-11


def require(condition, message):
    if not condition:
        raise ValueError(message)


def values(array):
    return np.round(array, 12).tolist()


def layout(car):
    distance = shortest_path(car.adjacency(), unweighted=True, directed=False)
    far = distance.argmax(axis=1)
    require(np.all((distance == 3).sum(axis=1) == 1), 'ambiguous antipode')
    north, south = 0, int(far[0])
    first = min(car.neighbours[north])
    step = min(g for g in rotations() if g[north] == north and g[south] == south
               and g[first] in car.neighbours[first])
    source = [first]
    for _ in range(4):
        source.append(step[source[-1]])
    receiver = [int(far[p]) for p in source]
    require(set(source) == set(car.neighbours[north]) and set(receiver) == set(car.neighbours[south])
            and not set(source).intersection(receiver), 'rings must be disjoint and complete')
    for ring in (source, receiver):
        require([step[p] for p in ring] == ring[1:] + ring[:1], 'inconsistent axis orientation')
    return north, south, step, source, receiver


def preparations(car, source):
    # Only the fixed launch port connects the preparation tree to the rest of the cell.
    # Branches outside the source ring carry zero target mass and receive no splitters.
    ring = set(source)
    root = source[0]
    active = [e for e in car.seams if len(ring.intersection(e)) != 1 or root in e]
    board = SimpleNamespace(car=car, active=active)
    azimuth = 2 * np.pi * np.arange(5) / 5
    profiles = {f'circular_{m:+d}': np.exp(1j*m*azimuth)/np.sqrt(5) for m in (1, -1, 2, -2)}
    for m in (1, 2):
        for shift in range(5):
            profiles[f'standing_{m}_{shift}'] = np.sqrt(2/5)*np.cos(m*(azimuth-azimuth[shift]))
    states, pulses = {}, {}
    for label, profile in profiles.items():
        target = np.zeros(car.n, dtype=complex)
        target[source] = profile
        states[label], pulses[label] = prepare(board, target, root=root)
        require(pulses[label]['splitters'] == 4 and pulses[label]['routing_swaps'] == 0,
                'source acted outside its pentagon')
    return states, pulses


def position(state, receiver):
    population = abs(state)**2
    quiet = np.delete(population, receiver).sum()
    p = np.r_[population[receiver], quiet]
    require(abs(p.sum()-1) < TOL, 'incomplete receiver')
    return p


def phase_read(state, receiver, splitter):
    analyzed = state.copy()
    seam = receiver[:2]
    analyzed[seam] = splitter @ analyzed[seam]
    return position(analyzed, receiver)


def comparison(first, second):
    average = (first+second)/2
    info = 0.
    for row in (first, second):
        occupied = row > 0
        info += .5 * (row[occupied] @ np.log2(row[occupied]/average[occupied]))
    tv = abs(first-second).sum()/2
    return {'total_variation': values(tv), 'optimal_accuracy': values((1+tv)/2),
            'mutual_information_bits': values(max(0., info))}


def sample(first, second, stream):
    # Resolve roundoff in analytically tied positions using the frozen tie rule.
    decision = (second-first > TOL).astype(int)
    counts = []
    for label, row in enumerate((first, second)):
        rng = np.random.default_rng(np.random.SeedSequence([20261123, stream, label]))
        counts.append(rng.multinomial(SHOTS, row/row.sum()).tolist())
    correct = int(np.asarray(counts)[decision, np.arange(len(decision))].sum())
    return {'counts_by_label': counts, 'decision_by_outcome': decision.tolist(),
            'correct': correct, 'total': 2*SHOTS}


def inspect_modes(states, propagator, receiver, splitter, exposure):
    released = {label: propagator @ state for label, state in states.items()}
    records = {label: position(state, receiver) for label, state in released.items()}
    modes = {}
    for m in (1, 2):
        labels = [f'circular_{sign*m:+d}' for sign in (1, -1)]
        plain = [records[label] for label in labels]
        analyzed = [phase_read(released[label], receiver, splitter) for label in labels]
        standing = [records[f'standing_{m}_{shift}'] for shift in range(5)]
        require(np.max(abs(plain[0]-plain[1])) < TOL, 'plain position resolved circular sign')
        for row in plain:
            require(np.ptp(row[:5]) < TOL, 'angularly nonuniform circular intensity')
        for shift, row in enumerate(standing):
            require(np.max(abs(np.roll(standing[0][:5], shift)-row[:5])) < TOL and
                    abs(row[-1]-standing[0][-1]) < TOL, 'standing orientation did not propagate covariantly')
        density = abs(np.fft.fft(standing[0][:5]))
        allowed = {0, (2*m) % 5, (-2*m) % 5}
        forbidden = max(density[q] for q in range(1, 5) if q not in allowed)
        require(forbidden < TOL and (exposure == 0 or density[(2*m) % 5] > TOL),
                'unexpected angular density')
        comparisons = {'standing_orientation_0_vs_1': standing[:2], 'circular_sign_plain': plain,
                       'circular_sign_after_seam_exchange': analyzed}
        summary = {name: comparison(*pair) for name, pair in comparisons.items()}
        if exposure == 1:
            for index, (name, pair) in enumerate(comparisons.items()):
                summary[name]['sample'] = sample(*pair, stream=10*m+index)
        modes[str(m)] = {'plain_circular_plus': values(plain[0]), 'plain_circular_minus': values(plain[1]),
                         'standing_0': values(standing[0]), 'standing_1': values(standing[1]),
                         'standing_0_density_dft_magnitudes': values(density),
                         'phase_analyzed_plus': values(analyzed[0]), 'phase_analyzed_minus': values(analyzed[1]),
                         'comparisons': summary}
    return modes


def separate_message_replay():
    folder = Path(message_headline.__file__).parent
    protocol = json.loads((folder/'protocol.json').read_text())
    receipt = json.loads((folder/'receipt.json').read_text())
    error, all_windows = message_headline.check(protocol, receipt)
    endpoints = []
    for station in receipt['stations']:
        for window in station['windows']:
            if window['exposure'] in (protocol['early'], protocol['late']):
                endpoints.append({'cells': station['cells'], 'exposure': window['exposure'],
                                  'trials': window['native']['total'], 'native_correct': window['native']['correct'],
                                  'no_encoding_correct': window['no_encoding']['correct']})
    return {'checked_windows': len(all_windows), 'max_abs_probability_replay_error': float(error),
            'endpoints': endpoints}


def run():
    car = cell()
    north, south, rotation, source, receiver = layout(car)
    states, pulses = preparations(car, source)
    h = candidate.generator(car).toarray()
    require(candidate.basis(car) == [('port', p) for p in range(12)], 'expected scalar position carrier')
    perm = np.eye(car.n)[:, rotation]
    require(np.max(abs(h @ perm-perm @ h)) < TOL, 'native generator breaks chosen rotation')
    reflection = list(range(car.n))
    for ring in (source, receiver):
        for j, p in enumerate(ring):
            reflection[p] = ring[-j]
    mirror = np.eye(car.n)[:, reflection]
    require(np.array_equal(h @ mirror, mirror @ h), 'native generator breaks the ring reflection')
    for m in (1, -1, 2, -2):
        state = states[f'circular_{m:+d}']
        require(np.max(abs(perm @ state-np.exp(-2j*np.pi*m/5)*state)) < TOL, 'wrong rotation character')
        require(np.max(abs(mirror @ state-states[f'circular_{-m:+d}'])) < TOL,
                'circular sign sources must be mirror images')
    require(receiver[1] in car.neighbours[receiver[0]], 'readout splitter is not a native seam')
    _, blank = local_source(1, candidate.KAPPA)
    splitter = expm(-1j*np.pi/2*blank)
    windows = []
    no_hop_hits = 0.
    for exposure in EXPOSURES:
        release = expm(-1j*exposure*h)
        parked = np.exp(-1j*exposure*np.diag(h))
        for state in states.values():
            no_hop_hits = max(no_hop_hits, position(parked*state, receiver)[:5].sum())
        windows.append({'exposure': exposure, 'orders': inspect_modes(states, release, receiver, splitter, exposure)})
    require(no_hop_hits == 0, 'no-hop source reached receiver')
    return {'law': candidate.NAME, 'sector': 'one unit load, blank cell, one position amplitude per port',
            'north': north, 'south': south, 'proper_rotation': list(rotation),
            'commuting_ring_reflection': reflection,
            'northern_source': source, 'southern_outcomes_before_quiet': receiver,
            'source_pulses': pulses, 'collapsed_outcome': 'quiet comprises both poles and the entire northern ring',
            'no_hop_max_receiver_probability': values(no_hop_hits), 'windows': windows,
            'separate_intercell_message_replay': separate_message_replay()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'results.json')
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False)+'\n')
