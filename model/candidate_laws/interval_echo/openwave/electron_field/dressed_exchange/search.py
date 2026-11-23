"""Return interference of two identical dressed charges on a native four-face gateway."""
import argparse
from collections import deque
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0, str(Path(__file__).resolve().parents[6]))
from model.executable_model import federation


LEFT = (14, 18, 15, 16, 23, 22)
RIGHT = (24, 31, 25, 29, 35, 34)
HEAD_COUNTS = (3, 1, 2, 0, 2, 1)


def edges_of(face):
    return list(zip(face, (*face[1:], face[0])))


class Fixture:
    def __init__(self, sign, head_at_left):
        self.sign = sign
        self.board = federation(0)
        faces = [tuple(fan[i] for i in (0, j, j+1)) for fan in (LEFT, RIGHT) for j in range(1, 5)]
        assert all(tuple(sorted(face)) in self.board.faces() for face in faces)
        self.bridge = (LEFT[0], RIGHT[0])
        assert self.bridge in self.board.external
        self.edges = sorted({tuple(sorted(e)) for face in faces for e in edges_of(face)} | {self.bridge})
        self.bits = {e: 1 << j for j, e in enumerate(self.edges)}
        self.mask = (1 << len(self.edges)) - 1
        self.rings = [self.orientation(face) for face in faces]
        words = []
        for k, fan in enumerate((LEFT, RIGHT)):
            local_edges = [e for e in self.edges if set(e).issubset(fan)]
            allowed = []
            for local in range(1 << len(local_edges)):
                counts = [0]*6
                word = 0
                for j, edge in enumerate(local_edges):
                    bit = local & (1 << j)
                    counts[fan.index(edge[0] if bit else edge[1])] += 1
                    if bit:
                        word |= self.bits[edge]
                if tuple(counts) == HEAD_COUNTS:
                    allowed.append(word)
            assert len(allowed) == 4
            local_order = sorted(tuple(sorted((fan.index(a), fan.index(b)))) for a,b in local_edges)
            def fan_coordinate(w):
                value = 0
                for j, (u, v) in enumerate(local_order):
                    e = tuple(sorted((fan[u], fan[v])))
                    if (fan[u] == e[0]) == bool(w & self.bits[e]):
                        value |= 1 << j
                return value
            base = min(allowed, key=fan_coordinate)
            first, last = self.rings[4*k][0], self.rings[4*k+3][0]
            square = [base ^ (u*first) ^ (v*last) for u, v in product((0, 1), repeat=2)]
            assert set(square) == set(allowed)
            words.append(square)
        bridge_word = self.bits[self.bridge] if head_at_left else 0
        self.seeds = [(self.bridge, (a | b | bridge_word) ^ (self.mask if sign < 0 else 0))
                      for a in words[0] for b in words[1]]
        assert len(set(self.seeds)) == 16
        # These are the two prepared buffer words in the formal four-contact exchange.
        if head_at_left:
            centre, peer, x, y, z = 24, 14, 35, 34, 29
        else:
            centre, peer, x, y, z = 14, 24, 23, 22, 16
        touched = [tuple(sorted((centre, v))) for v in (x, y, peer, z)]
        eligible = {tuple(e[0] if (sign > 0) == bool(w & self.bits[e]) else e[1]
                          for e in touched) for _, w in self.seeds}
        assert eligible == {(centre, y, peer, centre), (x, centre, peer, centre)}
        self.reserves = self.injection_routes(self.seeds[0][1])

    def injection_routes(self, seed):
        # Reverse an enabled outward material path in the sharp source configuration.
        # Each terminal inactive seam supplies a fixed opposite compensating load.
        chosen = []
        support = set(LEFT + RIGHT)
        for fan in (LEFT, RIGHT):
            queue = deque([(fan[0],)])
            while queue:
                path = queue.popleft()
                p = path[-1]
                destinations = sorted(set(self.board.neighbours[p]) - support - {x[-1] for x in chosen})
                if destinations:
                    chosen.append((*path, destinations[0]))
                    break
                for v in self.board.neighbours[p]:
                    e = tuple(sorted((p, v)))
                    if v in fan and v not in path and e in self.bits:
                        if (self.sign == 1) == ((e[0] == p) == bool(seed & self.bits[e])):
                            queue.append((*path, v))
            else:
                raise ValueError('Missing charge-compensation route')
        return chosen

    def orientation(self, face):
        e = edges_of(face)
        return (sum(self.bits[tuple(sorted(pair))] for pair in e),
                sum(self.bits[a, b] for a, b in e if a < b))

    def moves(self, occupants, word):
        for edge, bit in self.bits.items():
            a, b = edge
            if (a in occupants) == (b in occupants):
                continue
            source, dest = (a, b) if a in occupants else (b, a)
            if (self.sign == 1) != ((source == a) == bool(bit & word)):
                continue
            other = occupants[1 - occupants.index(source)]
            crossing = -1 if a < other < b else 1
            yield (tuple(sorted((other, dest))), word ^ bit), 4., crossing
        for mask, cw in self.rings:
            if word & mask in (cw, cw ^ mask):
                yield (occupants, word ^ mask), .5, 1


def six_pulse_exchange(model):
    """Enumerate shortest sharp-word exchanges then replay each fixed sequence on every code word."""
    steps = sorted(model.edges) + [("ring", k) for k in range(len(model.rings))]
    def apply(state, step):
        first, second, word = state
        if step[0] == "ring":
            mask, cw = model.rings[step[1]]
            return (first, second, word ^ mask) if word & mask in (cw, cw ^ mask) else state
        a, b = step
        if (a in (first, second)) == (b in (first, second)):
            return state
        src, dest = (a, b) if a in (first, second) else (b, a)
        bit = model.bits[step]
        if (model.sign == 1) != ((src == a) == bool(bit & word)):
            return state
        return (dest if first == src else first, dest if second == src else second, word ^ bit)
    a, b = model.bridge
    initial = a, b, model.seeds[0][1]
    final = b, a, initial[2]
    previous = {initial: []}
    front = {initial}
    for _ in range(6):
        new = {}
        for state in sorted(front):
            for k, step in enumerate(steps):
                target = apply(state, step)
                if target not in previous:
                    new.setdefault(target, []).append((state, k))
        previous.update(new)
        front = new
    assert final in front
    def sequences(target):
        if target == initial:
            yield ()
        for prior, step in previous[target]:
            for prefix in sequences(prior):
                yield (*prefix, step)
    records = []
    for sequence in sequences(final):
        outcomes = []
        for _, word in model.seeds:
            end = (a, b, word)
            for k in sequence:
                end = apply(end, steps[k])
            outcomes.append(end)
        exchanged = sum(end == (b, a, w) for end, (_, w) in zip(outcomes, model.seeds))
        direct = sum(end == (a, b, w) for end, (_, w) in zip(outcomes, model.seeds))
        records.append((exchanged, direct, sequence))
    exchanged, direct, best = max(records)
    # A second buffer at the same junction can accept whichever circulating head was blocked.
    # Try its two touches in all four relative orders. No control depends on a measurement.
    extensions = []
    for _, _, sequence in records:
        for k in sorted(set(sequence)):
            occurrences = [i for i, act in enumerate(sequence) if act == k]
            if len(occurrences) != 2 or steps[k] == model.bridge or steps[k][0] == "ring":
                continue
            for m, extra in enumerate(steps[:len(model.edges)]):
                if m in sequence or not set(steps[k]).intersection(extra):
                    continue
                for order in product((0, 1), repeat=2):
                    extended = []
                    for i, act in enumerate(sequence):
                        if act != k:
                            extended.append(act)
                        else:
                            extended.extend((m,k) if order[occurrences.index(i)] else (k,m))
                    matches = 0
                    for _, w in model.seeds:
                        end = (a, b, w)
                        for j in extended:
                            end = apply(end, steps[j])
                        matches += end == (b, a, w)
                    extensions.append((matches, tuple(extended)))
    best_count, extended = max(extensions)
    assert best_count == 16
    native, control = np.zeros((16,16)), np.zeros((16,16))
    word_ids = {w: i for i, (_,w) in enumerate(model.seeds)}
    for i, (_, w) in enumerate(model.seeds):
        end = (a, b, w)
        for k in extended:
            end = apply(end, steps[k])
        first, second, out_word = end
        native[word_ids[out_word], i] = 1
        control[word_ids[out_word], i] = 1 if first < second else -1
    square = np.array([[(-1)**u/2, (-1)**v/2] for u,v in product((0,1),repeat=2)])
    both = np.kron(square, square)
    np.testing.assert_array_equal(both.T @ native @ both, np.eye(4))
    np.testing.assert_array_equal(both.T @ control @ both, -np.eye(4))
    eight = {'native_doublet_return_matrix': (both.T @ native @ both).tolist(),
             'odd_permutation_doublet_return_matrix': (both.T @ control @ both).tolist(),
             'distinct_schedules_tested': len(set(s for _, s in extensions)),
             'max_sharp_words_preserved_and_exchanged': best_count,
             'different_schedules_exchanging_all_16': len(set(s for count,s in extensions if count == 16)),
             'best_schedule': [steps[k] for k in extended]}
    return {'distinct_shortest_schedules_for_first_sharp_input': len(records),
            'best_schedule': [steps[k] for k in best], 'exchanges_returning_same_heads_of_16': exchanged,
            'identity_returning_same_heads_of_16': direct, 'other_outputs_of_16': 16-exchanged-direct,
            'schedules_exchanging_all_16_inputs': sum(e == 16 for e, _, _ in records),
            'paired_buffer_eight_pulse_search': eight}


def simulate(sign, head_at_left):
    model = Fixture(sign, head_at_left)
    basis = list(model.seeds)
    ids = {state: i for i, state in enumerate(basis)}
    row, col, native, odd = [], [], [], []
    adjacency = []
    baseline = model.seeds[0][1]
    for i, (occupants, word) in enumerate(basis):
        delta = {v: 0 for v in LEFT + RIGHT}
        for (a, b), bit in model.bits.items():
            change = int(bool(word & bit)) - int(bool(baseline & bit))
            delta[a] += change
            delta[b] -= change
        assert all(delta[v] == sign*(int(v in occupants)-int(v in model.bridge)) for v in delta)
        diag = 8.*sum(6 - int(model.board.degree[v]) for v in occupants)
        neighbours = []
        for target, strength, parity in model.moves(occupants, word):
            j = ids.get(target)
            if j is None:
                j = len(basis)
                basis.append(target)
                ids[target] = j
            row.append(i)
            col.append(j)
            native.append(-strength)
            odd.append(-strength*parity)
            diag += strength
            neighbours.append((j, parity))
        adjacency.append(neighbours)
        row.append(i)
        col.append(i)
        native.append(diag)
        odd.append(diag)
    n = len(basis)
    parking = np.array([i for i, (occ, _) in enumerate(basis) if occ == model.bridge])
    assert set(parking) == set(range(16))
    square = np.array([[(-1)**u/2, (-1)**v/2] for u,v in product((0,1),repeat=2)])
    prepared = np.zeros((n, 4))
    prepared[:16] = np.kron(square, square)
    np.testing.assert_array_equal(prepared.T @ prepared, np.eye(4))
    cases = {'first_first': [1,0,0,0], 'second_second': [0,0,0,1],
             'symmetric_mixed': np.array([0,1,1,0])/np.sqrt(2),
             'antisymmetric_mixed': np.array([0,1,-1,0])/np.sqrt(2)}
    data, matrices = {}, []
    for key, weights in (('native',native), ('odd_worldline_control',odd)):
        h = csr_matrix((weights, (row, col)), shape=(n,n))
        assert (h-h.T).nnz == 0
        data[key] = []
        amp = []
        for k, state in enumerate(expm_multiply(-1j*h, prepared, start=0, stop=8, num=65)):
            np.testing.assert_allclose(state.conj().T @ state, np.eye(4), atol=1e-11)
            transfer = prepared.T @ state
            amp.append(transfer)
            by = {}
            for name, v in cases.items():
                gateway = float(np.linalg.norm(state[parking] @ v)**2)
                doublets = float(np.linalg.norm(transfer @ v)**2)
                by[name] = {'both_parked_doublets': doublets, 'bare_gateways_inclusive': gateway,
                            'gateways_other_field': gateway - doublets, 'other_positions': 1-gateway}
            data[key].append({'time': k/8, 'inputs': by})
        matrices.append(np.array(amp))
    # The virtual path lift keeps material origin only to check whether exchanges are reachable.
    # The observer reads the unlabelled charges and local parked modes above.
    distances = {(0, 1): 0}
    queue = deque([(0, 1)])
    while queue:
        i, parity = queue.popleft()
        for j, factor in adjacency[i]:
            target = j, parity*factor
            if target not in distances:
                distances[target] = distances[i, parity] + 1
                queue.append(target)
    exchange = (matrices[0] - matrices[1])/2
    direct = (matrices[0] + matrices[1])/2
    swap = np.eye(4)[[0,2,1,3]]
    exchange_success = np.linalg.svd(exchange, compute_uv=False)[:,-1]**2
    k = int(np.argmax(exchange_success))
    norm = float(np.sum(abs(exchange[k])**2))
    sector = {'best_worst_input_time': k/8, 'worst_input_odd_path_weight': float(exchange_success[k]),
              'mean_odd_path_weight': norm/4, 'swap_process_overlap_at_that_time':
              float(abs(np.trace(swap @ exchange[k]))**2/(4*norm)) if norm else None,
              'identity_process_overlap_at_that_time':
              float(abs(np.trace(exchange[k]))**2/(4*norm)) if norm else None}
    delta_read = {}
    for name in cases:
        p = np.array([entry['inputs'][name]['both_parked_doublets'] for entry in data['native']])
        f = np.array([entry['inputs'][name]['both_parked_doublets'] for entry in data['odd_worldline_control']])
        k = int(np.argmax(abs(p-f)))
        delta_read[name] = {'time': k/8, 'native': float(p[k]), 'control': float(f[k]),
                            'absolute_difference': float(abs(p[k]-f[k]))}
    return {'sign': sign, 'reference_bridge_head_at_left_before_conjugation': head_at_left,
            'compensation_paths': model.reserves, 'six_pulse_exchange': six_pulse_exchange(model), 'dimension': n, 'label_lift_dimension_from_one_seed': len(distances),
            'shortest_exchange_restoring_charges_and_full_sharp_word': distances.get((0,-1)),
            'max_operator_norm_odd_path_return_on_two_doublets': float(np.linalg.norm(exchange, 2, axis=(1,2)).max()),
            'max_operator_norm_even_path_return_on_two_doublets': float(np.linalg.norm(direct, 2, axis=(1,2)).max()),
            'largest_probability_differences': delta_read, 'odd_history_field_transfer': sector, 'readings': data}


def run():
    cases = [simulate(s, b) for b in (True, False) for s in (1, -1)]
    for a, b in ((cases[0], cases[1]), (cases[2], cases[3])):
        for key in a['readings']:
            for x, y in zip(a['readings'][key], b['readings'][key]):
                for name in x['inputs']:
                    np.testing.assert_allclose(list(x['inputs'][name].values()),
                                               list(y['inputs'][name].values()), atol=1e-11)
    return {'physical_fans':[LEFT, RIGHT], 'cases':cases}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    path = parser.parse_args().output
    result = run()
    with path.open('x', encoding='utf-8') as f:
        json.dump(result, f, allow_nan=False, indent=2)
        f.write('\n')
