"""Controlled winding and free release of a dual-path parity state on physical seams."""
from __future__ import annotations

import argparse
from collections import Counter, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply

REPO = Path(__file__).resolve().parents[6]
sys.path.insert(0, str(REPO))
from model.executable_model import cell, federation


def around(cycle):
    return list(zip(cycle, (*cycle[1:], cycle[0])))


class Patch:
    def __init__(self):
        self.cell = cell()
        self.edges = self.cell.seams
        self.bit = {e: 1 << i for i, e in enumerate(self.edges)}
        self.faces = self.cell.faces()
        self.rings = [self.orientation(r) for r in self.faces]
        self.dual = nx.Graph()
        for i, f in enumerate(self.faces):
            for j, g in enumerate(self.faces[:i]):
                common = tuple(sorted(set(f).intersection(g)))
                if len(common) == 2:
                    self.dual.add_edge(i, j, bit=self.bit[common])
        disk = [j for j, face in enumerate(self.faces) if set(face) & set(self.faces[0])]
        counts = Counter(tuple(sorted(e)) for j in disk for e in around(self.faces[j]))
        boundary, = nx.cycle_basis(nx.Graph(sorted(e for e, count in counts.items() if count == 1)))
        self.disk = frozenset(disk)
        self.boundary = min(tuple(seq[i:] + seq[:i]) for seq in (boundary, boundary[::-1])
                            for i in range(len(seq)))
        distance = nx.single_source_shortest_path_length(self.dual, 0)
        self.opposite, = [j for j, d in distance.items() if d == max(distance.values())]
        self.seed = sum(self.bit[a, b] for a, b in self.edges
                        if sha256(f"{a},{b}".encode()).digest()[0] & 1)
        for loop in (self.faces[0], self.faces[self.opposite], self.boundary):
            mask, cw = self.orientation(loop)
            self.seed = (self.seed & ~mask) | cw
        assert len(self.disk) == 10 and len(self.boundary) == 6
        assert all(set(self.boundary).isdisjoint(self.faces[i]) for i in (0, self.opposite))

    def orientation(self, cycle):
        mask = sum(self.bit[tuple(sorted(e))] for e in around(cycle))
        cw = sum(self.bit[a, b] for a, b in around(cycle) if a < b)
        return mask, cw

    def ring_moves(self, word):
        for i, (mask, cw) in enumerate(self.rings):
            if word & mask in (cw, mask ^ cw):
                yield i, word ^ mask

    def path(self, a, b):
        path = nx.shortest_path(self.dual, a, b)
        return sum(self.dual[u][v]['bit'] for u, v in zip(path, path[1:]))

    def move(self, pos, word, edge, sign):
        a, b = edge
        if pos in edge and ((pos == a) == bool(word & self.bit[edge])) == (sign == 1):
            return (b if pos == a else a), word ^ self.bit[edge]
        return pos, word

    def charge_moves(self, pos, word, sign):
        for other in self.cell.neighbours[pos]:
            new = self.move(pos, word, tuple(sorted((pos, other))), sign)
            if new[0] != pos:
                yield new

    def compensator(self, anchor, sign):
        board = federation(0)
        queue = deque([(anchor,)])
        exits = {p: q for p, q in board.external if p < 12}
        while queue:
            path = queue.popleft()
            if path[-1] in exits:
                return (*path, exits[path[-1]])
            for other, _ in self.charge_moves(path[-1], self.seed, sign):
                if other not in path:
                    queue.append((*path, other))
        raise ValueError('No external injection path for this charge')

    def field_component(self):
        seen = {self.seed}
        queue = deque([self.seed])
        while queue:
            for _, word in self.ring_moves(queue.popleft()):
                if word not in seen:
                    seen.add(word)
                    queue.append(word)
        return sorted(seen)


def phase(word, mask):
    return -1 if (word & mask).bit_count() % 2 else 1


def tour(patch, words, sign):
    at = patch.boundary[0]
    steps = [tuple(sorted(e)) for e in around(patch.boundary)]
    perimeter, oriented = patch.orientation(patch.boundary)
    records = Counter()
    winding_words = []
    word_set = set(words)
    for word in words:
        pos, final = at, word
        for e in steps + steps[-2::-1]:
            pos, final = patch.move(pos, final, e, sign)
        if pos != at:
            records['away'] += 1
        elif final == word:
            records['returned_idle'] += 1
        else:
            assert final == word ^ perimeter and word & perimeter in (oriented, perimeter ^ oriented)
            assert final in word_set
            records['returned_after_winding'] += 1
            winding_words.append((word, final))
        pos, final = at, word
        for e in steps + steps[::-1]:
            pos, final = patch.move(pos, final, e, sign)
        assert (pos, final) == (at, word)
    return records, winding_words


def source_scan(patch, words):
    edges = [[] for _ in patch.faces]
    ids = {w: i for i, w in enumerate(words)}
    for i, word in enumerate(words):
        for f, other in patch.ring_moves(word):
            if word < other:
                edges[f].append((i, ids[other]))
    edge_arrays = [np.asarray(pairs, dtype=int).reshape(-1, 2).T for pairs in edges]
    tours = {s: tour(patch, words, s) for s in (-1, 1)}
    diagonal_energy = 0.5 * sum(2*len(i) for i, _ in edge_arrays) / len(words)
    trials = []
    for a, b in combinations(range(len(patch.faces)), 2):
        path = patch.path(a, b)
        phases = np.array([phase(w, path) for w in words])
        energy = [0.5 * float(np.sum((phases[i] - phases[j]) ** 2)) / len(words)
                  for i, j in edge_arrays]
        assert all(e == 0 for f, e in enumerate(energy) if f not in (a, b))
        # Matched diagonal ring law counts eligible words and ignores relative phases.
        diagonal_by_face = [0.5 * float(np.sum(phases[i]**2 + phases[j]**2)) / len(words)
                            for i, j in edge_arrays]
        np.testing.assert_allclose(sum(diagonal_by_face), diagonal_energy, atol=1e-13)
        np.testing.assert_array_equal(diagonal_by_face, [len(i) / len(words) for i, _ in edge_arrays])
        assert np.array_equal(phases*phases, np.ones(len(words), dtype=int))
        # Different routes between the endpoints differ by a closed dual cycle.
        # Its parity is constant on this fixed-divergence pure-field component.
        first_edge = next(iter(zip(nx.shortest_path(patch.dual, a, b),
                                   nx.shortest_path(patch.dual, a, b)[1:])))
        cut = patch.dual.copy()
        cut.remove_edge(*first_edge)
        alternate = nx.shortest_path(cut, a, b)
        alt_mask = sum(patch.dual[u][v]['bit'] for u,v in zip(alternate, alternate[1:]))
        assert len({phase(w, path ^ alt_mask) for w in words}) == 1
        inside = int(a in patch.disk) + int(b in patch.disk)
        by_sign = {}
        for sign, (counts, paths) in tours.items():
            signed_count = sum(phase(w, path) * phase(final, path) for w, final in paths)
            assert len(paths) > 0 and signed_count == (-1)**inside * len(paths)
            overlap = (signed_count + counts['returned_idle']) / len(words)
            by_sign[str(sign)] = {'conditional_winding_amplitude': signed_count / len(paths),
                                  'inclusive_reference_overlap': overlap,
                                  'inclusive_bright_port': (1 + overlap) / 2}
        trials.append({'ends': [a, b], 'enclosed': inside, 'field_energy': energy,
                       'both_ends_active': min(energy[a], energy[b]) > 0, 'tour': by_sign})
    return {'trials': trials, 'branch_counts': {str(s): dict(c) for s, (c, _) in tours.items()},
            'vacuum_bright_port': {str(s): (1 + (c['returned_idle'] + len(p)) / len(words)) / 2
                                   for s, (c, p) in tours.items()},
            'diagonal_ring_control_excess_at_any_face': 0,
            'diagonal_ring_control_total_energy_both_preparations': diagonal_energy,
            'alternate_dual_route_equivalent_up_to_global_phase_for_all_trials': True,
            'compensation_paths_by_sign': {str(s): patch.compensator(patch.boundary[0], s) for s in tours},
            'zero_winding_bright_port_all_inputs': 1}


def displaced_endpoint_filters(patch, basis, positions):
    """Select a receiving field with its near endpoint moved, orthogonal to bare and old fields."""
    filters = []
    source, remote = 0, patch.opposite
    excluded = set(patch.faces[source]) | set(patch.faces[remote])
    baseline_mask = patch.path(source, remote)
    for p in sorted(set(positions) - excluded):
        idx = np.flatnonzero(positions == p)
        words = [basis[i][1] for i in idx]
        bare = np.ones(len(idx)) / np.sqrt(len(idx))
        old = bare * [phase(w, baseline_mask) for w in words]
        blockers = np.stack([bare, old], axis=1)
        u, sv, _ = np.linalg.svd(blockers, full_matrices=False)
        u = u[:, sv > 1e-12]
        for face, ports in enumerate(patch.faces):
            if p not in ports or set(ports) & set(patch.faces[source]) or face == remote:
                continue
            target_mask = patch.path(face, remote)
            target = bare * [phase(w, target_mask) for w in words]
            unique = target - u @ (u.T @ target)
            size = np.linalg.norm(unique)
            if size > 1e-12:
                filters.append((p, face, idx, unique/size, float(size**2)))
    return filters


def release(patch, neutral, sign):
    anchor = patch.faces[0][0]
    start = [(anchor, w) for w in neutral]
    basis = list(start)
    ids = {state: i for i, state in enumerate(basis)}
    ring_pairs = [[] for _ in patch.faces]
    left, right, weights = [], [], []

    def link(i, state, weight, face=None):
        j = ids.get(state)
        if j is None:
            j = len(basis)
            ids[state] = j
            basis.append(state)
        if i < j:
            left.append(i)
            right.append(j)
            weights.append(weight)
            if face is not None:
                ring_pairs[face].append((i, j))

    for i, (pos, word) in enumerate(basis):
        for f, other in patch.ring_moves(word):
            link(i, (pos, other), 0.5, f)
        if sign:
            for other in patch.charge_moves(pos, word, sign):
                link(i, other, 4.)
    left, right, weights = np.array(left), np.array(right), np.array(weights)
    n = len(basis)
    graph = coo_matrix((weights, (left, right)), shape=(n, n)).tocsr()
    graph = graph + graph.T
    pos = np.array([p for p, w in basis])
    if sign:
        for x, w in basis:
            delta = [0]*12
            for (a, b), bit in patch.bit.items():
                change = int(bool(w & bit)) - int(bool(patch.seed & bit))
                delta[a] += change
                delta[b] -= change
            assert all(d == sign*(int(p == x)-int(p == anchor)) for p, d in enumerate(delta))
    degree = federation(0).degree[:12]
    local_potential = 8. * (6 - degree[pos]) if sign else np.zeros(n)
    h = diags(np.asarray(graph.sum(axis=1)).ravel() + local_potential) - graph
    path = patch.path(0, patch.opposite)
    inputs = np.zeros((n, 2), dtype=complex)
    inputs[:len(start), 0] = 1 / np.sqrt(len(start))
    inputs[:len(start), 1] = [phase(w, path) / np.sqrt(len(start)) for _, w in start]
    pairs = [np.array(e, dtype=int).reshape(-1, 2).T for e in ring_pairs]
    filters = displaced_endpoint_filters(patch, basis, pos) if sign else []
    joint_records = []
    filtered = []
    for k, states in enumerate(expm_multiply(-1j*h, inputs, start=0., stop=2., num=17)):
        for p, face, idx, detector, orthogonal_weight in filters:
            probs = abs(detector @ states[idx])**2
            filtered.append({'time': k/8, 'port': int(p), 'new_face': face,
                             'distinct_fraction_of_target': orthogonal_weight,
                             'uniform_input_probability': float(probs[0]),
                             'string_input_probability': float(probs[1])})
        reads = []
        for state in states.T:
            location = np.bincount(pos, weights=abs(state)**2, minlength=12)
            np.testing.assert_allclose(location.sum(), 1, atol=1e-11)
            joint = np.array([np.bincount(pos[i], weights=0.5*abs(state[i]-state[j])**2, minlength=12)
                              for i, j in pairs])
            by_face = joint.sum(axis=1)
            total = by_face.sum()
            endpoints = by_face[0] + by_face[patch.opposite]
            incident = sum(joint[f, patch.faces[f]].sum() for f in range(20))
            product = sum(by_face[f] * location[list(patch.faces[f])].sum() for f in range(20))
            reads.append({'charge_port_probabilities': location.tolist() if sign else None,
                          'energy_by_face': by_face.tolist(),
                          'endpoint_energy_fraction': float(endpoints/total) if total > 1e-15 else None,
                          'energy_incident_to_charge_fraction': float(incident/total) if sign and total>1e-15 else None,
                          'independent_marginal_fraction': float(product/total) if sign and total>1e-15 else None,
                          'initial_state_survival': float(abs(np.vdot(inputs[:,len(reads)], state))**2)})
        delta = (np.abs(np.array(reads[0]['charge_port_probabilities'])
                        - reads[1]['charge_port_probabilities']).sum()/2) if sign else None
        joint_records.append({'time': k/8, 'uniform': reads[0], 'string': reads[1],
                              'charge_total_variation': delta})
    return {'dimension': n, 'sign': sign, 'anchor': anchor,
            'compensation_path': patch.compensator(anchor, sign) if sign else None,
            'source_overlap_uniform_string': float(np.vdot(inputs[:,0], inputs[:,1]).real),
            'records': joint_records, 'displaced_endpoint_filters': filtered}


def run():
    patch = Patch()
    words = patch.field_component()
    return {'construction': 'positive written Laplacians, J=0.5, hop=4, one-load reference',
            'seed': patch.seed, 'face_ports': patch.faces, 'opposite_face': patch.opposite,
            'boundary': patch.boundary, 'inside_faces': sorted(patch.disk),
            'neutral_component_size': len(words), 'scan': source_scan(patch, words),
            'release': [release(patch, words, s) for s in (0, 1, -1)]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    target = parser.parse_args().output
    result = run()
    with target.open('x', encoding='utf-8') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
