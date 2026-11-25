"""Compare native spatial exchange with an alternating core's frame reflection."""

import argparse
import gzip
import hashlib
from itertools import combinations, product
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.openwave.row_13_baryon.arm_symmetry import measure as source
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import carrier
from model.executable_model.symmetry import rotations
from ..projected_turns.measure import compare

HERE = Path(__file__).resolve().parent
SHOTS, SEED = 2**26, 807331
EDGES6 = list(combinations(range(6), 2))


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def draw(probabilities, rng):
    p = np.asarray(probabilities, float)
    check(np.all(np.isfinite(p)) and min(p) >= 0, 'Invalid probability')
    check(abs(p.sum() - 1) < 5e-11, 'Incomplete record')
    order = np.argsort(p, kind='stable')
    positive = p[order].copy()
    positive[-1] = 1 - positive[:-1].sum()
    check(positive[-1] >= 0, 'Negative explicit complement')
    counts = np.empty(len(p), dtype=np.int64)
    counts[order] = rng.multinomial(SHOTS, positive)
    return dict(probabilities=p.tolist(), counts=counts.tolist(),
                sampling_complement_l1=float(abs(positive[-1] - p[order[-1]])))


def pair_lookup():
    """Exactly the six-vertex graphs with two connected components of size three."""
    lookup = np.zeros(2**15, bool)
    for rest in combinations(range(1, 6), 2):
        first = (0,) + rest
        second = tuple(x for x in range(6) if x not in first)
        codes = []
        for group in (first, second):
            indices = [EDGES6.index(tuple(sorted(e))) for e in combinations(group, 2)]
            codes.append([sum(1 << indices[j] for j in range(3) if bits >> j & 1)
                          for bits in (3, 5, 6, 7)])
        for a, b in product(*codes):
            lookup[a | b] = True
    check(lookup.sum() == 160, 'Unexpected component census')
    return lookup


def connected_pairs(positions, adjacency, lookup):
    code = np.zeros(len(positions), dtype=np.uint16)
    for i, (a, b) in enumerate(EDGES6):
        code |= adjacency[positions[:, a], positions[:, b]].astype(np.uint16) << i
    return lookup[code]


def endpoint_word(graph, a, b, variant=0):
    paths = nx.shortest_simple_paths(graph, a, b)
    path = next(paths)
    if variant:
        path = next(paths)
    edges = list(zip(path, path[1:]))
    return edges + edges[-2::-1]


def permutation_word(graph, permutation, variant=0):
    seen, word = set(), []
    for root in range(len(permutation)):
        if root in seen:
            continue
        cycle, x = [], root
        while x not in seen:
            seen.add(x)
            cycle.append(x)
            x = int(permutation[x])
        for x in cycle[1:]:
            word += endpoint_word(graph, root, x, variant)
    return word


def geometry(sec, arms):
    board = carrier.federation(0)
    graph = nx.Graph()
    graph.add_nodes_from(range(board.n))
    graph.add_edges_from(board.seams)
    proper = set(rotations())
    gates = set(np.flatnonzero(board.degree[:12] == 6))
    face = tuple(arms[:, 2])
    group = [p for p in carrier.automorphisms(carrier.cell().seams)
             if {p[i] for i in face} == set(face) and {p[i] for i in gates} == gates]
    group_data = []
    reflection = None
    for p in group:
        perm = [face.index(p[x]) for x in face]
        parity = (-1) ** sum(a > b for a, b in combinations(perm, 2))
        determinant = 1 if tuple(p) in proper else -1
        check(parity == determinant, 'Arm parity differs from spatial orientation')
        group_data.append(dict(ports=list(p), spatial_determinant=determinant, arm_sign=parity))
        if determinant == -1 and reflection is None:
            reflection = p
    check(len(group) == 6 and sum(tuple(p) in proper for p in group) == 3, 'Wrong native stabilizer')
    choice = None
    for owner in range(1, board.cells):
        target_gates = set(np.flatnonzero(board.degree[12 * owner:12 * owner + 12] == 6))
        frame = next((r for r in sorted(proper) if {r[i] for i in gates} == target_gates), None)
        if frame is None:
            continue
        distances = nx.multi_source_dijkstra_path_length(graph, range(12))
        distance = min(distances[12 * owner + i] for i in range(12))
        if distance >= 4:
            choice = owner, frame, distance
            break
    check(choice is not None, 'No separated proper-frame source')
    owner, frame, distance = choice
    # The actual second guarded source has the same gateway diagonal in this frame.
    target_count = np.array([sum(board.degree[12 * owner + frame[x]] == 6 for x in w) for w in sec.states])
    native_count = np.array([sum(board.degree[x] == 6 for x in w) for w in sec.states])
    check(np.array_equal(target_count, native_count), 'Second native source differs')
    exchange = np.arange(board.n)
    for x in range(12):
        exchange[x] = 12 * owner + frame[x]
        exchange[12 * owner + frame[x]] = x
    reflected = np.arange(board.n)
    reflected[:12] = reflection
    # Pick a genuine +72-degree proper template turn, not a 144-degree generator.
    pts, seams = carrier.icosahedron_from_coordinates()
    iso = carrier.isomorphism(carrier.cell().seams, seams)
    xyz = pts[iso]
    turn = None
    for r in sorted(proper):
        power = np.arange(12)
        for _ in range(5):
            power = np.array(r)[power]
        matrix = np.linalg.lstsq(xyz, xyz[list(r)], rcond=None)[0]
        if np.array_equal(power, np.arange(12)) and 1.5 < np.trace(matrix) < 1.7:
            turn = np.arange(board.n)
            turn[:12] = r
            break
    check(turn is not None, 'No 72-degree template turn')
    return board, graph, owner, frame, exchange, reflected, turn, dict(
        second_owner=owner, second_frame=list(frame), source_support_min_distance=distance,
        arm_stabilizer=group_data, proper_turn_ports=turn[:12].tolist(),
        improper_reflection_ports=list(reflection))


def run():
    rng = np.random.default_rng(SEED)
    sec, arms, _, _, _, ideal, finite, summary = source.source()
    board, graph, owner, frame, exchange, reflection, turn, geometry_data = geometry(sec, arms)
    words = np.asarray(sec.states)
    inverse = np.argsort(frame)
    positions = np.concatenate([np.repeat(words, 220, axis=0),
                               np.tile(12 * owner + np.array(frame)[words], (220, 1))], axis=1)
    amplitude = np.outer(finite, finite).ravel()
    finite_weights = abs(amplitude)**2
    ideal_weights = np.outer(ideal**2, ideal**2).ravel()
    lookup = pair_lookup()
    adjacency = board.adjacency().toarray().astype(bool)
    index = np.full((12, 12, 12), -1, dtype=int)
    index[tuple(words.T)] = np.arange(220)
    masses = np.asarray(summary['filter_acceptances'])
    outcomes = np.r_[masses[:-1] - masses[1:], masses[-1]]
    fail = np.outer(outcomes, outcomes).ravel()[:-1]
    acceptance = float(masses[-1]**2)
    initial_pair = float(ideal_weights @ connected_pairs(positions, adjacency, lookup))
    exchange0 = permutation_word(graph, exchange, 0)
    exchange1 = permutation_word(graph, exchange, 1)
    reflection_word = permutation_word(graph, reflection)
    turn_word = permutation_word(graph, turn)
    protocols = [('exchange_shortest', exchange0, exchange),
                 ('exchange_alternative_paths', exchange1, exchange),
                 ('exchange_then_improper_frame', exchange0 + reflection_word, reflection[exchange]),
                 ('five_proper_72_degree_turns', turn_word * 5, np.arange(board.n))]
    calibrations = []
    for contrast in (8., 4., 7.5):
        gate = expm(-1j * np.pi * law.seam_term(contrast))
        expected = np.zeros((9, 9), complex)
        for j, (a, b) in enumerate(law.SRC):
            expected[law.SRC.index((b, a)), j] = np.exp(-1j * contrast * np.pi * (a - b)**2)
        check(np.max(abs(gate - expected)) < 2e-13, 'Native contact phase calibration fails')
        calibrations.append(dict(contrast=contrast, ternary_gate_error=float(np.max(abs(gate - expected))),
                                 unequal_binary_swap_phase=[float(expected[law.SRC.index((0, 1)), law.SRC.index((1, 0))].real),
                                                            float(expected[law.SRC.index((0, 1)), law.SRC.index((1, 0))].imag)]))
    reports = []
    for label, pulse_word, target in protocols:
        mapping = np.arange(board.n)
        unequal_count = np.zeros(220**2, dtype=int)
        stages = []
        for stage in range(len(pulse_word) + 1):
            physical = mapping[positions]
            probability = float(ideal_weights @ connected_pairs(physical, adjacency, lookup))
            finite_probability = float(finite_weights @ connected_pairs(physical, adjacency, lookup))
            record = draw([acceptance * finite_probability, acceptance * (1 - finite_probability), 1 - acceptance], rng)
            stages.append(dict(stage=stage, ideal_two_separate_connected_triples=probability,
                               finite_two_separate_connected_triples=finite_probability,
                               public_component_read=record))
            if stage == len(pulse_word):
                break
            a, b = pulse_word[stage]
            unequal_count += np.any(physical == a, axis=1) ^ np.any(physical == b, axis=1)
            mapping = np.where(mapping == a, b, np.where(mapping == b, a, mapping))
        check(np.array_equal(mapping, target), 'Physical contact word implements another spatial map')
        final = mapping[positions]
        left = np.sort(np.where(final < 12, final, board.n), axis=1)[:, :3]
        right = np.sort(np.where(final // 12 == owner, final % 12, 12), axis=1)[:, :3]
        check(np.max(left) < 12 and np.max(right) < 12, 'Final cores leave their two cells')
        right = np.sort(inverse[right], axis=1)
        destination = index[tuple(left.T)] * 220 + index[tuple(right.T)]
        check(len(np.unique(destination)) == 220**2, 'Final occupation words collide')
        phases = []
        for contrast in (8., 4., 7.5):
            phase = np.ones(220**2, complex) if contrast in (8., 4.) else (1j)**(unequal_count % 4)
            output = np.zeros_like(amplitude, complex)
            output[destination] = phase * amplitude
            plus = float(np.linalg.norm((output + amplitude) / 2)**2)
            minus = float(np.linalg.norm((output - amplitude) / 2)**2)
            bank = draw(np.r_[fail, acceptance * plus, acceptance * minus], rng)
            phases.append(dict(contrast=contrast, conditional_plus=plus, conditional_minus=minus,
                               conditional_overlap=[float(np.vdot(amplitude, output).real), float(np.vdot(amplitude, output).imag)],
                               complete_source_histories_and_phase=bank))
        reports.append(dict(protocol=label, physical_pi_contacts=[list(e) for e in pulse_word],
                            stages=stages, minimum_public_two_core_probability=min(x['ideal_two_separate_connected_triples'] for x in stages),
                            phases=phases))
    return dict(shots=SHOTS, seed=SEED, source=summary, geometry=geometry_data,
                pair_source_acceptance=acceptance,
                source_refusal_outcome_order='coreA first refusal1/2/3/4 or success, then coreB same; omit final success-success',
                phase_outcome_order='24 source-refusal histories, success-plus, success-minus',
                component_outcome_order='source-success and two connected triples, source-success other, source-refusal',
                initial_two_core_probability=initial_pair,
                finite_pair_vector_error_upper=2 * summary['finite_to_ideal_vector_upper'],
                stage_record_scope='Counts use the actual finite source at each endpoint; ideal-line convex interpolation bounds each full contact with the declared finite-source error',
                contact_calibrations=calibrations, protocols=reports,
                source_sha256={str(p.relative_to(HERE.parents[4])): hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__), Path(source.__file__), Path(law.__file__))},
                scope='Calibrated guided material transport and ordinary finite spatial representation; no intact free-core exchange or spinor assignment')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--receipt', type=Path, default=HERE / 'results.json.gz')
    args = parser.parse_args()
    data = run()
    serialized = json.dumps(data, indent=1, allow_nan=False).encode()
    if args.output:
        encoded = gzip.compress(serialized, mtime=0) if args.output.suffix == '.gz' else serialized
        with args.output.open('xb') as stream:
            stream.write(encoded)
    else:
        data_old = json.loads(gzip.decompress(args.receipt.read_bytes()))
        compare(data, data_old, 'mobile_core_exchange')
        print('Verified selected-core source, physical contacts, component histories and exchange phases.')


if __name__ == '__main__':
    main()
