"""A native extended source pulse probes the same freely arrived neutral pair."""

import argparse
from collections import deque
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from .. import measure as flight
from ..cell_field_probe import measure as previous
from . import circuit, physical

HERE = Path(__file__).resolve().parent
SHOTS, SEED, TIME = 2**26, 191619, 1_000_000


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def routes(board, region):
    choices = []
    for b in range(board.n):
        if b // 12 in region or board.degree[b] != 5:
            continue
        for a, j in product(board.neighbours[b], repeat=2):
            if a == j or a // 12 != b // 12 or j // 12 != b // 12 or board.degree[j] != 5:
                continue
            choices.append((b, a, j))
    b, a, j = min(choices)
    parents = {j: None}
    queue = deque([j])
    while queue:
        current = queue.popleft()
        for other in board.neighbours[current]:
            if other not in (a, b) and other not in parents:
                parents[other] = current
                queue.append(other)
    target_ports = [p for p in range(board.n) if p // 12 in region]
    contact_pulses = 0
    longest = 0
    paths = []
    for target in target_ports:
        check(target in parents, 'source cannot reach a target without disturbing meter')
        path = [target]
        while path[-1] != j:
            path.append(parents[path[-1]])
        length = len(path) - 1
        # Path forward and inverse internal edges implement an endpoint transposition;
        # repeat it after the contact to restore every unknown intermediate load.
        contact_pulses += 2 * (2 * length - 1)
        longest = max(longest, length)
        paths.append(path)
    return dict(source_park=a, source_cross=b, interaction_port=j,
                target_ports=target_ports, endpoint_transposition_paths=paths,
                routing_contact_pulses=contact_pulses,
                diagonal_contacts=4 * len(target_ports), diagonal_references=6 * len(target_ports),
                maximum_path_edges=longest)


def regions(board):
    adjacency = flight.sparse.csr_matrix((np.ones(2 * len(board.external)),
                      (np.r_[[a // 12 for a, b in board.external], [b // 12 for a, b in board.external]],
                       np.r_[[b // 12 for a, b in board.external], [a // 12 for a, b in board.external]])),
                      shape=(board.cells, board.cells))
    distance = shortest_path(adjacency, directed=False, unweighted=True)
    result = []
    for radius in (0, 1, 2):
        for owner in range(board.cells):
            selected = np.flatnonzero(distance[owner] <= radius).tolist()
            result.append(dict(center=owner, radius=radius, owners=selected, apparatus=routes(board, selected)))
    return result


def sample(probabilities, rng):
    raw = np.asarray(probabilities, dtype=float)
    check(np.all(np.isfinite(raw)) and min(raw) >= 0, 'negative complete detector probability')
    check(abs(sum(raw) - 1) < 4e-11, 'complete detector normalization failed')
    order = np.argsort(raw, kind='stable')
    positive = raw[order].copy()
    positive[-1] = 1 - sum(positive[:-1])
    check(positive[-1] >= 0, 'negative sampling complement')
    counts = np.empty(len(raw), np.int64)
    counts[order] = rng.multinomial(SHOTS, positive)
    return dict(probabilities=raw.tolist(), counts=counts.tolist(),
                sampling_l1_complement_bound=float(abs(positive[-1] - raw[order[-1]])))


def panel(weights, x, y, region, admission, meters, rng):
    q = np.isin(x // 12, region['owners']).astype(int) - np.isin(y // 12, region['owners']).astype(int)
    apparatus = region['apparatus']
    a, b = apparatus['source_park'], apparatus['source_cross']
    na = 1 + (x == a).astype(int) - (y == a).astype(int)
    nb = 1 + (x == b).astype(int) - (y == b).astype(int)
    contrast = na - nb
    check(np.all((contrast == 0) == ((na == 1) & (nb == 1))),
          'eligibility zero is not uniquely a unit-filled meter')
    eligibility = np.bincount(contrast + 2, weights=weights, minlength=5)
    admitted = contrast == 0
    law = np.bincount(q[admitted] + 1, weights=weights[admitted], minlength=3)
    all_charge = np.bincount(q + 1, weights=weights, minlength=3)
    records = []
    for multiplier in (1, 2):
        for axis in ('X', 'Y'):
            meter = sum(law[i] * meters[multiplier][qvalue][axis] for i, qvalue in enumerate((-1, 0, 1)))
            # Four nonzero eligibility records, nine final meter words and source-filter miss.
            complete = np.r_[admission * eligibility[[0, 1, 3, 4]], admission * meter, 1 - admission]
            records.append(dict(multiplier=multiplier, axis=axis, **sample(complete, rng)))
    return dict(eligible_probability=float(eligibility[2]), eligibility_by_contrast=eligibility.tolist(),
                admitted_region_charge=law.tolist(), all_region_charge=all_charge.tolist(),
                excluded_region_charged_probability=float(sum(all_charge[[0, 2]]) - sum(law[[0, 2]])),
                accepted_charged_fraction=float(sum(law[[0, 2]]) / eligibility[2]),
                complete_records=records)


def run():
    meters, circuit_audit = circuit.audit()
    board = flight.federation(0)
    x, y, _, kinetic, h, cost = flight.defect.generators(board)
    original = json.loads((flight.HERE / 'results.json').read_text())
    reduction = flight.PairReduction(h, cost)
    energies, basis, residual, _ = reduction.lowest(h, 24)
    check(reduction.inertia(flight.CUT)['upper_lines_below'] == 24, 'original pair band count differs')
    projector_error = float(max(original['cluster_projector_error'], residual / (flight.CUT - energies[-1])))
    profiles = regions(board)
    rng = np.random.default_rng(SEED)
    calibrations = []
    for multiplier in (1, 2):
        for qvalue in (-1, 0, 1):
            for axis in ('X', 'Y'):
                calibrations.append(dict(multiplier=multiplier, region_charge=qvalue, axis=axis,
                                         **sample(meters[multiplier][qvalue][axis], rng)))
    outputs = []
    for launch in original['flights']:
        if launch['owner'] not in (0, 5, 10) or launch['admission'] < 0.001:
            continue
        index = flight.defect.ordered_index(*launch['ports'], board.n)
        overlap = basis[index]
        admission = float(overlap @ overlap)
        state = basis @ (np.exp(-1j * (energies - energies[0]) * TIME) * overlap / np.sqrt(admission))
        rival = expm_multiply(-4j * kinetic, state)
        weights, norm_error = previous.normalized_weights(state)
        changed, changed_norm = previous.normalized_weights(rival)
        matched, matched_info = previous.hard_core_control(weights, x, y, board.n)
        effects = {}
        for label, probability in (('native_arrival', weights), ('matched_marginal_model', matched),
                                    ('exchange_only_four_more_units', changed)):
            effects[label] = [panel(probability, x, y, region, admission, meters, rng) for region in profiles]
        error = float(2 * projector_error / np.sqrt(admission) + (TIME + 4) * residual
                      + max(norm_error, changed_norm) + 1e-9)
        outputs.append(dict(owner=launch['owner'], initial_ports=launch['ports'],
                            initial_filter_admission=admission, filter_refusal=1 - admission,
                            native_event_probability_error=error, matched_model=matched_info,
                            panels=effects))
    check(len(outputs) == 3, 'not every native pair launch was retained')
    return dict(shots=SHOTS, seed=SEED, free_flight=TIME, source_circuit=circuit_audit,
                physical_compensated_standards=physical.audit(board, profiles, meters),
                pointer_decoder=previous.pointer_decoder(), profiles=profiles,
                charged_standards=calibrations, launches=outputs,
                complete_record_order='eligibility contrast -2,-1,+1,+2; meter words0..8; initial filter refusal',
                source_prep='After a native signed-zero host herald, supply fresh neutral20 input with fixed total load',
                scope='Native source-induced coarse-charge Ramsey measurement; source routing/profile/clock supplied, electromagnetism and inertial mass unestablished')


def compare(a, b, path='native_source_probe'):
    if isinstance(b, dict):
        check(isinstance(a, dict) and a.keys() == b.keys(), path)
        for key in b:
            compare(a[key], b[key], path + '.' + key)
    elif isinstance(b, list):
        check(isinstance(a, list) and len(a) == len(b), path)
        for i, (left, right) in enumerate(zip(a, b, strict=True)):
            compare(left, right, f'{path}[{i}]')
    elif isinstance(b, float):
        check(np.isfinite(a) and np.isclose(a, b, rtol=1e-9, atol=2e-10), path)
    else:
        check(type(a) is type(b) and a == b, path)


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output', type=Path)
    cli.add_argument('--receipt', type=Path, default=HERE / 'results.json')
    args = cli.parse_args()
    result = run()
    if args.output:
        payload = json.dumps(result, indent=2, allow_nan=False) + '\n'
        with args.output.open('x') as stream:
            stream.write(payload)
    else:
        compare(result, json.loads(args.receipt.read_text()))
        print('Replayed native-source field records, host eligibility, routing and charged standards.')
