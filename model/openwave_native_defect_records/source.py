"""Conservation and charge tests on the literal native repair process.

Positions are decoded from a declared preparation and complete current records.
No hidden source load is supplied to the observer's decoding function. The
modulo-six history statistic is separated from a protected source charge.
"""
from __future__ import annotations

from collections import defaultdict, deque
from fractions import Fraction as F
from itertools import combinations
from math import comb

from model.openwave_native_process import source as process


def load_positions(mask):
    return tuple(i for i in range(mask.bit_length()) if (mask >> i) & 1)


def record_current_chain(records):
    """Observer function of committed records only, without terminal source."""
    currents = defaultdict(int)
    for record in records:
        if record.kind == 'hold':
            continue
        if record.kind != 'repair' or record.current is None:
            raise ValueError('a positive-probability complete native record required')
        currents[record.seam] += record.current
    return {edge: value for edge, value in currents.items() if value}


def decode_positions(initial_positions, records):
    """Complete known-preparation position decoder, with retained failure."""
    if len(set(initial_positions)) != len(initial_positions) or any(type(i) is not int or i < 0 for i in initial_positions):
        raise ValueError('distinct nonnegative prepared positions required')
    mask = sum(1 << i for i in initial_positions)
    currents = defaultdict(int)
    catalogue = {row.record: row for row in process.raw_local_branches()}
    for record in records:
        key = (record.kind, record.law, record.coin, record.raw, record.label)
        row = catalogue.get(key)
        if row is None or not row.weight or record.current != row.current:
            return {'status': 'invalid_record', 'positions': (), 'currents': dict(currents)}
        try:
            target = process.apply_branch(mask, record.seam, row)
        except (ValueError, TypeError):
            target = None
        if target is None:
            return {'status': 'incompatible_preparation', 'positions': (), 'currents': dict(currents)}
        currents[record.seam] += record.current
        mask = target
    return {'status': 'identified_load_positions',
            'positions': load_positions(mask),
            'currents': {edge: value for edge, value in currents.items() if value}}


def cut_flux(currents, subset):
    """Signed outward load current, in integer native load units."""
    return sum(value*(int(a in subset)-int(b in subset))
               for (a, b), value in currents.items())


def all_cut_certificate(initial, final, currents):
    """Verify the vertex identity, which implies every global cut identity."""
    active = set(initial) | set(final) | {v for edge in currents for v in edge}
    for v in active:
        assert int(v in initial)-int(v in final) == cut_flux(currents, {v})
    # Enumerate every cut on the actual active support, including empty/full.
    active = tuple(sorted(active))
    for mask in range(1 << len(active)):
        subset = {v for i, v in enumerate(active) if (mask >> i) & 1}
        assert sum(v in subset for v in initial)-sum(v in subset for v in final) == cut_flux(currents, subset)
    return len(active), 1 << len(active)


def _distances(g, origin):
    distances, queue = {origin: 0}, deque([origin])
    while queue:
        x = queue.popleft()
        for y in g.neighbours[x]:
            if y not in distances:
                distances[y] = distances[x]+1
                queue.append(y)
    return distances


def endpoint_motion(level):
    g = process.graph(level)
    a, b = g.inter_seams[0]
    distance_a, distance_b = _distances(g, a), _distances(g, b)
    spectator = min(v for v in range(g.ports) if min(distance_a[v], distance_b[v]) >= 3)
    destination = min(v for v in g.neighbours[b] if v not in (a, spectator))
    seams = ((a, b), tuple(sorted((b, destination))), (a, b))
    initial_positions = (a, spectator)
    initial_mask = (1 << a) | (1 << spectator)
    initial = {(initial_mask, 0, initial_mask, 0): F(1)}
    program, history = process.execute_event_tape(level, 2, seams, initial)
    histogram, modulus = defaultdict(F), [F(0)]*6
    failures, decoded_mass, cuts_checked = F(0), F(0), 0
    for branch in history.branches:
        probability = process.native.trace(branch.density)
        if not probability:
            continue
        decoded = decode_positions(initial_positions, branch.archive+branch.memory)
        if decoded['status'] != 'identified_load_positions':
            failures += probability
            continue
        positions = decoded['positions']
        # The hidden state appears only in the independent verifier, never
        # in decode_positions or any selector/readout input.
        assert set(branch.density) == {(sum(1 << v for v in positions), 0,
                                      sum(1 << v for v in positions), 0)}
        assert spectator in positions and len(positions) == 2
        _, count = all_cut_certificate(initial_positions, positions, decoded['currents'])
        cuts_checked += count
        histogram[positions] += probability
        modulus[cut_flux(decoded['currents'], {a}) % 6] += probability
        decoded_mass += probability
    assert decoded_mass == 1 and failures == 0 and sum(modulus) == 1
    return {'level': level, 'ports': g.ports, 'global_fixed_total': 2,
            'initial_separated_positions': initial_positions,
            'minimum_initial_endpoint_separation': distance_a[spectator],
            'event_seams': seams, 'program': program,
            'complete_raw_history_blocks': len(history.branches),
            'all_raw_histories_sha256': process.compact_history(history)['all_raw_histories_sha256'],
            'exactly_zero_raw_blocks': sum(not b.density for b in history.branches),
            'complete_decoded_position_distribution': [[list(k), str(v)] for k, v in sorted(histogram.items())],
            'all_six_cut_current_mod6_probabilities': list(map(str, modulus)),
            'identification_failure_probability_exact': str(failures),
            'all_other_position_outcomes_probability_exact': '0',
            'complete_reference_support': 'declared known two-load preparation; no unknown-input position inference',
            'native_current_unit': 'one integer occupation; no empirical electric-charge calibration',
            'active_cut_checks': cuts_checked,
            'all_global_cuts_follow_from_vertex_identity': True}


def local_source_algebra():
    """Construct all rank-one projectors and a native transport tree on H66."""
    masks = tuple(sum(1 << i for i in pair) for pair in combinations(range(12), 2))
    index = {mask: i for i, mask in enumerate(masks)}
    seams = process.native.local.SEAMS
    signatures, edges = [], {mask: {} for mask in masks}
    for mask in masks:
        signature = []
        for seam in seams:
            labels = []
            for row in process.raw_local_branches():
                target = process.apply_branch(mask, seam, row)
                if target is None or not row.weight:
                    continue
                if row.record[1:3] == ('E', 0):
                    labels.append(row.record[-1])
                if target != mask:
                    edges[mask].setdefault(target, [list(seam), list(row.record), str(row.weight)])
            assert len(labels) == 1
            signature.append(labels[0])
        signatures.append(tuple(signature))
    assert len(set(signatures)) == len(masks) == 66
    tree, reached, queue = [], {masks[0]}, deque([masks[0]])
    while queue:
        a = queue.popleft()
        for b, record in sorted(edges[a].items()):
            if b not in reached:
                reached.add(b); queue.append(b)
                tree.append([index[a], index[b], record])
    assert len(reached) == 66 and len(tree) == 65
    return {'source': 'complete local H66 binary total-two carrier',
        'dimension': 66, 'actual_E_effect_signatures': [list(row) for row in signatures],
        'configuration_masks': list(masks), 'native_transport_spanning_tree': tree,
        'all_rank_one_projectors_from_E_effect_products': True,
        'generated_star_algebra': 'M66(C)', 'commutant': 'scalar multiples of I66',
        'exact_conserved_sharp_source_observables': 'constants on this fixed-number sector',
        'extra_coherent_memory_included': False}


def general_sector_certificate(level):
    g = process.graph(level)
    distances = _distances(g, 0)
    assert len(distances) == g.ports
    # All actual one-load E-effect signatures are independently distinct.
    signatures = tuple(tuple(int(x == a)-int(x == b) for a, b in g.seams)
                       for x in range(g.ports))
    assert len(set(signatures)) == g.ports
    return {'level': level, 'ports': g.ports, 'connected_ports': len(distances),
        'all_one_load_E_effect_signatures_distinct': True,
        'one_load_dimension': g.ports, 'two_load_dimension': comb(g.ports, 2),
        'sector_theorem': 'For each fixed 0<N<V, connected native swaps and all E differences generate M_binomial(V,N); its source-only commutant is scalar.',
        'two_load_all_configurations_enumerated': False,
        'proof_shortcut': 'edge differences determine a configuration at fixed total; edge transpositions on a connected graph act transitively on N-subsets'}


def winding_control():
    """An actual closed native load history creates a nonzero cycle record."""
    seams = ((0, 1), (1, 5), (0, 5))
    initial = {(1, 0, 1, 0): F(1)}
    program, history = process.execute_event_tape(0, 1, seams, initial)
    target_chain = {(0, 1): 1, (1, 5): 1, (0, 5): -1}
    loop_mass, stationary_mass = F(0), F(0)
    classes, returned_classes = [F(0)]*6, [F(0)]*6
    witnesses = {}
    for branch in history.branches:
        weight = process.native.trace(branch.density)
        if not weight:
            continue
        currents = record_current_chain(branch.archive+branch.memory)
        decoded = decode_positions((0,), branch.archive+branch.memory)
        assert decoded['status'] == 'identified_load_positions'
        all_cut_certificate((0,), decoded['positions'], currents)
        k = currents.get((0, 1), 0) % 6
        classes[k] += weight
        if decoded['positions'] == (0,):
            returned_classes[k] += weight
        if currents == target_chain:
            assert decoded['positions'] == (0,)
            loop_mass += weight
            witnesses.setdefault('cycle', [process.asdict(r) for r in branch.memory])
        if not currents:
            assert decoded['positions'] == (0,)
            stationary_mass += weight
            witnesses.setdefault('hold_position', [process.asdict(r) for r in branch.memory])
    assert loop_mass == F(1, 8) and stationary_mass == F(1, 4)
    assert sum(classes) == 1 and classes[1] > 0
    # An explicit extra repair is a perturbation of the history statistic.
    single = process.native.execute(process.admitted_program(0, 1,
                [process.fixed_command((0, 1))]), initial)
    changed_class = sum(process.native.trace(b.density) for b in single.branches
                       if b.record.current == 1)
    assert changed_class == F(1, 2)
    return {'program': program, 'complete_raw_history_blocks': len(history.branches),
        'all_raw_histories_sha256': process.compact_history(history)['all_raw_histories_sha256'],
        'initial_and_returned_position': 0,
        'nonzero_cycle_chain': [[list(edge), value] for edge, value in sorted(target_chain.items())],
        'cycle_probability_exact': str(loop_mass), 'zero_chain_probability_exact': str(stationary_mass),
        'all_six_memory_class_probabilities': list(map(str, classes)),
        'returned_source_joint_six_class_probabilities': list(map(str, returned_classes)),
        'raw_positive_history_witnesses': witnesses,
        'one_repair_changes_prepared_zero_history_class_probability_exact': str(changed_class),
        'protected_Z6_source_charge_identified': False,
        'centre_flux': 'no native seam-link register or Wilson-loop observable supplied',
        'spatial_winding': 'the recorded triangle is a graph cycle and a boundary face of the full cell complex; no nontrivial sphere winding is inferred',
        'electric_charge': 'integer load current; no Standard Model charge normalization or species dictionary'}


def run():
    return {'verdict': 'PASS_NATIVE_CURRENT_RECORDS_AND_SOURCE_CHARGE_CLASSIFICATION',
        'local_source_algebra': local_source_algebra(),
        'connected_federation_sectors': [general_sector_certificate(level) for level in (0, 1)],
        'complete_endpoint_motion': [endpoint_motion(level) for level in (0, 1)],
        'six_outcome_history_control': winding_control(),
        'memory_instrument': {'residue_outcomes': list(range(6)),
            'seventh_failure_outcome': 'invalid_or_incompatible_records',
            'effect': 'sum of orthogonal complete native-record projectors whose declared integer current statistic has that residue',
            'completeness': 'the six residue sets partition every valid history; invalid raw-register histories retain a separate identification failure of probability zero for the declared preparation',
            'backaction': 'classical postprocessing of records committed by the original instrument; no extra hidden-load measurement',
            'work': 'native attempts, prepared blank records and terminal classical computation counted as services; physical energy is not supplied'},
        'scope': 'Exact source-only conservation classification and native load motion. The modulo-six history output is not a protected centre charge, fractional electric charge or physical quark. A coherent protected memory extension changes the carrier and must specify its own preparation and operations.'}
