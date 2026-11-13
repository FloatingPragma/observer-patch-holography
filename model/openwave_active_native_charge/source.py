"""Complete native load repairs with an active six-state memory.

The additional transport voltage is declared architecture. Its protected
characters are distinct from an electromagnetic unit or a particle assignment.
"""
from __future__ import annotations
from collections import deque
from fractions import Fraction as F
from itertools import combinations
from math import comb
import hashlib
import json

from model.openwave_native_process import source as native
from model.exchange_sector import statistics as exchange

PAIRS = exchange.PAIRS
INDEX = exchange.INDEX
SEAMS = exchange.SEAMS
MEMORY = 6


def voltage(before, after, seam, extra=True):
    """Transported-order parity plus a declared spectator-sensitive voltage."""
    if type(extra) is not bool:
        raise ValueError('an explicit extra-voltage switch is required')
    if len(before) != 2 or len(after) != 2 or tuple(sorted(set(before))) != before or tuple(sorted(set(after))) != after:
        raise ValueError('two distinct sorted source positions are required')
    if before == after:
        return 0
    lost, gained = set(before)-set(after), set(after)-set(before)
    if len(lost) != 1 or len(gained) != 1 or tuple(sorted(lost | gained)) != seam:
        raise ValueError('one collision-free motion on the declared seam is required')
    origin, destination = next(iter(lost)), next(iter(gained))
    transported = tuple(destination if x == origin else x for x in before)
    parity = int(transported != after)
    current = 1 if origin == seam[0] else -1
    spectator = next(x for x in before if x != origin)
    return (3*parity + 2*current*int(extra and seam == (0, 5) and spectator == 1)) % 6


def target(pair, seam, row):
    mask = sum(1 << v for v in pair)
    out = native.apply_branch(mask, seam, row)
    if out is None:
        return None
    result = tuple(v for v in set(pair) | set(seam) if (out >> v) & 1)
    return tuple(sorted(result))


def lifted_column(pair, memory, seam, row, extra=True):
    if type(memory) is not int or memory not in range(6):
        raise ValueError('one of six active memory addresses is required')
    after = target(pair, seam, row)
    return None if after is None else (after, (memory + voltage(pair, after, seam, extra)) % 6)


def _digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), sort_keys=True).encode()).hexdigest()


def local_complete_certificate():
    """Exact effects and shift commutation on every column and raw slot."""
    counts = dict(instruments=0, raw_slots=0, active_columns=0,
                  zero_weight_slots=0, trace_matrix_units_implied_per_instrument=396**2)
    adjacency = [set() for _ in range(396)]
    map_hashes = []
    for law in ('C', 'E', 'mixture'):
        rows = native.raw_local_branches(law)
        assert len(rows) == 25
        for seam in SEAMS:
            effects = [F(0)]*396
            recorded = []
            for row in rows:
                counts['raw_slots'] += 1
                counts['zero_weight_slots'] += int(row.weight == 0)
                images = set()
                columns = []
                for i, pair in enumerate(PAIRS):
                    after = target(pair, seam, row)
                    if after is None:
                        continue
                    j = INDEX[after]
                    shift = voltage(pair, after, seam)
                    if row.weight:
                        assert voltage(after, pair, seam) == (-shift) % 6
                    for m in range(6):
                        before = 6*i+m
                        out = 6*j+(m+shift) % 6
                        assert out not in images
                        images.add(out)
                        effects[before] += row.weight
                        assert (m+1+shift) % 6 == ((m+shift) % 6+1) % 6
                        if row.weight:
                            adjacency[before].add(out)
                        columns.append((before, out))
                        counts['active_columns'] += 1
                recorded.append((row.record, str(row.weight), columns))
            assert effects == [F(1)]*396
            counts['instruments'] += 1
            map_hashes.append([law, list(seam), _digest(recorded)])
    reached, queue = {0}, deque([0])
    while queue:
        for end in adjacency[queue.popleft()]:
            if end not in reached:
                reached.add(end)
                queue.append(end)
    assert len(reached) == 396
    return {**counts, 'full_active_carrier_dimension':396,
        'universal_raw_slots_per_event':25, 'all_effect_sums_exact_identity':True,
        'every_branch_is_partial_injection':True, 'every_branch_commutes_with_X6':True,
        'classical_lift_reachable_configurations':len(reached),
        'complete_raw_map_hashes':map_hashes,
        'arbitrary_reference_extension':'Exact column identities and partial-injection orthogonality imply every source/reference matrix-unit identity; no reference dimension restriction.'}


def federation_extension_certificate():
    """Complete source columns for selected actual seams on the 240-port graph."""
    graph = native.graph(0)
    pairs = tuple(combinations(range(graph.ports), 2))
    # Local moving, inter-cell moving, and untouched/padding controls share one rule.
    seams = ((0, 5), graph.inter_seams[0], graph.inter_seams[-1])
    records = []
    for seam in seams:
        assert seam in graph.seams
        slots = native.raw_local_branches()
        effects = [F(0)]*len(pairs)
        checks = 0
        nonzero = 0
        for row in slots:
            images = set()
            for i, pair in enumerate(pairs):
                after = target(pair, seam, row)
                if after is None:
                    continue
                assert after not in images
                images.add(after)
                effects[i] += row.weight
                shift = voltage(pair, after, seam)
                assert shift in range(6)
                if row.weight:
                    assert voltage(after, pair, seam) == (-shift) % 6
                    nonzero += 1
                # Tensor translation supplies all six memory columns exactly.
                checks += 1
        assert effects == [F(1)]*len(pairs)
        records.append({'seam':list(seam),'source_basis_columns':len(pairs),
            'universal_slots':25,'partial_source_columns':checks,'nonzero_weight_source_columns':nonzero,
            'active_memory_columns_implied':6*checks,'effect_sum_exact_identity':True})
    return {'ports':graph.ports,'two_load_source_dimension':comb(graph.ports,2),
        'joint_dimension':6*comb(graph.ports,2),'complete_source_column_checks':records,
        'all_seam_extension_proof':'For any seam, native endpoint partial injections tensor the unchanged spectators. The configuration-dependent cyclic memory translation is bijective at every nonzero source column. Its inverse voltage is antisymmetric. Completeness and X6 commutation hold for every seam and any finite source/reference carrier.',
        'all_global_seams_numerically_enumerated':False}


def loop(start, moves, law='mixture'):
    pair = start
    total = 0
    probability = F(1)
    witness = []
    for seam, after in moves:
        found = [r for r in native.raw_local_branches(law) if r.weight and target(pair,seam,r)==after]
        assert found
        step_probability = sum(r.weight for r in found)
        shift = voltage(pair,after,seam)
        witness.append({'seam':list(seam),'before':list(pair),'after':list(after),
            'voltage':shift,'probability_exact':str(step_probability),
            'raw_records':[list(r.record)for r in found]})
        probability *= step_probability
        total = (total+shift)%6
        pair = after
    assert pair == start
    return {'initial_pair':list(start),'steps':witness,'voltage_mod6':total,
        'probability_exact':str(probability),'other_histories_retained':True}


def charge_boundaries():
    halves = {str(k):[j for j in range(6)if 2*j%6==k]for k in range(6)}
    assert halves['3']==[] and halves['2']==[1,4]
    # Every additive real-valued charge on Z6 obeys 6 q(1)=q(0)=0.
    selected = [F(0),F(1,3),F(2,3),F(1),F(-2,3),F(-1,3)]
    assert selected[3]+selected[3] != selected[(3+3)%6]
    return {'equal_class_split_solutions_mod6':halves,
        'class_three_has_no_equal_modular_half':True,
        'electron_assignment_to_class_three_is_supplied':True,
        'class_two_has_equal_half_control':[1,1],
        'additive_real_charge_from_Z6_alone':'Every homomorphism Z6 to additive R is zero; finite centre class alone cannot supply a nonzero additive electric charge.',
        'optional_charge_representatives_in_supplied_unit':list(map(str,selected)),
        'representative_nonadditivity_control':{'class_three_plus_class_three_mod6':0,
            'sum_assigned_charge_exact':'2','representative_of_sum_exact':'0'},
        'integer_lift_Gauss_constraint_or_electric_unit_derived':False,
        'spatial_partial_winding_or_physical_quark_string_derived':False,
        'physical_charge_criterion_satisfied':False,'physical_quark_criterion_satisfied':False}


def localization_and_memory_error():
    # Characters have the same diagonal source distribution and all native
    # records for diagonal preparation; the phase-sensitive service differs.
    probabilities = [F(1,6)]*6
    return {'source_reference':'I66/66',
        'sharp_joint_states':'I66/66 tensor |k><k|, where X6|k>=exp(2*pi*i*k/6)|k>',
        'different_character_joint_trace_distance_exact':'1',
        'source_reduced_trace_distance_exact':'0',
        'all_classical_native_history_laws_equal_for_diagonal_source':True,
        'source_localization_from_memory_charge_proved':False,
        'transported_source_density':'n_v Q/2 transports the total named charge but is a chosen local-density assignment; individual half-density eigenvalues do not establish two independently protected daughter charges.',
        'basis_memory_dephasing_character_distribution_exact':list(map(str,probabilities)),
        'basis_memory_reset_character_distribution_exact':list(map(str,probabilities)),
        'exact_preservation_condition':'Every admitted memory, readout and apparatus Kraus operator must commute with I tensor X6.',
        'energy_barrier_lifetime_or_mass_derived':False}


def run():
    loops = [loop((0,11), [((0,5),(5,11)),((0,11),(0,5)),((5,11),(0,11))]),
             loop((0,1), [((0,5),(1,5)),((5,11),(1,11)),((0,11),(0,1))])]
    assert [x['voltage_mod6']for x in loops]==[3,2]
    assert all(x['probability_exact']=='1/8'for x in loops)
    return {'verdict':'PASS_ACTIVE_NATIVE_Z6_MEMORY_AND_NORMALIZATION_BOUNDARY',
        'architecture':{'source':'unchanged binary two-load native C/E carrier',
            'memory_dimension':6,'voltage':'3*transported_order_parity +2*current*1{seam=(0,5),spectator=1} mod6',
            'voltage_selected_as_architecture':True,'global_physical_clock_supplied':False},
        'complete_local_instruments':local_complete_certificate(),
        'larger_native_carrier':federation_extension_certificate(),
        'actual_closed_native_loops':loops,
        'algebra':{'generated':'M66(C) tensor C[X6]','commutant':'I66 tensor C[X6]',
            'commutant_complex_dimension':6,'six_protected_character_projectors':True,
            'source_only_charge_theorem_bypassed_by_explicit_memory':True,
            'coherence_between_distinct_characters_protected':False},
        'charge_normalization_and_splitting':charge_boundaries(),
        'localization_and_noise':localization_and_memory_error(),
        'scope':{'physical_electric_charge_identified':False,'physical_quark_identified':False,
            'rest_mass_GeV':None,'OpenWave_criterion_awards':0,
            'actual_memory_service_execution':'separate native_charge_readout package; this receipt verifies exact complete target maps and structural tensor transfer'}}
