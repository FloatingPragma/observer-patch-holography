"""Exact repair algebra on the existing complete binary occupation carrier.

Certificates use actual signed partial injections, diagonal effects and token
trees; no dense4096-square matrix or numerical rank enters.  The direct-sum
carrier is the many_body component domain, not an amendment of class V2's
fixed-global-number admission.  Its extremal noiseless qubit requires coherence
between numbers zero and twelve and has distance one against local errors.
"""
from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, prod

from ..exchange_sector import many_body as native


DIMENSION = native.DIMENSION
FULL = native.FULL
MODES = native.MODES
EXTREMES = (0, FULL)
SCHEMA = 'oph.native_occupation_algebra.v1'


def _domain(law, character):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('a literal C or E law required')
    if type(character) is not int or character not in (-1, 1):
        raise ValueError('an exact exchange character +1 or -1 required')


@lru_cache(None)
def _rows(law, character):
    _domain(law, character)
    return native.branches(law, character)


def _weights(port_weights):
    if type(port_weights) is not tuple or len(port_weights) != MODES:
        raise ValueError('twelve strictly positive exact rational port weights required')
    result = []
    for weight in port_weights:
        if type(weight) is int:
            numerator, denominator = weight, 1
        elif type(weight) is F:
            numerator, denominator = weight.numerator, weight.denominator
        else:
            raise ValueError('strictly positive exact rational port weights required')
        if (type(numerator) is not int or type(denominator) is not int
                or numerator <= 0 or denominator <= 0):
            raise ValueError('valid positive rational numerators and denominators required')
        # Detach even Fraction payloads: Fraction(existing_fraction) may alias.
        result.append(F(numerator, denominator))
    return tuple(result)


def _raw_weights(port_weights):
    weights = _weights(port_weights)
    return tuple(prod((weights[p] for p in range(MODES) if mask >> p & 1), start=F(1))
                 for mask in range(DIMENSION))


def weighted_reference(port_weights=(1,)*MODES):
    """One faithful proof reference; no unique sector population is selected."""
    weights = _weights(port_weights)
    normalizer = prod((1+w for w in weights), start=F(1))
    return tuple(w/normalizer for w in _raw_weights(weights))


def _tree_seams():
    reached, selected = {0}, []
    while len(reached) < MODES:
        for index, (a, b) in enumerate(native.SEAMS):
            if (a in reached) != (b in reached):
                reached.update((a, b))
                selected.append(index)
                break
        else:
            raise ArithmeticError('the actual native seam graph is disconnected')
    return tuple(selected)


def _effect_labels(law):
    # The two oriented moving effects are exactly n_a(1-n_b), (1-n_a)n_b.
    offset = 4 if law == 'E' else 0
    return [[edge, coin, raw+offset] for edge in _tree_seams()
            for coin, raw in ((0, 1), (1, 2))]


def _signatures(rows, labels):
    lookup = {(native.SEAMS.index(r.seam), r.coin, r.raw): r for r in rows}
    signatures = [0]*DIMENSION
    for bit, label in enumerate(labels):
        if (type(label) is not list or len(label) != 3
                or any(type(v) is not int for v in label)
                or tuple(label) not in lookup):
            raise ValueError('a valid raw native effect label required')
        for mask in lookup[tuple(label)].inputs:
            signatures[mask] |= 1 << bit
    return signatures


def _inventory(rows):
    """Check the upper algebra inclusion, including the identical extremes."""
    if not rows or any(type(row) is not native.Branch for row in rows):
        raise ValueError('a complete native Branch inventory required')
    law = rows[0].law
    _domain(law, 1)
    labels = []
    for row in rows:
        if (type(row.law) is not str or row.law != law
                or type(row.seam) is not tuple or len(row.seam) != 2
                or any(type(p) is not int for p in row.seam) or row.seam not in native.SEAMS
                or type(row.coin) is not int or row.coin not in (0, 1)
                or type(row.raw) is not int or not 0 <= row.raw < (1 << native.circuit.RECORD_BITS[law])
                or type(row.abort) is not bool
                or row.abort != (row.raw not in native.circuit.record_values(law, row.coin))
                or any(type(values) is not tuple for values in (row.inputs, row.outputs, row.signs))):
            raise ValueError('malformed native raw branch or abort label')
        labels.append((row.seam, row.coin, row.raw))
    expected = {(seam, coin, raw) for seam in native.SEAMS for coin in (0, 1)
                for raw in range(1 << native.circuit.RECORD_BITS[law])}
    if len(labels) != len(expected) or set(labels) != expected:
        raise ArithmeticError('complete raw menu omitted or duplicated a slot, including zero/abort')
    extrema, incoming = [], Counter()
    for index, row in enumerate(rows):
        if len(set(row.inputs)) != len(row.inputs) or len(set(row.outputs)) != len(row.outputs):
            raise ArithmeticError('an actual branch is not a partial injection')
        if row.abort and row.inputs:
            raise ArithmeticError('a reserved abort slot acquired support')
        restricted = []
        for source, target, phase in zip(row.inputs, row.outputs, row.signs, strict=True):
            if (type(source) is not int or type(target) is not int
                    or not 0 <= source < DIMENSION or not 0 <= target < DIMENSION
                    or type(phase) is not int or phase not in (-1, 1)
                    or source.bit_count() != target.bit_count()):
                raise ArithmeticError('native transport changed its number/sign domain')
            incoming[row.seam, row.coin, source] += 1
            if source in EXTREMES:
                restricted.append((source, target, phase))
        if restricted not in ([], [(0, 0, 1), (FULL, FULL, 1)]):
            raise ArithmeticError('the two extremal representations differ')
        if restricted:
            extrema.append(index)
    if (len(incoming) != 2*len(native.SEAMS)*DIMENSION
            or set(incoming.values()) != {1}):
        raise ArithmeticError('the complete coin/raw instrument omitted a column')
    return {'raw_slots': len(rows), 'zero_slots': sum(not row.inputs for row in rows),
            'abort_slots': sum(row.abort for row in rows),
            'columns': sum(len(row.inputs) for row in rows),
            'extremal_identity_slots': extrema}


def algebra_certificate(law='C', character=-1):
    """Produce an explicit effect-separation and signed transport witness."""
    _domain(law, character)
    rows = _rows(law, character)
    labels = _effect_labels(law)
    neighbors = [dict() for _ in range(DIMENSION)]
    for index, row in enumerate(rows):
        for a, b, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            if a != b:
                neighbors[a].setdefault(b, (index, sign))
    roots, forest = [], []
    for number in range(1, MODES):
        root = (1 << number)-1
        roots.append(root)
        reached, queue = {root}, deque([root])
        while queue:
            a = queue.popleft()
            for b, (index, sign) in sorted(neighbors[a].items()):
                if b not in reached:
                    reached.add(b)
                    queue.append(b)
                    forest.append([a, b, index, sign])
        if len(reached) != comb(MODES, number):
            raise ArithmeticError('an actual number-sector token graph is disconnected')
    result = {'schema': SCHEMA, 'law': law, 'character': character,
              'effect_labels': labels, 'effect_signatures': _signatures(rows, labels),
              'roots': roots, 'transport_forest': forest, 'inventory': _inventory(rows),
              'sector_dimensions': [comb(MODES, n) for n in range(MODES+1)],
              'algebra_dimension': comb(2*MODES, MODES)-1,
              'commutant_dimension': MODES-1+4}
    verify_algebra(law, character, result)
    return result


def verify_algebra(law, character, certificate):
    """Replay supplied witnesses independently of the tree/signature search."""
    _domain(law, character)
    keys = {'schema', 'law', 'character', 'effect_labels', 'effect_signatures',
            'roots', 'transport_forest', 'inventory', 'sector_dimensions',
            'algebra_dimension', 'commutant_dimension'}
    if (type(certificate) is not dict or set(certificate) != keys
            or certificate['schema'] != SCHEMA or certificate['law'] != law
            or type(certificate['character']) is not int or certificate['character'] != character):
        raise ValueError('a complete occupation algebra certificate required')
    rows = _rows(law, character)
    labels = certificate['effect_labels']
    if (type(labels) is not list or len(labels) != 2*(MODES-1)
            or any(type(label) is not list or len(label) != 3 for label in labels)):
        raise ValueError('a bounded list of22 native effect witnesses required')
    signatures = _signatures(rows, labels)
    supplied = certificate['effect_signatures']
    if (type(supplied) is not list or any(type(v) is not int for v in supplied)
            or supplied != signatures):
        raise ArithmeticError('effect signatures differ from actual branch supports')
    classes = {}
    for mask, signature in enumerate(signatures):
        classes.setdefault(signature, []).append(mask)
    if (len(classes) != DIMENSION-1
            or [group for group in classes.values() if len(group) > 1] != [[0, FULL]]):
        raise ArithmeticError('effects fail to separate every nonextremal occupation')
    roots = certificate['roots']
    if (type(roots) is not list or len(roots) != MODES-1
            or any(type(x) is not int or not 0 < x < FULL for x in roots)
            or sorted(x.bit_count() for x in roots) != list(range(1, MODES))):
        raise ValueError('exactly one root in every intermediate number sector required')
    forest = certificate['transport_forest']
    if type(forest) is not list or len(forest) != DIMENSION-MODES-1:
        raise ValueError('a complete bounded transport forest required')
    lookups = [dict(zip(r.inputs, zip(r.outputs, r.signs))) for r in rows]
    reached = set(roots)
    for edge in forest:
        if (type(edge) is not list or len(edge) != 4
                or any(type(v) is not int for v in edge)):
            raise ValueError('an exact signed transport edge required')
        a, b, branch, sign = edge
        if (a not in reached or b in reached or not 0 < b < FULL
                or a.bit_count() != b.bit_count() or not 0 <= branch < len(rows)
                or lookups[branch].get(a) != (b, sign)):
            raise ArithmeticError('transport witness is not a new actual signed matrix unit')
        reached.add(b)
    if reached != set(range(1, FULL)):
        raise ArithmeticError('transport forest does not span all intermediate sectors')
    expected = {'inventory': _inventory(rows),
                'sector_dimensions': [comb(MODES, n) for n in range(MODES+1)],
                'algebra_dimension': comb(2*MODES, MODES)-1,
                'commutant_dimension': MODES-1+4}
    for key, value in expected.items():
        if json.dumps(certificate[key], sort_keys=True, allow_nan=False) != json.dumps(value, sort_keys=True, allow_nan=False):
            raise ArithmeticError('occupation algebra conclusion or inventory changed: '+key)
    if any(type(certificate[key]) is not int for key in ('algebra_dimension', 'commutant_dimension')):
        raise ValueError('exact integer algebra dimensions required')
    return True


def faithful_certificate(port_weights=(1,)*MODES, law='C', character=-1):
    """Exact full-column stationarity, separately for every native seam."""
    _domain(law, character)
    theta = _weights(port_weights)
    mass = _raw_weights(theta)
    _inventory(_rows(law, character))
    by_seam = {edge: [] for edge in native.SEAMS}
    for row in _rows(law, character):
        by_seam[row.seam].append(row)
    for (a, b), rows in by_seam.items():
        incoming, outgoing = [F(0)]*DIMENSION, [F(0)]*DIMENSION
        for row in rows:
            weight = theta[a if row.coin else b]/(theta[a]+theta[b])
            for source, target in zip(row.inputs, row.outputs, strict=True):
                incoming[source] += weight
                outgoing[target] += weight*mass[source]
        if incoming != [F(1)]*DIMENSION or tuple(outgoing) != mass:
            raise ArithmeticError('a native seam lost completeness or product invariance')
    normalizer = sum(mass, F(0))
    if normalizer != prod((1+w for w in theta), start=F(1)):
        raise ArithmeticError('full-number product normalization failed')
    conditional = [sum((w for mask, w in enumerate(mass) if mask.bit_count() == n), F(0))
                   for n in range(MODES+1)]
    return {'law': law, 'character': character, 'port_weights': list(map(str, theta)),
            'seams_checked': len(native.SEAMS), 'columns_checked': 2*len(native.SEAMS)*DIMENSION,
            'reference_normalizer': str(normalizer),
            'sector_normalizers': list(map(str, conditional)),
            'smallest_reference_weight': str(min(mass)/normalizer),
            'largest_reference_weight': str(max(mass)/normalizer),
            'all_seams_complete_and_faithfully_invariant': True,
            'sector_populations_or_chemical_potential_selected': False,
            'weighted_source_service_executed': False}


def peripheral_certificate(law='C', character=-1):
    """A strictly positive diagonal in the supported Kraus linear span."""
    _domain(law, character)
    _inventory(_rows(law, character))
    diagonal, labels = [0]*DIMENSION, []
    for index, row in enumerate(_rows(law, character)):
        if row.inputs and row.raw & 3 == 0:
            if row.inputs != row.outputs or any(sign != 1 for sign in row.signs):
                raise ArithmeticError('a holding branch is not its positive diagonal effect')
            labels.append(index)
            for mask in row.inputs:
                diagonal[mask] += 1
    if min(diagonal) <= 0:
        raise ArithmeticError('holding Kraus span lacks a strictly positive diagonal witness')
    return {'holding_branch_indices': labels, 'positive_diagonal': diagonal,
            'smallest_diagonal': min(diagonal), 'largest_diagonal': max(diagonal),
            'only_peripheral_eigenvalue': 1,
            'convergence_rate_claimed': False,
            'scope': 'one fixed positive all-seam selector with fixed endpoint weights; no adaptive-history mixing assertion'}


def locality_certificate():
    """Exact finite witnesses for the count-overlap proof of proper locality."""
    overlaps = []
    for size in range(1, MODES):
        # Adjacent local counts k,k+1 share total k+1 by outside counts1,0.
        overlaps.append({'core_size': size, 'outside_size': MODES-size,
                         'adjacent_count_witnesses': [[k, 1, k+1, 0] for k in range(size)]})
    return {'proper_nonempty_core_sizes': list(range(1, MODES)), 'overlap_witnesses': overlaps,
            'proper_port_local_fixed_algebra': 'C I',
            'proper_CAR_local_fixed_observable_algebra': 'C I',
            'fixed_number_source_algebra': 'C I_N for each N=0,...,12',
            'number_superselected_fixed_algebra': 'span{P_0,...,P_12}, complex dimension13',
            'moving_particles_with_global_charge_excluded': False,
            'controlled_lifetimes_or_refining_defects_excluded': False}


def extremal_certificate():
    """The complete carrier's noiseless code and its one-port error control."""
    one_port = []
    for port in range(MODES):
        if ((0 >> port) & 1, (FULL >> port) & 1) != (0, 1):
            raise ArithmeticError('extremal local occupation control failed')
        leaked = (1 << port, FULL ^ (1 << port))
        if any(mask in EXTREMES for mask in leaked):
            raise ArithmeticError('one-port bit flip did not leave the extremal code')
        one_port.append({'port': port, 'compressed_number_diagonal': [0, 1],
                         'compressed_Z_diagonal': [1, -1], 'X_leakage_masks': list(leaked)})
    return {'code_masks': list(EXTREMES), 'quantum_dimension': 2,
            'native_branch_action': 'scalar identity or zero on the entire code',
            'one_port_controls': one_port, 'local_error_code_distance': 1,
            'cat_plus_minus_proper_marginals_equal': True,
            'full_code_locally_indistinguishable': False,
            'phase_flip_probability_p_X_expectation_after_t': '(1-2*p)^t times initial expectation',
            'added_phase_or_bit_error_is_native_C_E_branch': False,
            'cross_global_number_coherence_required': True,
            'class_V2_fixed_global_number_admission_satisfied': False,
            'localized_topological_defect_constructed': False}


def certificate():
    """Compact custody of fully replayed witnesses; raw APIs stay reproducible."""
    def digest(value):
        return sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                 allow_nan=False).encode()).hexdigest()

    cases = {}
    for law in ('C', 'E'):
        for character in (1, -1):
            algebra = algebra_certificate(law, character)
            verify_algebra(law, character, algebra)
            peripheral = peripheral_certificate(law, character)
            references = {}
            for label, weights in (('fair', (1,)*MODES), ('biased_control', (1, 2)*6)):
                witness = faithful_certificate(weights, law, character)
                references[label] = {
                    'witness_sha256': digest(witness),
                    'seams_checked': witness['seams_checked'],
                    'columns_checked': witness['columns_checked'],
                    'normalizer': witness['reference_normalizer'],
                    'smallest_reference_weight': witness['smallest_reference_weight'],
                    'largest_reference_weight': witness['largest_reference_weight'],
                    'all_seams_complete_and_faithfully_invariant': True}
            inventory = algebra['inventory']
            cases[f'{law}:{character:+d}'] = {
                'algebra': {'witness_sha256': digest(algebra),
                            'effect_count': len(algebra['effect_labels']),
                            'signature_count': len(algebra['effect_signatures']),
                            'effect_class_count': len(set(algebra['effect_signatures'])),
                            'transport_roots': len(algebra['roots']),
                            'transport_edges': len(algebra['transport_forest']),
                            'raw_slots': inventory['raw_slots'],
                            'zero_slots': inventory['zero_slots'],
                            'abort_slots': inventory['abort_slots'],
                            'columns': inventory['columns'],
                            'extremal_identity_slots': len(inventory['extremal_identity_slots']),
                            'algebra_dimension': algebra['algebra_dimension'],
                            'commutant_dimension': algebra['commutant_dimension'],
                            'full_witness_replayed': True},
                'references': references,
                'peripheral': {'witness_sha256': digest(peripheral),
                               'holding_branches': len(peripheral['holding_branch_indices']),
                               'diagonal_entries': len(peripheral['positive_diagonal']),
                               'smallest_diagonal': peripheral['smallest_diagonal'],
                               'largest_diagonal': peripheral['largest_diagonal'],
                               'only_peripheral_eigenvalue': 1}}
    locality, extremal = locality_certificate(), extremal_certificate()
    local_summary = {key: value for key, value in locality.items()
                     if key not in ('overlap_witnesses', 'proper_nonempty_core_sizes')}
    local_summary.update(witness_sha256=digest(locality), proper_core_sizes_checked=MODES-1,
                         adjacent_count_witnesses=sum(len(row['adjacent_count_witnesses'])
                                                     for row in locality['overlap_witnesses']))
    extremal_summary = {key: value for key, value in extremal.items() if key != 'one_port_controls'}
    extremal_summary.update(witness_sha256=digest(extremal), one_port_controls_checked=MODES)
    return {'schema': 'oph.native_occupation_catalogue.v1', 'cases': cases,
            'sector_dimensions': [comb(MODES, n) for n in range(MODES+1)],
            'effect_tree_seam_indices': list(_tree_seams()),
            'effect_labels': {law: [':'.join(map(str, label)) for label in _effect_labels(law)]
                              for law in ('C', 'E')},
            'effect_label_format': 'native seam index:coin:raw',
            'fair_port_weights': 'twelve weights equal1',
            'biased_control_port_weights': '(1,2) repeated six times',
            'branch_algebra': '(direct_sum N=1..11 M_binomial(12,N)) direct_sum C I_extremal',
            'fixed_algebra': 'C^11 direct_sum M2_extremal',
            'stationary_states': 'arbitrary positive extremal2x2 block plus sum N=1..11 p_N*pi_N; total trace1',
            'fixed_selector_Schrodinger_limit': 'P_ext*rho*P_ext + sum N=1..11 Tr(P_N*rho)*pi_N',
            'pi_N': 'diagonal product weights conditioned on total N; sector populations are initial data',
            'maximum_passive_quantum_code_dimension': 2,
            'nontrivial_passive_quantum_code_only_on_extremal_multiplicity': True,
            'fixed_number_quantum_multiplicity_dimension': 1,
            'locality': local_summary, 'extremal': extremal_summary,
            'scope': {'complete_binary_component': True, 'class_V2_carrier_revised': False,
                      'all_source_quantum_observables_classified': True,
                      'joint_controller_history_charges_classified': False,
                      'active_recovery_codes_classified': False,
                      'refinement_or_physical_particle_identification': False,
                      'convergence_rate_or_physical_clock_calibration': False}}
