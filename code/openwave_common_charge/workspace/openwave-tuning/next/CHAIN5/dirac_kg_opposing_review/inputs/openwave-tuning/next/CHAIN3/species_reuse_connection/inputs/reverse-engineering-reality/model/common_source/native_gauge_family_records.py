"""Complete current-mediated records for a declared geometric family protocol.

The trace frame is lifted from the actual native port currents.  With
``C = sum(F_a.H F_a)``, the thirteen Kraus operators are
``M_0 = cos(t sqrt(C))`` and ``M_a = i F_a sin(t sqrt(C))/sqrt(C)``.
Their Gram sum is identity on the whole selected fifteen-state matter space.
The twelve jump indices transform in the adjoint representation; measuring
one fixed adjoint basis is not an invariant classical species record.  The
coarse jump/no-jump instrument is gauge covariant and keeps every outcome.

The declared default protocol is t=pi/5, the half angle of an
order-five rotation.  The complete small catalog also retains pi/n for every
nonidentity rotation order in the native geometry; no best entry is selected.
This is a declared catalog, not a prospective prediction or a certificate of
target-independent design history.
A native C event advances family memory only with a
left-current jump; a native E event only with a right-current jump.  Left
and right supports are recovered from both components of each native Higgs
line, without a species-to-Casimir table.  Their jump effects are scalar, so
the two binary records are independent even on an entangled state supported
on the corresponding left/right tensor product.

The interaction, its duration, two separate input legs, and its feedback into
family memory are declared architecture, not consequences of unchanged C/E
repair.  The producer loads no empirical mass packet; this is runtime
separation, not a claim about design ancestry.  The stationary coefficients are
not physical pole masses, and no implementation by native source words or
detector history is asserted.  No cycle affinity is identified with a U(1)
phase.
"""
from __future__ import annotations

from functools import lru_cache
import json

import numpy as np
from scipy.linalg import block_diag, expm
import sympy as s

from . import native_mass_dressing as dressing
from . import native_matter_attachment as matter
from . import native_repair_mass as repair
from . import native_gauge_object as objects


PROTOCOL_ANGLE = s.pi / 5
ADMISSION = 2e-11


def _zero(matrix):
    return all(s.expand(x) == 0 or s.simplify(x) == 0 for x in matrix)


def _angle(value):
    value = s.sympify(value)
    if value.free_symbols or value.is_real is not True or value.is_finite is not True:
        raise ValueError('a finite real interaction angle is required')
    return value


def protocol_catalog():
    """Every native nonidentity rotation order supplies one declared half angle.

    This is a declared geometry catalog, not a prospective prediction or a
    target-independent design certificate.  The producer loads no empirical
    mass packet and performs no empirical ranking of its entries.  The
    default is the smallest angle, pi/5; these are model interaction angles,
    not calibrated durations in an observer clock.
    """
    orders = sorted({objects.permutation_order(g) for g in objects.geometry()['rotations']}
                    - {1}, reverse=True)
    return {f'half_order_{order}_rotation': s.pi / order for order in orders}


@lru_cache(maxsize=1)
def current_certificate():
    """Exact same-current identities underlying the functional-calculus proof.

    Centrality implies that every spectral function of C commutes with every
    F_a.  Consequently sum(M_a.H M_a)=sin(t sqrt(C))^2, including the zero
    eigenspace by the continuous value sin(t sqrt(c))/sqrt(c)=t.
    Adding M_0.H M_0 gives identity, for every finite real t.  The source
    adjoint intertwiner transports all twelve jump outcomes together.
    """
    data = dressing.current_data()
    frame, casimir = data['matter_frame'], data['casimir15']
    gram = sum((g.H * g for g in frame), s.zeros(15))
    if not _zero(gram - casimir):
        raise ArithmeticError('complete native current Gram identity failed')
    if any(not _zero(casimir * g - g * casimir) for g in
           (*frame, *data['representation']['matter_currents'])):
        raise ArithmeticError('native trace Casimir is not central')
    if not casimir.is_diagonal() or any(x.is_nonnegative is not True for x in casimir.diagonal()):
        raise ArithmeticError('native matter Casimir is not positive diagonal')
    return data


def validate_kraus(kraus):
    """Admit the entire finite instrument, never a conditioned jump branch."""
    kraus = np.asarray(kraus)
    if (kraus.shape != (13, 15, 15) or kraus.dtype.kind != 'c'
            or not np.isfinite(kraus).all()):
        raise ValueError('all thirteen finite complex15-by15 Kraus operators are required')
    gram = sum((k.conj().T @ k for k in kraus), np.zeros((15, 15), complex))
    error = float(np.linalg.norm(gram - np.eye(15), 2))
    if error > ADMISSION:
        raise ValueError('complete current instrument is not trace preserving')
    return error


def instrument(angle=PROTOCOL_ANGLE):
    """All thirteen Kraus operators on all fifteen input matter states."""
    angle = _angle(angle)
    data = current_certificate()
    c = np.asarray(data['casimir15'].diagonal(), float).reshape(15)
    root = np.sqrt(c)
    duration = float(angle)
    factor = np.full_like(root, duration)
    np.divide(np.sin(duration * root), root, out=factor, where=root > 0)
    zero = np.diag(np.cos(duration * root)).astype(complex)
    jump = [1j * np.asarray(g, complex) @ np.diag(factor) for g in data['matter_frame']]
    kraus = np.stack([zero, *jump])
    error = validate_kraus(kraus)
    return {'angle': angle, 'kraus': kraus, 'completeness_error': error,
            'labels': ('no_jump', *(f'current_{a}' for a in range(12))),
            'casimir_exact': data['casimir15']}


def covariance_errors(kraus=None, angle=PROTOCOL_ANGLE):
    """Dilation covariance for all twelve finite native current actions.

    With V formed by stacking all M_a and O_j the real adjoint action,
    (diag(1,O_j) tensor U_j) V = V U_j.  Tracing the jump index yields
    the covariant binary instrument.  A fixed individual jump label is not
    claimed to be invariant under continuous adjoint rotations.
    """
    data = current_certificate()
    if kraus is None:
        kraus = instrument(angle)['kraus']
    validate_kraus(kraus)
    V = np.asarray(kraus).reshape(195, 15)
    errors = []
    for rho, adjoint in zip(data['representation']['matter_currents'], data['adjoints'], strict=True):
        U = expm(float(PROTOCOL_ANGLE) * np.asarray(rho, complex))
        O = expm(float(PROTOCOL_ANGLE) * np.asarray(adjoint, float))
        errors.append(float(np.linalg.norm(np.kron(block_diag(1, O), U) @ V - V @ U, 2)))
    return np.asarray(errors)


@lru_cache(maxsize=1)
def higgs_supports():
    """Recover full invariant doublet/singlet supports from both Higgs components.

    A neutral Higgs background alone selects only half a weak doublet and is
    not invariant under the full current.  Both tensor components are used
    here; the source weak Casimir orients the two legs, and every projector
    is checked against every actual native current.
    """
    data = current_certificate()
    rep = data['representation']
    selected = list(rep['selected'])
    weak = rep['operators']['weak_C2'].extract(selected, selected)
    result = []
    for line in matter.invariant_higgs_lines():
        support = {i for component in line['components'] for (i, _), value in component.todok().items()
                   if value != 0}
        left = sorted(i for i in support if weak[i, i] != 0)
        right = sorted(i for i in support if weak[i, i] == 0)
        if not left or not right:
            raise ArithmeticError('a native Higgs line lost one of its two legs')
        row = {'origin': line['origin'], 'conjugate_higgs': line['conjugate'],
               'components': line['components']}
        for name, indices in (('left', left), ('right', right)):
            P = s.diag(*(int(i in indices) for i in range(15)))
            c = s.simplify(s.trace(P * data['casimir15']) / len(indices))
            if (not _zero(P * P - P) or not _zero(data['casimir15'] * P - c * P)
                    or any(not _zero(P * g - g * P) for g in rep['matter_currents'])):
                raise ArithmeticError('Higgs support is not an invariant scalar-Casimir domain')
            row[name + '_projector'] = s.ImmutableMatrix(P)
            row[name + '_dimension'] = len(indices)
            row[name + '_casimir'] = c
        result.append(row)
    return tuple(result)


def native_event_rates():
    """Derive both native probabilities and independence from all128 raw labels."""
    records = repair.native_records()
    p = sum(r[-1] for r in records if r[0] == 0 and r[1] == 1)
    q = sum(r[-1] for r in records if r[2] == 0 and r[3] in (4, 5))
    joint = sum(r[-1] for r in records if r[0] == 0 and r[1] == 1
                and r[2] == 0 and r[3] in (4, 5))
    if len(records) != 128 or sum(r[-1] for r in records) != 1 or joint != p * q:
        raise ArithmeticError('complete native C/E event factorization failed')
    return p, q


def sector_gates(angle=PROTOCOL_ANGLE):
    """One declared current-gated law per Higgs line, with all four binary outcomes."""
    angle = _angle(angle)
    p_native, q_native = native_event_rates()
    result = []
    for support in higgs_supports():
        left = s.sin(angle * s.sqrt(support['left_casimir'])) ** 2
        right = s.sin(angle * s.sqrt(support['right_casimir'])) ** 2
        row = dict(support)
        row.update(angle=angle, left_jump=left, right_jump=right,
                   p=p_native * left, q=q_native * right,
                   binary_probabilities=((1-left)*(1-right), (1-left)*right,
                                         left*(1-right), left*right))
        result.append(row)
    return tuple(result)


def complete_binary_records(origin, angle=PROTOCOL_ANGLE):
    """All512 native-label/binary-current records, including zero-weight labels.

    The two binary labels coarse-grain the complete169 current-outcome pairs.
    No binary outcome is postselected: failed gates reset family memory.
    Tuple entries are (cc,cr,ec,er,left_jump,right_jump,C_event,E_event,weight).
    """
    matches = [row for row in sector_gates(angle) if tuple(row['origin']) == tuple(origin)]
    if len(matches) != 1:
        raise ValueError('one native Higgs charge origin is required')
    row = matches[0]
    probabilities = row['binary_probabilities']
    result = []
    for cc, cr, ec, er, _, weight in repair.native_records():
        for left in (False, True):
            for right in (False, True):
                result.append((cc, cr, ec, er, left, right,
                               cc == 0 and cr == 1 and left,
                               ec == 0 and er in (4, 5) and right,
                               weight * probabilities[2 * int(left) + int(right)]))
    return tuple(result)


def family_models(angle=PROTOCOL_ANGLE):
    """Feed source-derived rates into the same fixed family-memory rule.

    This handoff uses reset2 and the same unconjugated native triple in every
    channel.  Rates are exactly specified analytic expressions; numerical
    family matrices remain numerical evaluations, not exact certificates.
    """
    from . import native_family_record_memory as memory
    result = []
    for row in sector_gates(angle):
        family = memory.architecture_from_factors(row['p'], row['q'], reset=2, conjugate=False)
        result.append({'gate': row, 'family': family})
    return tuple(result)


def candidate_channels(angle=PROTOCOL_ANGLE):
    """Numerical field handoff, with names matched by the exact native intertwiner.

    Naming happens after the rates and matrices have been constructed.  The
    dictionary changes no channel's protocol, reset, orientation, or scale.
    """
    _, matches = matter.action_intertwiner()
    names = {tuple(row['charge_origin']): row['existing_action_channel'] for row in matches}
    result = {}
    for model in family_models(angle):
        row, family = model['gate'], model['family']
        origin = tuple(map(str, row['origin']))
        name = names[origin]
        if name in result:
            raise ArithmeticError('native Higgs line has a duplicate action channel')
        result[name] = {'p': float(row['p']), 'q': float(row['q']),
                        'p_exact': str(row['p']), 'q_exact': str(row['q']),
                        'reset': 2, 'conjugate': False, 'angle': float(row['angle']),
                        'charge_origin': origin, 'Y': np.asarray(family['Y'], float),
                        'Higgs_components': np.asarray(row['components'], complex),
                        'left_projector': np.asarray(row['left_projector'], float),
                        'right_projector': np.asarray(row['right_projector'], float),
                        'left_casimir_exact': str(row['left_casimir']),
                        'right_casimir_exact': str(row['right_casimir'])}
    if len(result) != len(higgs_supports()):
        raise ArithmeticError('native Higgs/action matching lost a channel')
    return result


def _family_spectrum(family):
    """Keep exact geometric zero directions and undefined hierarchy ratios.

    Positive stationary weights make the range of their projector average
    equal to the span of the supported projector ranges.  Its rank is checked
    in exact native geometry, without deciding rank from floating eigenvalues.
    """
    support = [P for weight, P in zip(family['stationary'], family['projectors'], strict=True)
               if weight != 0]
    if any(float(weight) < 0 for weight in family['stationary']):
        raise ArithmeticError('family law has a negative stationary weight')
    rank = sum(support, s.zeros(3)).rank()
    values = np.linalg.eigvalsh(np.asarray(family['Y'], float))
    if values.min() < -ADMISSION or (rank == 3 and values[0] <= 0):
        raise ArithmeticError('family spectrum is inconsistent with its positive geometric rank')
    values[:3-rank] = 0
    ratios = (values/values[0]).tolist() if rank == 3 else None
    return {'family_rank_exact': rank, 'family_eigenvalues': values.tolist(),
            'family_eigenvalue_ratios': ratios,
            'hierarchy_ratio_status': 'defined' if rank == 3 else 'undefined_exact_zero_direction'}


@lru_cache(maxsize=1)
def higgs_tensor_certificate():
    """Exact neutral tensor normalizations; color counts copies, not masses.

    Each directed Higgs component is a partial isometry with initial space R.
    Its nonzero singular values are one.  The symmetric neutral component has
    twice rank(R) unit singular values.  The dictionary into the existing
    action changes signs only.  Thus a common canonical kinetic form and a
    neutral Higgs norm v give v times the family singular values, with no
    additional sector-dependent tensor modulus.
    """
    _, mapping = matter.action_intertwiner()
    matches = {tuple(row['charge_origin']): row for row in mapping}
    result = []
    for row in higgs_supports():
        origin = tuple(map(str, row['origin']))
        match = matches[origin]
        L, R = row['left_projector'], row['right_projector']
        for component in row['components']:
            D = L * component * R
            if not _zero(D.H * D - R) or not _zero(D * D.H * D - D):
                raise ArithmeticError('native Higgs tensor is not a unit partial isometry')
        neutral = row['components'][1]
        gram = neutral.H * neutral
        rank = int(s.trace(gram))
        scale = s.sympify(match['same_current_tensor_scale_exact'])
        if (not _zero(gram * gram - gram) or rank != match['neutral_rank']
                or abs(scale) != 1):
            raise ArithmeticError('native neutral tensor normalization or dictionary modulus changed')
        result.append({'charge_origin': origin, 'action_channel': match['existing_action_channel'],
                       'neutral_rank': rank,
                       'directed_nonzero_singular_value_exact': '1',
                       'directed_nonzero_singular_multiplicity': row['right_dimension'],
                       'neutral_nonzero_singular_value_exact': '1',
                       'neutral_Hilbert_Schmidt_norm_squared_exact': str(s.trace(gram)),
                       'component_nonzero_values_exact': [sorted(set(map(str, C.todok().values())))
                                                          for C in row['components']],
                       'dictionary_scale_exact': str(scale), 'dictionary_scale_is_sign_only': True,
                       'canonical_common_kinetic_form': 'I45', 'neutral_Higgs_norm': '1',
                       'relative_sector_coupling_selected_by_covariance': False})
    return tuple(result)


@lru_cache(maxsize=1)
def duty_coefficients():
    """All four unconditional stationary pulse generators at the fixed pi/5.

    For every old memory m and every one of the512 complete records r, keep
    K_(r,m)=sqrt(w_r)|F(m,r)><m| tensor exp(-i h chi_D(r) X tensor P_F).
    Disabled pulses use identity; their records remain in the channel.  The
    first-order stationary mean is Z_D=sum(pi_m w_r chi_D(r) P_F(m,r)).
    Every old memory input is retained.  Aggregating equal C/E flags is exact
    because destination and pulse policy depend only on those flags.

    Tr(Z_D) is 1,p,q,pq for always,C,E,CE.  No division by the duty factor is
    performed.  Z_C=Y-(1-p)P_reset, generally different from pY.  These are
    stationary first-order generators, not exact finite-step unitary channels
    or interacting particle poles.
    """
    from . import native_family_record_memory as memory
    policies = ('always', 'C', 'E', 'CE')
    output = {policy: [] for policy in policies}
    names = {row['charge_origin']: row['action_channel'] for row in higgs_tensor_certificate()}
    for model in family_models():
        gate, family = model['gate'], model['family']
        records = complete_binary_records(gate['origin'])
        if len(records) != 512 or len(family['stationary']) != 12:
            raise ArithmeticError('complete duty domain changed')
        flags = ((False, False), (False, True), (True, False), (True, True))
        weights = {flag: s.Integer(0) for flag in flags}
        for record in records:
            weights[(record[-3], record[-2])] += record[-1]
        weights = {flag: s.cancel(weight) for flag, weight in weights.items()}
        p, q = gate['p'], gate['q']
        expected = ((1-p)*(1-q), (1-p)*q, p*(1-q), p*q)
        if any(s.cancel(weights[flag]-value) != 0 for flag, value in zip(flags, expected, strict=True)):
            raise ArithmeticError('complete current/native record aggregation lost a probability')
        traces = {'always': s.Integer(1), 'C': p, 'E': q, 'CE': p*q}
        for policy in policies:
            projector_weights = [s.Integer(0)] * 12
            for old in range(12):
                for c, e in flags:
                    enabled = (policy == 'always' or policy == 'C' and c
                               or policy == 'E' and e or policy == 'CE' and c and e)
                    if enabled:
                        new = memory._destination(old, c, e, family['reset'], family['gC'], family['gE'])
                        projector_weights[new] += family['stationary'][old] * weights[(c, e)]
            projector_weights = tuple(s.cancel(weight) for weight in projector_weights)
            Z = s.ImmutableMatrix(sum((weight*P for weight, P in
                                      zip(projector_weights, family['projectors'], strict=True)), s.zeros(3)))
            if s.cancel(sum(projector_weights)-traces[policy]) != 0:
                raise ArithmeticError('unconditional duty trace identity failed')
            spectrum = _family_spectrum({'stationary': projector_weights,
                                          'projectors': family['projectors'], 'Y': Z})
            origin = tuple(map(str, gate['origin']))
            output[policy].append({'charge_origin': origin, 'action_channel': names[origin],
                                   'coefficient': Z, 'trace': traces[policy],
                                   'projector_weights': projector_weights,
                                   'supported_memory_directions': tuple(i for i, w in enumerate(projector_weights) if w != 0),
                                   'raw_record_count': 512, 'old_memory_input_count': 12,
                                   'disabled_pulses_are_identity_and_records_retained': True,
                                   **spectrum})
    return {policy: tuple(rows) for policy, rows in output.items()}


def duty_catalog():
    """Compact complete duty receipt, with a common-scale separation bound.

    For any positive3-by3 Z with trace r, r/3<=lambda_max(Z)<=r.
    Therefore the largest sector-max ratio is at most3 max(r)/min(r).
    The common physical multiplier cancels from this bound and is not derived.
    In the declared unit-Higgs free action, trace-one coefficients have equal
    sums of the three family mass parameters in every sector.  Tensor multiplicities give
    Tr|K45|=6 TrZ_up+6 TrZ_down+2 TrZ_lepton and three neutrino zero directions.
    """
    catalog = []
    for policy, rows in duty_coefficients().items():
        traces = [float(row['trace']) for row in rows]
        maxima = [row['family_eigenvalues'][-1] for row in rows]
        channels = [{'action_channel': row['action_channel'], 'trace_exact': str(row['trace']),
                     'trace_diagnostic': float(row['trace']),
                     'family_rank_exact': row['family_rank_exact'],
                     'family_eigenvalues': row['family_eigenvalues'],
                     'supported_memory_directions': list(row['supported_memory_directions'])}
                    for row in rows]
        catalog.append({'policy': policy, 'channels': channels,
                        'global_heaviest_ratio_bound_formula': '3 max_s Tr(Z_s) / min_s Tr(Z_s)',
                        'global_heaviest_ratio_bound_diagnostic': 3*max(traces)/min(traces),
                        'actual_heaviest_span_diagnostic': max(maxima)/min(maxima),
                        'pulse_duty_renormalized': False,
                        'common_physical_mass_scale_derived': False})
    return catalog


def catalog_summaries():
    """Retain every declared geometry-catalog entry, without empirical ranking."""
    rows = []
    for name, angle in protocol_catalog().items():
        channels = []
        for model in family_models(angle):
            gate = model['gate']
            channels.append({'charge_origin': list(map(str, gate['origin'])),
                             'p': float(gate['p']), 'q': float(gate['q']),
                             **_family_spectrum(model['family'])})
        rows.append({'protocol': name, 'angle_exact': str(angle), 'Higgs_lines': channels})
    return rows


def summary(*, include_family=True):
    """A bounded, non-writing receipt; no empirical reference packet is loaded."""
    result = instrument()
    covariance = covariance_errors(result['kraus'])
    rows = []
    for row in sector_gates():
        rows.append({'charge_origin': list(map(str, row['origin'])),
                     'left_dimension': row['left_dimension'], 'right_dimension': row['right_dimension'],
                     'left_casimir_exact': str(row['left_casimir']),
                     'right_casimir_exact': str(row['right_casimir']),
                     'left_jump': float(row['left_jump']), 'right_jump': float(row['right_jump']),
                     'p_exact': str(row['p']), 'q_exact': str(row['q']),
                     'p': float(row['p']), 'q': float(row['q'])})
    if include_family:
        for row, model in zip(rows, family_models(), strict=True):
            row.update(_family_spectrum(model['family']))
    return {'protocol_angle_exact': str(PROTOCOL_ANGLE),
            'catalog_angle_formulas': {name: str(angle) for name, angle in protocol_catalog().items()},
            'all_catalog_entries': catalog_summaries() if include_family else [],
            'duty_catalog': duty_catalog() if include_family else [],
            'native_Higgs_tensor_certificate': list(higgs_tensor_certificate()),
            'duration_parameter': 'dimensionless declared interaction angle; no physical clock calibration',
            'kraus_count': 13, 'matter_dimension': 15,
            'two_leg_raw_current_outcome_pairs': 169, 'coarse_records_per_Higgs_line': 512,
            'completeness_error': result['completeness_error'],
            'finite_dilation_covariance_errors': covariance.tolist(), 'Higgs_lines': rows,
            'scope': {'source_current_Casimir_and_Higgs_supports_derived': True,
                      'interaction_protocol_and_feedback_are_declared': True,
                      'individual_adjoint_classical_labels_are_gauge_invariant': False,
                      'coarse_binary_instrument_is_gauge_covariant': True,
                      'all_outcomes_retained_without_postselection': True,
                      'empirical_mass_inputs_consumed': False,
                      'prospective_prediction_or_target_independent_design_claimed': False,
                      'autonomous_source_selection_of_protocol_proved': False,
                      'compiled_native_source_words_or_detector_history': False,
                      'physical_pole_mass_prediction': False}}


if __name__ == '__main__':
    print(json.dumps(summary(), indent=2, sort_keys=True, allow_nan=False))
