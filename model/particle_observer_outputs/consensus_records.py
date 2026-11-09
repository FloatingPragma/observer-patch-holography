"""Finite detector records, observer meanings and complete comparisons.

The preparation, apparatus and semantic dictionaries are declared program
data. Orthogonal records can be copied; an unknown quantum state is not
cloned. A shared interpretation pushes forward the actual instrument law.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product
import json

import sympy as s

from model.particle_atlas_horizontal import measurement
from . import readout_contract


SOURCE_INPUTS = (
    'model/particle_atlas_horizontal/measurement.py',
    'model/particle_observer_outputs/readout_contract.py',
    'Lean/EventAlgebra/OperationalBornBridge.lean',
    'Lean/EventAlgebra/PublicRecordAlgebra.lean',
    'Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean',
    'Lean/Geometry/RecordGluingPrinciple.lean',
    'Lean/QFT/ObserverNetDescent.lean',
    'Lean/QFT/OperationalOverlapEvidence.lean',
)
FAILURE = 4
ALICE_ENCODING = (0, 1, 2, 3)
BOB_ENCODING = (2, 0, 3, 1)


def _clean(value):
    if isinstance(value, s.MatrixBase):
        return s.ImmutableSparseMatrix(value.applyfunc(s.simplify))
    return s.simplify(value)


def menu(setting):
    """Reuse the exact owning four-bin preparation and projective apparatus."""
    _, psi, position, momentum = measurement.construction()
    if setting not in ('position', 'momentum'):
        raise ValueError('choose one declared apparatus setting')
    projectors = position if setting == 'position' else momentum
    result = {str(j): p for j, p in enumerate(projectors)}
    readout_contract.verify_projectors(result)
    return {'state': _clean(psi*psi.H), 'psi': psi, 'projectors': result}


def _encoding(encoding):
    encoding = tuple(encoding)
    if len(encoding) != 4 or set(encoding) != set(range(4)):
        raise ValueError('an injective four-event record encoding is required')
    return encoding


def dictionary(encoding):
    """Map an observer's local codes back to the common raw event number."""
    encoding = _encoding(encoding)
    return {code: event for event, code in enumerate(encoding)} | {FAILURE: None}


def copy_unitary(encoding):
    """All 20 inputs: |j,a> -> |j,(a+encode(j)) mod5>."""
    encoding = _encoding(encoding)
    return s.ImmutableSparseMatrix(20, 20,
        {(5*j+(a+encoding[j]) % 5, 5*j+a): 1 for j in range(4) for a in range(5)})


def pointer_unitary(setting):
    """The inherited complete pointer coupling, including nonblank inputs."""
    blocks = readout_contract.pointer_blocks(menu(setting)['projectors'])
    return s.ImmutableSparseMatrix(s.BlockMatrix(blocks).as_explicit())


def record_isometry(setting, alice=ALICE_ENCODING, bob=BOB_ENCODING):
    """Field -> pointer x Alice x Bob x field, retaining all 400 output rows.

    The input pointer and observer registers are specified blanks. Their
    complete operations are the pointer and controlled-copy unitaries.
    """
    alice, bob = _encoding(alice), _encoding(bob)
    entries = {}
    for j, projection in enumerate(menu(setting)['projectors'].values()):
        base = ((j*5+alice[j])*5+bob[j])*4
        for (field_out, field_in), value in projection.todok().items():
            entries[base+field_out, field_in] = value
    return s.ImmutableSparseMatrix(400, 4, entries)


def read_channel(encoding, failure_probability=0):
    """Declared independent classical read service; failure is retained."""
    encoding = _encoding(encoding)
    probability = s.sympify(failure_probability)
    if probability.has(s.Float) or probability.is_real is not True:
        raise ValueError('an exact real failure probability is required')
    if probability.is_nonnegative is not True or (1-probability).is_nonnegative is not True:
        raise ValueError('the failure probability must lie in [0,1]')
    return s.ImmutableMatrix(5, 4, lambda code, event:
        probability if code == FAILURE else
        1-probability if code == encoding[event] else 0)


def _channel(channel):
    channel = s.ImmutableMatrix(channel)
    if channel.shape != (5, 4) or any(x.has(s.Float) or x.is_nonnegative is not True for x in channel):
        raise ValueError('an exact nonnegative five-record/four-event channel is required')
    if any(sum(channel[:, j]) != 1 for j in range(4)):
        raise ValueError('each event must retain all original outcome weight')
    return channel


def joint_kraus(setting, alice_channel=None, bob_channel=None):
    """Every raw event/local-record triple; all 100 branches are retained."""
    a = _channel(read_channel(ALICE_ENCODING) if alice_channel is None else alice_channel)
    b = _channel(read_channel(BOB_ENCODING) if bob_channel is None else bob_channel)
    projections = tuple(menu(setting)['projectors'].values())
    return {(j, x, y): _clean(s.sqrt(a[x, j]*b[y, j])*projections[j])
            for j, x, y in product(range(4), range(5), range(5))}


def joint_law(setting, rho=None, alice_channel=None, bob_channel=None,
              alice_dictionary=None, bob_dictionary=None):
    """Original probabilities and interpretations; no agreement conditioning."""
    data = menu(setting)
    rho = data['state'] if rho is None else s.ImmutableMatrix(rho)
    if rho.shape != (4, 4) or _clean(s.trace(rho)) != 1 or rho != rho.H:
        raise ValueError('a trace-one Hermitian state on this four-bin source is required')
    if any(x.is_nonnegative is not True for x in rho.eigenvals()):
        raise ValueError('a positive source state is required')
    ad = dictionary(ALICE_ENCODING) if alice_dictionary is None else dict(alice_dictionary)
    bd = dictionary(BOB_ENCODING) if bob_dictionary is None else dict(bob_dictionary)
    if set(ad) != set(range(5)) or set(bd) != set(range(5)):
        raise ValueError('each observer must interpret every record, including failures')
    rows = []
    for (j, a, b), k in joint_kraus(setting, alice_channel, bob_channel).items():
        p = _clean(s.trace(k*rho*k.H))
        meaning_a, meaning_b = ad[a], bd[b]
        status = ('failure' if meaning_a is None or meaning_b is None else
                  'agreement' if meaning_a == meaning_b else 'disagreement')
        rows.append({'setting': setting, 'raw_event': j, 'alice_record': a, 'bob_record': b,
                     'alice_meaning': meaning_a, 'bob_meaning': meaning_b,
                     'comparison': status, 'probability': p})
    if _clean(sum(row['probability'] for row in rows)) != 1:
        raise ArithmeticError('the complete observer comparison lost event probability')
    return rows


def marginal(rows, keys):
    keys = (keys,) if isinstance(keys, str) else tuple(keys)
    result = {}
    for row in rows:
        key = tuple(row[k] for k in keys)
        result[key] = _clean(result.get(key, 0)+row['probability'])
    return result


def pushforward(probabilities, meanings):
    """A coarse interpretation sums event weights; it does not select them."""
    if len(probabilities) != len(meanings):
        raise ValueError('every original event requires a meaning')
    result = {}
    for probability, meaning in zip(probabilities, meanings, strict=True):
        result[meaning] = _clean(result.get(meaning, 0)+probability)
    return result


def glue_sections(events, sections):
    """Finite function-algebra descent on an actual covering family.

    The caller supplies local event functions, not probability marginals.
    Pairwise agreement plus literal event coverage produces the unique
    global function. Compatible marginal laws alone need not have this
    property.
    """
    events = set(events)
    result = {}
    if not events or not sections:
        raise ValueError('a nonempty event set and covering family are required')
    for section in sections:
        if not set(section) <= events:
            raise ValueError('an observer refers to an event outside the declared universe')
        for event, value in section.items():
            if event in result and result[event] != value:
                raise ValueError('local meanings disagree on their shared event')
            result[event] = value
    if set(result) != events:
        raise ValueError('the local records do not cover every declared event')
    return result


def coarse_error_bound(source_law, comparison_law, meanings, event_error=F(0), convention_error=F(0)):
    """Complete finite-event data processing and a declared approximation budget."""
    source_law, comparison_law = tuple(map(s.sympify, source_law)), tuple(map(s.sympify, comparison_law))
    if len(source_law) != len(comparison_law) or len(meanings) != len(source_law):
        raise ValueError('the comparison must use the same complete event alphabet')
    for law in (source_law, comparison_law):
        if sum(law) != 1 or any(p.is_nonnegative is not True for p in law):
            raise ValueError('complete normalized nonnegative event laws are required')
    event_error, convention_error = F(event_error), F(convention_error)
    if min(event_error, convention_error) < 0:
        raise ValueError('error budgets must be nonnegative')
    p, q = pushforward(source_law, meanings), pushforward(comparison_law, meanings)
    fine = _clean(sum(abs(a-b) for a, b in zip(source_law, comparison_law, strict=True))/2)
    coarse = _clean(sum(abs(p[key]-q[key]) for key in p)/2)
    if _clean(fine-coarse).is_nonnegative is not True:
        raise ArithmeticError('deterministic interpretation increased total variation')
    return {'complete_event_TV': fine, 'coarse_record_TV': coarse,
            'all_errors_upper': min(s.S.One, fine+s.Rational(event_error)+s.Rational(convention_error)),
            'continuum_interpretation_requires_separate_domain_and_approximation_bound': True}


def certificate():
    ideal = {setting: joint_law(setting) for setting in ('position', 'momentum')}
    wrong = dictionary(BOB_ENCODING)
    wrong[BOB_ENCODING[0]], wrong[BOB_ENCODING[1]] = 1, 0
    wrong_rows = joint_law('momentum', bob_dictionary=wrong)
    noisy = joint_law('momentum', alice_channel=read_channel(ALICE_ENCODING, s.Rational(1, 20)),
                     bob_channel=read_channel(BOB_ENCODING, s.Rational(1, 10)))
    complete = {}
    for setting, rows in ideal.items():
        data = menu(setting)
        original = measurement.outcome_probabilities(data['state'], tuple(data['projectors'].values()))
        complete[setting] = {
            'original_event_probabilities': list(map(str, original)),
            'all_100_joint_record_probabilities': [dict(row, probability=str(row['probability'])) for row in rows],
            'comparison_probabilities': {k[0]: str(v) for k, v in marginal(rows, 'comparison').items()},
            'same_apparatus_repeat_records': measurement.sequential_records(
                data['state'], tuple(data['projectors'].values()), tuple(data['projectors'].values()))}
    glued = glue_sections(range(4), [{0: 'bin 0', 1: 'bin 1', 2: 'bin 2'},
                                   {1: 'bin 1', 2: 'bin 2', 3: 'bin 3'}])
    return {'schema': 'oph.observer_outputs.consensus_records.v1',
        'source_inputs': list(SOURCE_INPUTS),
        'source_preparation': '(1,1,0,0)/sqrt(2)',
        'apparatus_choices': ['direct four-bin detector', 'F4 phase analyzer followed by the same detector'],
        'event_key_includes_apparatus_setting': True,
        'alice_encoding': list(ALICE_ENCODING), 'bob_encoding': list(BOB_ENCODING),
        'complete_ideal_experiments': complete,
        'wrong_dictionary_comparison': {k[0]: str(v) for k, v in marginal(wrong_rows, 'comparison').items()},
        'correct_dictionary_comparison': {k[0]: str(v) for k, v in marginal(ideal['momentum'], 'comparison').items()},
        'retained_failure_comparison': {k[0]: str(v) for k, v in marginal(noisy, 'comparison').items()},
        'noninjective_meaning': {'map': ['even', 'odd', 'even', 'odd'],
            'probabilities': {k: str(v) for k, v in pushforward(
                (s.Rational(1, 2), s.Rational(1, 4), 0, s.Rational(1, 4)),
                ('even', 'odd', 'even', 'odd')).items()},
            'original_bin_not_recoverable_from_parity_alone': True},
        'finite_cover_glued_section': glued,
        'graphic_nodes': [
            {'id': 'finite_process', 'title': 'A finite process runs',
             'text': 'The chosen state and apparatus determine a complete instrument.'},
            {'id': 'orthogonal_record', 'title': 'A pointer writes a record',
             'text': 'An individual record is one event; the histogram follows Tr(rho E).'},
            {'id': 'local_meanings', 'title': 'Observers interpret their copies',
             'text': 'Different local codes refer to the same event through explicit dictionaries.'},
            {'id': 'comparison', 'title': 'They compare every record',
             'text': 'Agreement, mismatch and failure retain their original weights.'},
            {'id': 'effective_description', 'title': 'A field describes the agreed statistics',
             'text': 'Position, momentum and mass require a tested observable map and approximation domain.'}],
        'physical_interpretation': {
            'one_event_is_not_the_whole_probability_distribution': True,
            'Born_rule_is_unchanged_by_copying_or_dictionary_repair': True,
            'orthogonal_record_copy_is_not_unknown_quantum_state_cloning': True,
            'all_quantum_references_retained': True,
            'record_agreement_proves_three_plus_one_continuum': False,
            'rest_mass_relation': 'E0=m c^2 only after the same effective sector has relativistic rest-energy and unit conventions',
            'simulation_generates_records_not_an_additional_objective_continuum': True}}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
