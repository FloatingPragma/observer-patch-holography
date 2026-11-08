from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest

from . import menus as m, packet


@pytest.fixture(scope='module')
def compiled():
    words = m.compile_menu(packet.example())
    for pair in words:
        for word in pair:
            m.s.certify_word(word)
    return words


def test_exact_noncommuting_joint_menu_and_complete_residual():
    stages = packet.example(); data = m.ideal_menu(stages)
    e1, e2, residual = data['effects']
    assert e1 == {(0, 0): (F(16, 25), F())}
    assert e2 == {(0, 0): (F(81, 1250), F()), (0, 65): (F(27, 250), F()),
                  (65, 0): (F(27, 250), F()), (65, 65): (F(9, 50), F())}
    assert m.matrix_sum(e1, e2, residual) == m.identity()
    assert residual[(32, 32)] == m.ONE
    assert m.product(e1, e2)[(0, 65)] == (F(216, 3125), F())
    assert m.product(e2, e1).get((0, 65), m.ZERO) == m.ZERO
    assert e2 != {key: m.scale(x, F(9, 25)) for key, x in m.projector(stages[1][0]).items()}


def test_coarse_overlapping_events_have_one_four_event_refinement():
    effects = m.ideal_menu(packet.example())['effects']
    a = m.coarse_effect(effects, (0, 1)); b = m.coarse_effect(effects, (1, 2))
    assert m.matrix_sum(effects[0], effects[1]) == a
    assert m.matrix_sum(effects[1], effects[2]) == b
    assert m.coarse_effect(effects, ()) == {}
    assert m.coarse_effect(effects, (0, 1, 2)) == m.identity()
    for bad in ((0, 0), (True,), (3,), (-1,)):
        with pytest.raises(ValueError): m.coarse_effect(effects, bad)


@pytest.mark.parametrize('survival', [F(), F(1), F(2, 7)])
def test_exact_complex_filter_includes_endpoints(survival):
    vector = tuple((F(3) if i == 0 else F(), F(4) if i == 65 else F()) for i in range(66))
    data = m.ideal_menu(((vector, survival),))
    assert m.matrix_sum(*data['effects']) == m.identity()
    assert data['effects'][0].get((0, 65), m.ZERO) == (F(), -F(12, 25)*(1-survival**2))
    if survival != 1:
        assert data['residual'][(0, 65)] == (F(), -F(12, 25)*(survival-1))
    else:
        assert (0, 65) not in data['residual']


def test_eight_stage_bound_and_exact_menu_algebra():
    stages = packet.example()*4
    data = m.ideal_menu(stages)
    assert len(data['effects']) == 9 and m.matrix_sum(*data['effects']) == m.identity()
    bound = m.menu_bounds(8)
    assert bound['native_words'] == 16 and bound['native_pulses'] == 934272
    assert F() < bound['complete_halfdiamond_upper'] < F(1, 50000)
    assert F(1) < bound['complete_gain_upper'] < F(10001, 10000)
    one = m.word_bounds(m.WORD_PULSES, 24)
    two = m.word_bounds(m.WORD_PULSES, 23)
    assert m.menu_bounds(1)['complete_halfdiamond_upper'] == one['halfdiamond']+one['gain']*two['halfdiamond']
    assert m.menu_bounds(8)['complete_halfdiamond_upper'] > 8*m.menu_bounds(1)['complete_halfdiamond_upper']


def test_all_words_are_separately_admitted_without_host_expansion(compiled):
    for pair in compiled:
        for word in pair:
            assert len(word) == 58392 < 65536
            admitted = m.s.interpreter.source_word(word, 9)
            assert admitted['declared_resources']['native_pulses'] == 58392
            assert sum(event['angle'] == 0 for event in word) > 0
            result = m.s.interpreter.execute_word(admitted, np.array([[2/3]]))
            assert result['status'] == 'RESOURCE_STOP' and result['pulses_executed'] == 0
            assert result['blocks'][('RESOURCE_STOP',)][0,0] == 2/3
    with pytest.raises(ValueError):
        m.s.interpreter.source_word(compiled[0][0]+compiled[0][1], 9)


def test_complete_reference_output_matches_derived_instrument(compiled):
    rng = np.random.default_rng(781)
    columns = rng.normal(size=(66, 3))+1j*rng.normal(size=(66, 3))
    columns /= np.linalg.norm(columns)
    aborts, residual = m.diagnostic_output(compiled, columns)
    ideal_aborts, ideal_residual = packet.ideal_output(packet.example(), columns)
    assert len(aborts) == 4 and all(x.shape == (3, 3) for x in aborts)
    for actual, ideal in zip(aborts, ideal_aborts, strict=True):
        assert np.linalg.norm(actual-ideal, 'nuc') < 1e-8
    assert np.linalg.norm(residual-ideal_residual) < 1e-8
    assert abs(sum(np.trace(x).real for x in aborts)+np.linalg.norm(residual)**2-1) < 1e-8
    assert np.linalg.norm(ideal_aborts[1]) == np.linalg.norm(ideal_aborts[3]) == 0


def test_phase_correction_is_required_for_the_claimed_state(compiled):
    columns = np.zeros((66, 1), complex); columns[0, 0] = columns[65, 0] = 1/np.sqrt(2)
    uncorrected, abort = m.successful_word(compiled[0][0], columns)
    corrected, correction_abort = m.successful_word(compiled[0][1], uncorrected)
    wanted = columns.copy(); wanted[0] *= .6
    assert np.linalg.norm(corrected-wanted) < 1e-8
    assert np.linalg.norm(uncorrected-wanted) > .3
    assert abs(np.trace(abort)-.32) < 1e-8 and np.linalg.norm(correction_abort) < 1e-8


def test_complex_source_phase_and_untouched_directions():
    vector = tuple((F(3) if i == 0 else F(), F(4) if i == 65 else F()) for i in range(66))
    words = (m.filter_words(vector, F(3, 5)),)
    for word in words[0]: m.s.certify_word(word)
    columns = np.zeros((66, 2), complex)
    columns[0] = (.3, .2j); columns[65] = (.1j, -.4); columns[31] = (.2, .1j)
    aborts, residual = m.diagnostic_output(words, columns)
    ideal_aborts, ideal_residual = packet.ideal_output(((vector, F(3, 5)),), columns)
    assert np.linalg.norm(residual-ideal_residual) < 1e-8
    assert np.linalg.norm(residual[31]-columns[31]) < 1e-8
    assert np.linalg.norm(aborts[0]-ideal_aborts[0]) < 1e-8


def test_full_failure_alphabet_and_padding_rows(compiled):
    catalogue = m.outcomes(compiled)
    assert len(catalogue['words']) == 4 and catalogue['final_quantum_residual_labels'] == 1
    assert catalogue['native_tree_measurements'] == 0 and catalogue['zero_failure_rows_retained']
    for row in catalogue['words']:
        word = compiled[row['stage']][row['component']]
        assert row['processor_first_failure_labels'] == sum(6-2**len(p['targets']) for p in word)
        assert row['observed_padding_abort_flags'] == 1 and row['internal_padding_kraus_rows'] == 446
    columns = np.zeros((66, 2), complex); columns[0] = (1, 2j)
    address = next(i for i in range(512) if i not in m.s.tomo.codes())
    row = m.padding_branch(compiled, 0, 0, address, columns)
    assert row.shape == (2,) and np.allclose(row[1], 2j*row[0])
    assert abs(abs(row[0])**2-F(16, 25)) < 1e-8
    with pytest.raises(ValueError): m.padding_branch(compiled, 0, 0, m.s.tomo.codes()[0], columns)


def test_first_processor_failure_preserves_reference_and_rejects_bad_gate(compiled, monkeypatch):
    columns = np.zeros((66, 2), complex); columns[0] = (1, 2j)
    first = compiled[0][0][0]; dimension = 1 << len(first['targets'])
    clean = m.failure_branch(compiled, 0, 0, 0, dimension, columns)
    assert clean.shape == (512, 2) and np.linalg.norm(clean) < 1e-8
    original = m.s.compiler.native_pulse
    def corrupt(kind, angle):
        array = original(kind, angle).copy()
        array[dimension, :dimension] += 1/1000
        return array
    monkeypatch.setattr(m.s.compiler, 'native_pulse', corrupt)
    bad = m.failure_branch(compiled, 0, 0, 0, dimension, columns)
    assert np.linalg.norm(bad) > 1e-4 and np.allclose(bad[:, 1], 2j*bad[:, 0])
    m.s.certify_primitive.cache_clear()
    with pytest.raises(ArithmeticError): m.s.certify_primitive(first['kind'], first['angle'])
    m.s.certify_primitive.cache_clear()


@pytest.mark.parametrize('bad', [(), [], True, [packet.example()[0]]*9,
                                [(packet.example()[0][0], True)],
                                [(packet.example()[0][0], F(2))],
                                [([(F(),F())]*66, F(1,2))]])
def test_fail_closed_menu_domain(bad):
    with pytest.raises(ValueError): m.specification(bad)


@pytest.mark.parametrize('bad', [True, 0, 9, F(1), 1.0])
def test_integer_menu_budgets_reject_aliases(bad):
    with pytest.raises(ValueError): m.menu_bounds(bad)


def test_dense_exact_target_compiles_without_source_normalization_guess():
    vector = tuple((F(i%3-1, 7), F(i%5-2, 11)) for i in range(66))
    word = m.phase_word(vector, F(1, 2), 0)
    assert len(word) == m.WORD_PULSES
    assert m.s.interpreter.source_word(word, 9)['qubits'] == 9
    p = m.projector(vector)
    assert sum((p.get((i, i), m.ZERO)[0] for i in range(66)), F()) == 1
    assert p == m.adjoint(p)


def test_receipt_is_typed_and_bound_to_owning_sources(monkeypatch, tmp_path):
    assert all(name in packet.source_bindings() for name in packet.SOURCE_FILES)
    value = {'schema': 'control', 'scope': {'derived': True}}
    monkeypatch.setattr(packet, 'produce', lambda: deepcopy(value))
    path = tmp_path/'receipt.json'; path.write_text(packet.canonical(value))
    assert packet.verify(path)
    changed = deepcopy(value); changed['scope']['derived'] = 1
    path.write_text(packet.canonical(changed))
    with pytest.raises(ValueError): packet.verify(path)
    path.write_text('{"schema":"control","schema":"control"}')
    with pytest.raises(ValueError): packet.verify(path)
