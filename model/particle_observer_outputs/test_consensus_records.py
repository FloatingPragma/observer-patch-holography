import sympy as s
import pytest

from . import consensus_records as c


@pytest.mark.parametrize('setting', ['position', 'momentum'])
def test_complete_pointer_copy_and_record_isometry(setting):
    pointer = c.pointer_unitary(setting)
    assert c._clean(pointer.H*pointer) == s.eye(16)
    for encoding in (c.ALICE_ENCODING, c.BOB_ENCODING):
        copy = c.copy_unitary(encoding)
        assert copy.H*copy == s.eye(20)
        for j in range(4):
            blank = s.eye(20)[:, 5*j]
            assert copy*blank == s.eye(20)[:, 5*j+encoding[j]]
    v = c.record_isometry(setting)
    assert v.shape == (400, 4)
    assert c._clean(v.H*v) == s.eye(4)


def test_copying_an_orthogonal_record_does_not_clone_a_superposition():
    pointer = (s.eye(4)[:, 0]+s.eye(4)[:, 1])/s.sqrt(2)
    blank = s.eye(5)[:, 0]
    output = c.copy_unitary(c.ALICE_ENCODING)*s.kronecker_product(pointer, blank)
    coefficient = s.Matrix(4, 5, list(output))
    reduced = c._clean(coefficient*coefficient.H)
    assert reduced == s.diag(s.Rational(1, 2), s.Rational(1, 2), 0, 0)
    assert reduced != pointer*pointer.H


@pytest.mark.parametrize('setting,expected', [
    ('position', (s.Rational(1, 2), s.Rational(1, 2), 0, 0)),
    ('momentum', (s.Rational(1, 2), s.Rational(1, 4), 0, s.Rational(1, 4))),
])
def test_complete_joint_record_law_and_both_observer_marginals(setting, expected):
    rows = c.joint_law(setting)
    assert len(rows) == 100
    assert c.marginal(rows, 'raw_event') == {(j,): p for j, p in enumerate(expected)}
    assert c.marginal(rows, 'comparison')['agreement',] == 1
    for observer in ('alice_meaning', 'bob_meaning'):
        marginal = c.marginal(rows, observer)
        assert all(marginal[j,] == p for j, p in enumerate(expected))
        assert marginal[None,] == 0
    projectors = tuple(c.menu(setting)['projectors'].values())
    repeat = c.measurement.sequential_records(c.menu(setting)['state'], projectors, projectors)
    assert len(repeat) == 16
    assert all(s.sympify(row['probability']) == (expected[row['first']] if row['first'] == row['second'] else 0)
               for row in repeat)


def test_mistaken_dictionary_is_exposed_and_repair_keeps_original_histogram():
    wrong = c.dictionary(c.BOB_ENCODING)
    wrong[c.BOB_ENCODING[0]], wrong[c.BOB_ENCODING[1]] = 1, 0
    incorrect = c.joint_law('momentum', bob_dictionary=wrong)
    corrected = c.joint_law('momentum')
    assert c.marginal(incorrect, 'comparison')['disagreement',] == s.Rational(3, 4)
    assert c.marginal(corrected, 'comparison')['agreement',] == 1
    assert c.marginal(incorrect, 'raw_event') == c.marginal(corrected, 'raw_event')
    assert [(r['alice_record'], r['bob_record'], r['probability']) for r in incorrect] == [
        (r['alice_record'], r['bob_record'], r['probability']) for r in corrected]


def test_complete_failure_channel_keeps_reference_and_no_postselection():
    a = c.read_channel(c.ALICE_ENCODING, s.Rational(1, 20))
    b = c.read_channel(c.BOB_ENCODING, s.Rational(1, 10))
    maps = c.joint_kraus('momentum', a, b)
    assert len(maps) == 100
    assert c._clean(sum((k.H*k for k in maps.values()), s.zeros(4))) == s.eye(4)
    rows = c.joint_law('momentum', alice_channel=a, bob_channel=b)
    comparison = c.marginal(rows, 'comparison')
    assert comparison['failure',] == s.Rational(29, 200)
    assert comparison['agreement',] == s.Rational(171, 200)
    # Untouched two-dimensional reference, including its conditional records.
    bell = (s.eye(8)[:, 0]+s.eye(8)[:, 3])/s.sqrt(2)
    rho = bell*bell.H
    outputs = [c._clean(s.kronecker_product(k, s.eye(2))*rho*s.kronecker_product(k.H, s.eye(2)))
               for k in maps.values()]
    total = c._clean(sum(outputs, s.zeros(8)))
    reference = s.Matrix(2, 2, lambda i, j: sum(total[2*x+i, 2*x+j] for x in range(4)))
    assert reference == s.eye(2)/2
    assert sum(s.trace(output) for output in outputs) == 1
    assert sum(s.trace(output) for record, output in zip(maps, outputs, strict=True)
               if c.FAILURE in record[1:]) == s.Rational(29, 200)


def test_noninjective_meanings_lose_information_and_cannot_create_a_histogram():
    meanings = ('even', 'odd', 'even', 'odd')
    assert c.pushforward((1, 0, 0, 0), meanings) == c.pushforward((0, 0, 1, 0), meanings)
    result = c.coarse_error_bound((s.Rational(1, 2), s.Rational(1, 4), 0, s.Rational(1, 4)),
                                  (0, s.Rational(1, 4), s.Rational(1, 2), s.Rational(1, 4)), meanings)
    assert result['complete_event_TV'] == s.Rational(1, 2)
    assert result['coarse_record_TV'] == 0
    with pytest.raises(ValueError, match='injective'):
        c.copy_unitary((0, 0, 2, 3))


def test_actual_cover_glues_functions_and_rejects_disagreement_or_missing_events():
    local = [{0: 'x0', 1: 'x1', 2: 'x2'}, {1: 'x1', 2: 'x2', 3: 'x3'}]
    assert c.glue_sections(range(4), local) == {j: 'x'+str(j) for j in range(4)}
    with pytest.raises(ValueError, match='disagree'):
        c.glue_sections(range(4), [local[0], {1: 'wrong', 2: 'x2', 3: 'x3'}])
    with pytest.raises(ValueError, match='cover'):
        c.glue_sections(range(4), [local[0]])


def test_incomplete_channels_are_rejected_without_renormalization():
    broken = s.MutableMatrix(c.read_channel(c.ALICE_ENCODING))
    broken[0, 0] = s.Rational(1, 2)
    with pytest.raises(ValueError, match='original outcome weight'):
        c.joint_kraus('position', broken)
    with pytest.raises(ValueError, match='positive'):
        c.joint_law('position', s.diag(2, -1, 0, 0))


def test_certificate_exposes_actual_records_and_description_boundary():
    receipt = c.certificate()
    assert receipt['wrong_dictionary_comparison']['disagreement'] == '3/4'
    assert receipt['retained_failure_comparison']['failure'] == '29/200'
    assert not receipt['physical_interpretation']['record_agreement_proves_three_plus_one_continuum']
    assert len(receipt['graphic_nodes']) == 5
