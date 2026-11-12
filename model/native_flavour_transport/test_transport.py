import pytest
import sympy as s
from . import transport as t


def test_source_word_unitarity_and_exact_nonzero_CP():
    v,cube,j=t.invariant(t.PRIMARY_WORD)
    assert t.zero(v.H*v-s.eye(3))
    assert cube!=0 and j!=0
    assert t.source()['discriminant'].is_positive
    assert s.simplify(t.invariant(t.PRIMARY_WORD,s.I)[2]+j)==0


def test_all_predeclared_short_words_and_identity_retained():
    rows=t.controls()
    assert len(rows)==40
    assert [sum(len(r['word'])==n for r in rows)for n in range(4)]==[1,3,9,27]
    assert rows[0]['J_exact']=='0'
    assert any(r['J_exact']=='0'and len(r['word'])==2 for r in rows)
    assert any(r['J_exact']!='0'and len(r['word'])==1 for r in rows)


def test_complete_full_source_mass_current_and_original_records_without_calibration(monkeypatch):
    def forbidden():raise AssertionError('empirical parameter packet read')
    monkeypatch.setattr(t.matter,'parameters',forbidden)
    monkeypatch.setattr(t.matter,'ckm_matrix',forbidden)
    r=t.full_source_records()
    assert r['full45_mass_current_max_numeric_residual']<2e-12
    assert len(r['complete18_original_records'])==18
    assert abs(sum(x['probability']for x in r['complete18_original_records'])-1)<2e-12
    pair={x['Y']:x['probability']for x in r['complete18_original_records']if x['pair']==[0,1]}
    assert abs(pair['minus']-pair['plus']-r['J_numeric'])<2e-13
    assert all('exact_probability_formula'in row for row in r['complete18_original_records'])


def test_invalid_word_and_nonunitary_phase_fail():
    with pytest.raises(ValueError):t.word_unitary((9,))
    with pytest.raises(ValueError):t.word_unitary((2,),s.Rational(1,2))


def test_minimal_single_pulse_shortcut_on_actual_native_projector():
    r=t.minimal_pulse_certificate()
    assert r['generic_unit_circle_polynomial_remainder']=='0'
    assert r['all_three_weights_nonzero_exact']
    assert r['minimum_nonzero_pulse_count_in_predeclared_alphabet']==1
    assert r['single_pulse_J_numeric']!=0
