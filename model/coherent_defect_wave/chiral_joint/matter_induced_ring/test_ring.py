"""Complete-domain controls and hostile projection/resolvent checks."""
import pytest
import sympy as s
from .model import (physical_split,norm_certificate,split_certificate,
                    ancestry_certificate,resolvent_certificate,eliminated_kernel)
from .records import complete_records


def test_complete_gauge_projector_and_actual_full_H_split():
    row=split_certificate()
    assert row['retained_dimension']+row['eliminated_dimension']==72
    assert row['full_CAR_projector_current_commutators']==24
    assert row['complete_embedding_endpoint_Gauss_identities']==24
    assert norm_certificate()['full_Hamiltonian_norm_upper_bound']=='30'
    block=physical_split()
    # Omitting the scalar/link P states changes the declared projection.
    assert len(block['p'])==5 and len(block['ep'].todok())==14
    assert block['b']!=s.zeros(67,5)


def test_irreducible_ring_uses_actual_edges_and_quantum_Yukawa():
    row=ancestry_certificate()
    histories=row['nonzero_minimal_term_histories']
    assert sum(s.Rational(h['coefficient']) for h in histories if h['target']=='singlet_loop')==s.Rational(1,16)
    assert sum(s.Rational(h['coefficient']) for h in histories if h['target']=='doublet_loop')==s.Rational(3,8)
    assert row['complete_Q_irreducible_moments'][0]['blank_return_self_energy']!='0'


def test_exact_energy_dependent_kernel_and_reconstruction():
    row=resolvent_certificate()
    assert all(s.Rational(v['blank_to_singlet_loop'])<0 and s.Rational(v['blank_to_doublet_loop'])<0
               for v in row['evaluations'])
    assert row['energy_dependent_loop_coefficients_verified']
    for z in (-30,-29,0,1,s.Rational(-61,2)+s.sqrt(2)):
        with pytest.raises(ValueError):eliminated_kernel(z)


def test_complete_positive_time_records_use_full_H_residual():
    row=complete_records()
    assert [r['effect_rank'] for r in row['outcomes']]==[1,1,70]
    assert all(0<s.Rational(r['probability_interval'][0])<=s.Rational(r['probability_interval'][1])<1
               for r in row['outcomes'])
    assert sum(s.Rational(r['probability_interval'][0]) for r in row['outcomes'])<=1
    assert sum(s.Rational(r['probability_interval'][1]) for r in row['outcomes'])>=1
    assert row['failure_and_all_other_states_retained']
