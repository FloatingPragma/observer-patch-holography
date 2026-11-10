"""Full algebra and same-H spectral bridge regression controls."""
from itertools import product
import sympy as s
from .model import (normal_order,monomial,multiply,adjoint,commutator,anticommutator,
    column,quantum_yukawa,screened_eom,composite,neutral_composite,gauge_current,
    tensor_certificate,compiler_certificate,noncanonical_certificate,
    composite_dynamics_certificate,add,scale)
from ..model import sequence


def test_car_normal_ordering_against_complete_three_mode_car():
    # Every word of length <=3 on six elementary operators is compared with
    # all eight occupation states. Adjoint products exercise length-six words.
    atoms=tuple(product(range(3),(False,True)))
    words=[()]+[(a,) for a in atoms]+list(product(atoms,repeat=2))+list(product(atoms,repeat=3))
    for word in words:
        for mask in range(8):
            raw=sequence(mask,word);out={}
            for normal,coefficient in normal_order(word):
                step=sequence(mask,normal)
                if step is not None:out[step[0]]=out.get(step[0],0)+coefficient*step[1]
            out={k:v for k,v in out.items() if v}
            assert out==({raw[0]:raw[1]} if raw is not None else {})
        op=monomial(word)
        assert adjoint(adjoint(op))==op
        assert adjoint(multiply(op,adjoint(op)))==multiply(op,adjoint(op))


def test_full_current_spin_eom_certificate():
    proof=tensor_certificate()
    assert proof['full_CAR_creation_EOM_identities']==30
    assert proof['screened_current_component_identities']==192
    assert proof['spin_component_identities']==48


def test_exact_compiler_ancestry():
    assert compiler_certificate()['exact_equality_to_periodic_full_family_Yukawa_word_compiler']


def test_neutral_and_charged_are_distinct_gauge_characters():
    for spin in (0,1):
        assert not commutator(gauge_current(11),neutral_composite(spin))
        assert commutator(gauge_current(11),composite('charged_lepton',9,spin))==composite('charged_lepton',9,spin)


def test_all_composites_are_noncanonical():
    rows=noncanonical_certificate()
    assert {row['channel'] for row in rows}=={'up','down','charged_lepton','neutral'}
    assert all(not row['canonical_CAR'] for row in rows)


def test_quantum_scalar_terms_obstruct_linear_dirac_closure():
    rows=composite_dynamics_certificate()['exact_quantum_composite_EOM']
    assert all(row['nonzero_cubic_CAR_terms']>0 for row in rows)


def test_wrong_eom_sign_is_rejected():
    parts,ys=quantum_yukawa()
    actual=commutator(parts['charged_lepton'],monomial(((9,True),)))
    correct=scale(screened_eom('charged_lepton',9,0),ys['charged_lepton'])
    assert actual==correct and actual!=scale(correct,-1)


def test_hopping_eom_identity_retains_boson_matrix_unit():
    matrix=s.zeros(5);matrix[1,0]=1
    term=monomial(((2,True),(29,False)),matrix)
    assert commutator(term,monomial(((29,True),)))==monomial(((2,True),),matrix)
    assert commutator(term,monomial(((2,False),)))==scale(monomial(((29,False),),matrix),-1)
    assert not commutator(term,monomial(((5,True),)))


def test_periodic_same_h_transfer_matrices():
    from .periodic_bridge import transfer_matrices
    for kind in ('add','remove'):
        target,t,j,k=transfer_matrices(kind)
        assert t.shape==j.shape==k.shape
        assert j and k


def test_local_neutral_composite_obeys_exact_star_selection():
    from .model import star_selection_certificate
    proof=star_selection_certificate()
    assert proof['local_composite_star_character_checks']==16
    assert not proof['neutral_current_singlet_evades_star_selection']
    assert not proof['retained_link_string_dressed_correlators_excluded']


def test_positive_same_H_first_composite_overlap():
    from .periodic_bridge import first_residue_enclosure
    for kind in ('add','remove'):
        proof=first_residue_enclosure(kind)
        lo,hi=map(s.Rational,proof['quantum_composite_residue_interval'])
        assert s.Rational(4230841,10000000)<lo<hi<s.Rational(4230842,10000000)


def test_transitive_custody_and_consumed_receipt():
    from .certificate import input_hashes,validate_pins,CONSUMED_RECEIPT
    pins=input_hashes()
    prefix='model/coherent_defect_wave/chiral_joint/higgs_composites/'
    assert all(prefix+name in pins for name in ('model.py','certificate.py','__main__.py','test_model.py','PROOF.md','README.md'))
    assert CONSUMED_RECEIPT in pins
    assert 'model/coherent_defect_wave/chiral_joint/field_limit/PROOF.md' in pins
    validate_pins(pins)


def test_replay_tolerance_never_weakens_certified_interval_or_custody():
    from .certificate import replay_differences
    assert not replay_differences({'target':[[0,-1],1]},{'target':[(0,-1),1]})
    assert not replay_differences({'diagnostic':.4},{'diagnostic':.4+1e-12})
    assert replay_differences({'interval':['1/3','1/2']},{'interval':['1/4','1/2']})
    assert replay_differences({'input_sha256':{'x':'abc'}},{'input_sha256':{'x':'def'}})
    assert replay_differences({'diagnostic':.4},{'diagnostic':.4+1e-7})
