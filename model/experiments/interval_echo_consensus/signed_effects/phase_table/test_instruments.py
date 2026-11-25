"""Native complex effect tables, retained poststates and alternative read records."""

from unittest.mock import patch
import json

import numpy as np
import sympy as s
from scipy.linalg import block_diag

from . import acquire, measure as native


def test_exact_eight_slots_four_axes_and_positive_square_roots():
    projectors = [acquire.exact_projector(row["context"]) for row in native.context_table()]
    assert len(projectors)==8 and len({tuple(p) for p in projectors})==4
    rows = []
    for p in projectors:
        assert p*p==p and p.conjugate().T==p
        q = s.eye(2)-p
        e = s.eye(2)/4+p/2
        k0,k1 = s.sqrt(3)*p/2+q/2,p/2+s.sqrt(3)*q/2
        assert s.simplify(k0*k0-e)==s.zeros(2)
        assert s.simplify(k1*k1-(s.eye(2)-e))==s.zeros(2)
        assert e*e != e
        assert s.simplify((e-s.eye(2)/4)*(e-3*s.eye(2)/4))==s.zeros(2)
        rows.append([s.trace(p),p[0,0]-p[1,1],s.re(p[0,1]),s.im(p[0,1])])
    assert s.Matrix(rows).rank()==4
    assert projectors[-1]==s.Matrix([[s.Rational(1,2),-s.I/2],[s.I/2,s.Rational(1,2)]])


def test_source_primitive_formulas_in_full_99_state_space():
    law=native.candidate.native_law
    h0,h1=law.local_term_parts()
    h=h0+native.candidate.KAPPA*h1
    for sign in (-1,1):
        indices=[law.SRC.index(pair) for pair in ((1+sign,1),(1,1+sign))]
        j=np.eye(law.DIM)[:,indices]
        expected=np.array([[8.5,-0.5],[-0.5,8.5]])
        np.testing.assert_allclose(h@j,j@expected,atol=3e-13,rtol=0)
    reference=np.eye(law.DIM)[:,[law.SRC.index((load,0)) for load in range(3)]]
    np.testing.assert_allclose(native.candidate.KAPPA*h1@reference,
                               reference@np.diag([0,8,32]),atol=3e-13,rtol=0)


def test_explicit_source_words_realize_maps_on_both_signs():
    for row in native.context_table():
        context = row["context"]
        for mode in ("sharp","coarse","fine","coin"):
            expected = acquire.ideal_branches(context,mode)
            minus,plus = native.raw_branch_maps(context,mode,-1),native.raw_branch_maps(context,mode,1)
            for label in range(17):
                target = expected.get(label,np.zeros((2,2)))
                for actual in (minus[label],plus[label]):
                    np.testing.assert_allclose(actual[:2],target,atol=5e-13,rtol=0)
                    np.testing.assert_allclose(actual[2:],0,atol=5e-13,rtol=0)
                np.testing.assert_allclose(minus[label],plus[label],atol=5e-13,rtol=0)
            np.testing.assert_allclose(sum(k.conj().T@k for k in minus.values()),np.eye(2),atol=5e-13,rtol=0)


def test_unread_logical_record_is_reversibly_removed_and_tags_untouched():
    for mode,kept,cleared in (("coarse",14,[0]),("fine",20,[])):
        result = acquire.cleared_record_check(mode)
        assert result["physical_code_inputs"]==384
        assert result["cleared_bits"]==cleared and result["retained_native_copy_calls"]==kept
        assert result["phase_comparison_error"]<5e-13
        endpoints = {p for _,a,b,_ in native.record_word(mode) for p in (a,b)}
        assert endpoints.isdisjoint({76,77,79,239})
    for row in native.context_table():
        for mode in ("sharp","coarse","fine","coin"):
            word = native.context_word(row["context"])+native.dilation_word(mode)
            for raw in range(4):
                word += native.post_word(row["context"],mode,raw)
            assert {p for _,a,b,_ in word for p in (a,b)} <= set(range(12))


def test_equal_effects_have_distinct_conditional_maps_and_union_rules():
    for row in native.context_table():
        c = row["context"]
        fine = native.fine_matrices(c)
        coarse = native.branch_matrices(c,"weak")
        np.testing.assert_allclose(sum(fine[r].conj().T@fine[r] for r in (0,1)),
                                   coarse[0].conj().T@coarse[0],atol=5e-13,rtol=0)
        np.testing.assert_allclose(sum(fine[r].conj().T@fine[r] for r in (0,3)),
                                   3*np.eye(2)/4,atol=5e-13,rtol=0)
    psi = np.array([1,1])/np.sqrt(2)
    rho = np.outer(psi,psi)
    c = native.branch_matrices("phase","weak")
    f = native.fine_matrices("phase")
    coarse0 = c[0]@rho@c[0].conj().T
    refined0 = sum(f[r]@rho@f[r].conj().T for r in (0,1))
    assert np.linalg.norm(coarse0-refined0,2)>0.2
    coin = native.branch_matrices("phase","coin")
    np.testing.assert_allclose(coin[0]@rho@coin[0].conj().T,3*rho/4,atol=5e-13,rtol=0)
    np.testing.assert_allclose(sum(f[r]@rho@f[r].conj().T for r in (0,3)),
                               3*np.eye(2)/8,atol=5e-13,rtol=0)


def test_reference_complete_maps_include_cross_sign_coherence():
    for context in ("web_conjugated_2","phase"):
        for mode in ("sharp","coarse","fine","coin"):
            minus,plus = native.raw_branch_maps(context,mode,-1),native.raw_branch_maps(context,mode,1)
            expected = acquire.ideal_branches(context,mode)
            # Four-dimensional Choi state includes all sign/rail matrix units.
            bell = np.eye(4).reshape(-1)/2
            marginal = np.zeros((4,4),complex)
            for label in range(17):
                k = block_diag(minus[label],plus[label])
                out = np.kron(k,np.eye(4))@bell
                reshaped = out.reshape(24,4)
                marginal += reshaped.T@reshaped.conj()
                if label in expected:
                    target = np.zeros((24,4),complex)
                    target[:2,:2]=target[12:14,2:]=expected[label]
                    target_out = np.kron(target,np.eye(4))@bell
                    np.testing.assert_allclose(np.outer(out,out.conj()),np.outer(target_out,target_out.conj()),
                                               atol=5e-13,rtol=0)
            np.testing.assert_allclose(marginal,np.eye(4)/4,atol=5e-13,rtol=0)


def test_native_phase_orientation_control_and_poststate_readback():
    original = native.context_word
    def changed(context):
        return [native.pulse("contact",0,1,native.Fraction(1,2))] if context=="phase" else original(context)
    psi = acquire.input_vector("plus_i",-1)
    ordinary = native.branch_matrices("phase","weak")
    with patch.object(native,"context_word",changed):
        rival = native.branch_matrices("phase","weak")
    assert abs(np.linalg.norm(ordinary[0]@psi)**2-0.75)<5e-13
    assert abs(np.linalg.norm(rival[0]@psi)**2-0.25)<5e-13
    for sign in (-1,1):
        for mode,yes in (("coarse",(2+np.sqrt(3))/4),("fine",0.5),("coin",1)):
            law = acquire.readback_probabilities("phase",mode,"plus",sign)
            assert abs(law.sum()-1)<5e-13 and abs(law[:,0].sum()-yes)<5e-13
            assert law[:,2:].sum()<1e-24


def test_complete_acquisition_banks_and_conditional_history_counts():
    saved = json.loads(acquire.RECEIPT.read_text())
    native.sequential.compare(saved,acquire.report())
    assert len(saved["banks"])==250
    assert sum(bank["shots"] for bank in saved["banks"])==1024000
    for bank in saved["banks"]:
        assert len(bank["all_17_counts"])==17
        assert sum(bank["all_17_counts"])==bank["shots"]
        assert sum(bank["all_17_counts"][4:])==0
    assert len(saved["readback_banks"])==6
    for bank in saved["readback_banks"]:
        bins=np.array(bank["all_17_by_17_counts"])
        assert bins.shape==(17,17) and bins.sum()==bank["shots"]
        measured=bins[:,0].sum()/bank["shots"]
        if bank["mode"]=="coarse":
            assert measured>0.9
        elif bank["mode"]=="fine":
            assert 0.48<measured<0.52
        else:
            assert measured==1
