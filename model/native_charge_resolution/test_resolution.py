from fractions import Fraction as F
import sympy as s

from . import bounds as b

def test_exact_roots_and_preserved_original_residue():
    value=b.calculate();rows=value['nested_window_frontier']
    same=next(r for r in rows if r['radius_in_units_of_full_hopping_norm']==128)
    assert F(same['unnormalized_bin_weight_interval'][0])>F(value['old_bin_lower_weight'])
    assert F(same['old_to_refined_width_ratio_lower'])>364
    sharp=rows[0]
    assert F(sharp['source_excitation_width'])<F('8.382e-8')
    assert F(sharp['unnormalized_bin_weight_interval'][0])>F('.31477')
    for row in rows:
        lo,hi=map(F,row['unnormalized_bin_weight_interval'])
        clo,chi=map(F,row['complete_outside_weight_interval'])
        assert lo+chi==hi+clo==F(1,2)
    assert not value['original_maximum_residue_eigenvalue_relocalized']

def test_full_neutral_shift_by_exact_positive_block_factorization():
    # For arbitrary scalar P/Q coordinates the complete-square form is
    # positive; matrix off-diagonal couplings are bounded by its norm.
    gap=s.Rational(7,3);v=s.Rational(1,19);c=gap-v;delta=v*v/c
    matrix=s.Matrix([[delta,-v],[-v,c+delta]])
    assert matrix[0,0]>0 and matrix.det()==delta**2>0
    # This verifies the sharp two-by-two lower eigenvalue is above-delta
    # using Sylvester, without evaluating any eigenvalue.
    assert s.Matrix([[0,-v],[-v,c]])+delta*s.eye(2)==matrix

def test_pythagorean_bound_and_complete_weight_with_dark_states():
    # A literal perturbed pair plus an unrelated dark state in the bin.
    # The selected atomic vector has a prescribed leaking norm.
    P=s.diag(1,0,1);psi=s.Matrix([s.Rational(12,13),s.Rational(5,13),0])
    assert (psi.T*psi)[0]==1
    leakage=(psi.T*(s.eye(3)-P)*psi)[0]
    assert (psi.T*P*psi)[0]==1-leakage==s.Rational(144,169)
    assert s.sqrt(1-leakage)>1-s.sqrt(leakage)
    for n in (2,4,8,128):
        r=b.root(F(1)-F(1,n*n))
        assert r*r<=F(1)-F(1,n*n)
        assert b.round_up(F(1,n+1))>=F(1,n+1)

def test_largest_residue_cannot_be_inferred_from_large_nested_weight():
    # Exact17-dimensional countermodel; all spectral information is known
    # by construction and checked without diagonalization.
    H=s.hadamard(16) if hasattr(s,'hadamard') else None
    if H is None:
        H=s.Matrix([[1]])
        for _ in range(4):H=H.row_join(H).col_join(H.row_join(-H))
    H=H/4;assert H*H.T==s.eye(16)
    U=s.diag(H,1);R=s.eye(17)
    R[0,0]=R[16,16]=s.sqrt(15)/4;R[0,16]=-s.Rational(1,4);R[16,0]=s.Rational(1,4)
    assert R*R.T==s.eye(17)
    energies=[s.Rational(2*j-15,150) for j in range(16)]+[s.Integer(3)]
    basis=R*U;actual=basis*s.diag(*energies)*basis.T
    assert actual*basis==basis*s.diag(*energies)
    weights=[s.simplify(basis[0,j]**2) for j in range(17)]
    assert weights[:16]==[s.Rational(15,256)]*16 and weights[16]==s.Rational(1,16)
    assert sum(weights)==1 and sum(weights[:16])==s.Rational(15,16)
    assert weights[-1]>max(weights[:-1])
    assert max(abs(e) for e in energies[:16])==s.Rational(1,10)
    assert 3*s.Rational(1,4)+s.Rational(1,10)==s.Rational(17,20)<1
    # The same construction at the actual128v bin scale retains almost
    # all weight and still puts its largest individual atom outside.
    count=2**19;outside=F(1,512**2)
    assert (1-outside)/count<outside and 1-outside>F('.999996')
    assert 256*F(1,512)+F(1,10)==F(3,5)<1
    assert 128*b.data()['v']<256*b.data()['v']<F(1,1024)

def test_public_receipt_preserves_refined_selection_distinction():
    import copy,pytest
    from .certificate import raw_certificate,compare_research,certificate,validate
    raw=raw_certificate();assert compare_research(raw)['all_scientific_fields_equal']
    wrong=copy.deepcopy(raw);wrong['original_maximum_residue_eigenvalue_relocalized']=True
    with pytest.raises(ValueError):compare_research(wrong)
    value=certificate();validate(value['input_sha256'])
    assert value['narrower_selection_is_explicit']
    assert not value['affine_or_Dirac_couplers_automatically_unchanged_by_refined_selection']
    assert 'model/native_charge_resolution/research_reference.json' in value['input_sha256']
    assert not value['construction_scope']['outer_affine_coupled_H_is_the_readout_source']
