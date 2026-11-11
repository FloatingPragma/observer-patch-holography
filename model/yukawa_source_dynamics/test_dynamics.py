"""Independent source replay and adversarial compression checks."""
import numpy as np
import sympy as s
from model.common_source import native_repair_mass
from model.common_source.hierarchy_diagnosis import native_family_dynamics
from model.observer import native
from . import complete_c3, schur


def test_native_branch_replay_gives_full_response():
    x=schur.exact(); states=native.states(binary=True)
    transition=s.zeros(66)
    for before,flag,seam,coin,current,after,prob in native_family_dynamics.proposal_branches():
        transition[after,before]+=s.Rational(prob)
    loads=s.Matrix(12,66,lambda i,j:states[j][i])
    assert loads*transition==(s.eye(12)-x['H'])*loads
    assert all(sum(transition[:,i])==1 for i in range(66))
    assert schur.zero(x['D']*x['H']-x['H']*x['D'])
    assert not schur.zero(x['F']*x['A']-x['A']*x['F'])
    assert schur.zero(x['projected_commutator']-x['virtual_commutator'])


def test_virtual_elimination_cancels_two_apparent_masses():
    x=schur.exact()
    assert x['A'].rank()==3 and x['S'].rank()==1
    assert x['R0'].rank()==2
    assert schur.zero(x['S']*x['R0'])
    assert schur.zero(x['harmonic'].T*x['H']*x['harmonic']/x['metric']-x['S'])
    assert schur.zero(x['harmonic'].T*x['harmonic']/x['metric']-x['Z'])
    # Two protected zero modes cannot be removed by a family coordinate change.
    change=s.Matrix([[1,2,0],[0,1,1],[0,0,1]])
    assert (change.T*x['S']*change).rank()==1
    old=np.linalg.eigvals(np.linalg.solve(np.array(x['Z'],float),np.array(x['S'],float)))
    other=np.linalg.eigvals(np.linalg.solve(np.array(change.T*x['Z']*change,float),np.array(change.T*x['S']*change,float)))
    assert np.allclose(np.sort(old),np.sort(other),atol=2e-15,rtol=1e-12)


def test_full_resolvent_and_all_spectral_residues():
    x=schur.exact(); H=np.array(x['H'],float); B=np.array(x['B'],float)/np.sqrt(float(x['metric']))
    Q=np.linalg.qr(np.column_stack((B,np.eye(12))),mode='complete')[0][:,3:]
    for z in (-1,-0.5+0.2j):
        whole=B.T@np.linalg.solve(z*np.eye(12)-H,B)
        a=B.T@H@B; v=B.T@H@Q; c=Q.T@H@Q
        effective=np.linalg.inv(z*np.eye(3)-a-v@np.linalg.solve(z*np.eye(9)-c,v.T))
        assert np.allclose(whole,effective,atol=2e-14,rtol=2e-13)
    data=schur.spectral_data()
    assert len(data['poles'])==9
    assert all(row['residue_trace']>0 for row in data['poles'])
    assert data['residue_sum_error']<5e-13


def test_independent_native_fine_records_give_exact_cyclic_channel():
    # Rebuild routing directly from native raw labels, without the grouping
    # helper used by the producer.
    x=complete_c3.exact(); r=x['R']; out=s.zeros(9); probabilities={-1:s.S.Zero,0:s.S.Zero,1:s.S.Zero}
    for cc,cr,ec,er,advance,w in native_repair_mass.native_records():
        plus=cc==0 and cr==1; minus=cc==1 and cr==2; event=ec==0 and er in (4,5)
        exponent=1 if plus and not event else -1 if minus and not event else 0
        probabilities[exponent]+=w
        u=r**exponent
        out+=w*s.kronecker_product(u.conjugate(),u)
    assert probabilities=={-1:s.Rational(5,72),0:s.Rational(31,36),1:s.Rational(5,72)}
    assert complete_c3.clean(out-x['channel'])==s.zeros(9)
    assert x['channel'].rank()==9
    assert complete_c3.clean(x['channel']**3-x['twirl']-x['contraction']**3*(s.eye(9)-x['twirl']))==s.zeros(9)


def test_analyzer_complements_remove_selected_spectrum():
    x=complete_c3.exact()
    selected=x['analyzer_effects'][0]
    full=complete_c3.clean(sum(x['analyzer_effects'],s.zeros(3)))
    assert full==s.eye(3)
    assert selected!=full
    small=s.Rational(5,72)-5*s.sqrt(2)/108
    assert complete_c3.clean(selected-small*s.eye(3)).rank()==1
    assert len(full.eigenvals())==1


def test_scalar_record_phases_cannot_generate_family_phase():
    x=complete_c3.exact(); phase=(s.Integer(3)+4*s.I)/5
    out=s.zeros(9)
    for raw,exponent,weight,rotation in x['records']:
        u=phase**exponent*rotation
        out+=weight*s.kronecker_product(u.conjugate(),u)
    assert complete_c3.clean(out-x['channel'])==s.zeros(9)


def test_retained_quantum_records_have_no_oriented_family_phase():
    from . import relational
    x=relational.exact()
    assert x['record_vectors'].rows==128
    assert x['gram'].det()>0
    assert x['bargmann']==s.Rational(19,24)**3
    assert x['weak_records'].shape==(8,6)
    assert all(s.simplify(s.im(v))==0 for v in x['weak_triples'])
    # Source-fixed Fourier readout is a common unitary on flag space, so its
    # relational overlap invariants cannot create an additional family phase.
    e=x['record_vectors']
    swapped=e.copy(); swapped.row_swap(0,1)
    assert complete_c3.clean(swapped.H*swapped-x['gram'])==s.zeros(3)
