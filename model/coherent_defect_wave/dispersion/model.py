"""Exact native triplet energy branches of the declared positive local clock.

The invariant reference uses uniform link words at every cell. It differs
from the fixed-face vacuum in the general cubic record experiment.
"""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import sympy as s

from ..model import JointObserver, SparseState
from ..gauss import gamma_amplitudes


def laplacian(observer):
    out = s.zeros(12)
    for a,b in observer.law.edges:
        out[a,a] += 1; out[b,b] += 1
        out[a,b] -= 1; out[b,a] -= 1
    return out


@lru_cache(maxsize=2048)
def exact(value):
    if isinstance(value,F):
        return value
    value=s.simplify(value)
    return F(int(value.p),int(value.q)) if value.is_Rational else value


def port_embed(observer, vector, norm):
    """Uniform links and addresses; a selected triplet component and parity0."""
    if len(vector)!=12 or type(norm) is not int or norm<=0:
        raise ValueError('twelve amplitudes and positive integer squared norm required')
    values=tuple(exact(v) for v in vector)
    return SparseState({(g,p,0,0,a):values[p] for g in observer.link_words
                        for p in range(12) if vector[p]!=0 for a in range(60)},600*norm)


def positive_local_action(observer,state):
    """The actual A=I−(U+U*)/2, without a logarithm or imported mass matrix."""
    out=state.as_dict()
    for evolved in (observer.step(state),observer.inverse_step(state)):
        for key,value in evolved.amplitudes:
            out[key]=out.get(key,F(0))-value/2
    return SparseState({key:exact(v) for key,v in out.items()},state.norm_denominator)


def same(left,right):
    if left.norm_denominator!=right.norm_denominator:
        return False
    a,b=left.as_dict(),right.as_dict()
    return all(a.get(key,0)==b.get(key,0) or s.expand(a.get(key,0)-b.get(key,0))==0
               for key in set(a)|set(b))


def rational_parts(vector):
    z=s.sqrt(5)
    split=[s.expand(v) for v in vector]
    return ([s.expand(v).coeff(z,0) for v in split],
            [s.expand(v).coeff(z,1) for v in split])


@lru_cache(maxsize=1)
def source_certificate():
    observer=JointObserver(); L=laplacian(observer)
    roots=(s.Integer(0),5-s.sqrt(5),s.Integer(6),5+s.sqrt(5))
    ranks=(1,3,5,3)
    eigenrows=[]
    projectors=[]
    for lam,rank in zip(roots,ranks,strict=True):
        p=s.eye(12)
        for other in roots:
            if other!=lam:
                p=(L-other*s.eye(12))*p/(lam-other)
        p=p.applyfunc(s.simplify)
        assert p*p==p and p.T==p and s.trace(p)==rank
        assert (L*p-lam*p).applyfunc(s.simplify)==s.zeros(12)
        projectors.append(p)
        vector=12*p[:,0]
        norm=s.simplify((vector.T*vector)[0])
        assert norm==12*rank
        state=port_embed(observer,vector,int(norm)); energy=lam/60
        assert state.norm_squared()==1
        assert gamma_amplitudes(state.as_dict())==state.as_dict()
        # U has rational entries. Execute both rational coefficient vectors
        # of Q(sqrt5), avoiding repeated symbolic expansion of 7200 entries.
        parts=rational_parts(vector)
        expected_parts=rational_parts(energy*vector)
        for coefficients,expected in zip(parts,expected_parts,strict=True):
            source=port_embed(observer,coefficients,int(norm))
            target=port_embed(observer,expected,int(norm))
            assert same(positive_local_action(observer,source),target)
        eigenrows.append({'source_id':'blank-triplet-L-'+str(lam),
                          'carrier_laplacian_eigenvalue':str(lam),
                          'zero_momentum_source_energy':str(s.simplify(energy)),
                          'carrier_multiplicity_per_triplet_component':rank,
                          'localized_port_spectral_probability':str(s.Rational(rank,12)),
                          'normalized_port_vector':list(map(str,vector)),
                          'port_vector_squared_norm':str(norm),
                          'actual_joint_A_eigenidentity_verified':True})
    assert sum(projectors,s.zeros(12))==s.eye(12)
    assert sum(s.Rational(row['localized_port_spectral_probability']) for row in eigenrows)==1
    m0,m1,m2=roots[1]/60,roots[2]/60,roots[3]/60
    ratios=[s.Integer(1),s.simplify(m1/m0),s.simplify(m2/m0)]
    # The same hopping coefficient gives different low-momentum speeds for
    # different positive source rest energies. Overall units cannot remove this.
    speeds=[s.simplify(s.pi*m/2) for m in (m0,m1,m2)]
    assert len(set(speeds))==3
    m,D,q=s.symbols('m D q',positive=True)
    energy=m+D*q**2
    residual=s.expand(energy**2-m**2-2*m*D*q**2)
    assert residual==D**2*q**4
    # Reuse the owned reflection-walk theorem on the exactly intertwined
    # uniform-link triplet sector; these are register multiplicities.
    prior=json.loads((Path(__file__).resolve().parents[2]/'coherent_seam_repair/results.json').read_text())['walk']
    assert prior['full_walk_eigenvalue_plus_one_multiplicity']==38
    assert prior['full_walk_eigenvalue_minus_one_multiplicity']==1356
    full=[('0',38),(str(s.simplify(m0)),12),('1/15',10),('1/10',10),
          (str(s.simplify(m2)),12),('1/6',2),('2',1356)]
    assert sum(mult for _,mult in full)==1440
    return {'source_clock':'H=sum_x [I-Re(U_F,x)] +(pi/4)sum_edges dGamma(I-S_edge)',
            'reference':'empty-wave state with uniform10linkwords at every site',
            'separate_from_fixed_face_cubic_reference':True,
            'native_icosahedron_laplacian':[[int(x) for x in row] for row in L.tolist()],
            'blank_triplet_source_bands':eigenrows,
            'full_uniform_link_triplet_wave_spectrum_per_internal_component':[
                {'source_energy':energy,'register_multiplicity':mult} for energy,mult in full],
            'complete_localized_port_spectral_weight':'1',
            'spatial_band':'E_lambda(k)=lambda/60+(pi/2)*sum_i(1-cos(k_i))',
            'three_nonzero_source_energy_ratios_exact':list(map(str,ratios)),
            'three_nonzero_source_energy_ratios_approximate':[str(s.N(x,14)) for x in ratios],
            'leading_speed_squared_in_unit_spacing_and_time':list(map(str,speeds)),
            'same_speed_for_three_positive_branches':False,
            'continuum_quartic_mass_shell_residual':'D^2*|p|^4',
            'scope':{'exact_populated_three_dimensional_source_energy_bands':True,
                     'physical_particle_or_Standard_Model_identification':False,
                     'energies_are_physical_rest_masses':False,
                     'source_derived_GeV_masses':False,
                     'Lorentz_dispersion_on_this_family':False}}
