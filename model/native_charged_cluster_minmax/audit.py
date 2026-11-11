"""Independent exact minmax, corner symmetry and residual upper controls."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import json
import sympy as s
from model.field_record_addenda.native_orthogonal_wilson_charged_band import model as band
HERE=Path(__file__).resolve().parent
root=band.root;enc=band.enc

@lru_cache(None)
def compute():
    bridge,_,_,_=band.admit();a=band.native_algebra();H,Js,minimum,powers,_=a.cube()
    alpha,beta=a.current.matrices();Gamma=alpha[0]*alpha[1]*alpha[2]*beta
    pauli=a.pauli();spin=((s.eye(2)-s.I*sum(pauli,s.zeros(2)))/2).applyfunc(s.expand)
    assert (spin.H*spin-s.eye(2)).applyfunc(s.simplify)==s.zeros(2)
    assert (spin**3+s.eye(2)).applyfunc(s.simplify)==s.zeros(2)
    assert all((spin*pauli[j]*spin.H-pauli[(j+1)%3]).applyfunc(s.simplify)==s.zeros(2) for j in range(3))
    vertices=a.field.cube()['vertices'];perm=s.zeros(8);lookup={x:j for j,x in enumerate(vertices)}
    for j,(x,y,z) in enumerate(vertices):perm[lookup[(z,x,y)],j]=1
    R=s.kronecker_product(perm,s.diag(spin,spin))
    assert (R*H*R.H-H).applyfunc(s.expand)==s.zeros(32)
    assert (R**3+s.eye(32)).applyfunc(s.expand)==s.zeros(32)
    # Exact lowest positive free rank-two projector from the actual polynomial.
    Z=(H*H).applyfunc(s.expand);z0=s.Rational(5,4);I=s.eye(32)
    Q=((4*Z-13*I)*(16*Z*Z-136*Z+241*I)*(16*Z*Z-264*Z+977*I)).applyfunc(s.expand)
    denominator=(4*z0-13)*(16*z0*z0-136*z0+241)*(16*z0*z0-264*z0+977)
    EZ=(Q/denominator).applyfunc(s.expand);EA=EZ/2;EB=(EZ*H).applyfunc(s.expand)
    E=EA+EB/s.sqrt(5)
    assert (EZ*EZ-EZ).applyfunc(s.expand)==s.zeros(32) and s.trace(EZ)==4
    assert (EA*EA+EB*EB/5-EA).applyfunc(s.expand)==s.zeros(32)
    assert (EA*EB+EB*EA-EB).applyfunc(s.expand)==s.zeros(32) and s.trace(E)==2
    character=s.simplify(s.trace(E*R));character2=s.simplify(s.trace(E*R*R))
    # Reflection x<->y followed by ordinary T fixes the corner root.
    mirror=s.zeros(8)
    for j,(x,y,z) in enumerate(vertices):mirror[lookup[(y,x,z)],j]=1
    normal_spin=-s.I*(pauli[0]-pauli[1])/s.sqrt(2)
    Mspin=beta*s.diag(normal_spin,normal_spin)
    Tspin=s.diag(s.I*pauli[1],s.I*pauli[1]);anti=s.kronecker_product(mirror,Mspin*Tspin)
    assert (anti*H.conjugate()*anti.H-H).applyfunc(s.simplify)==s.zeros(32)
    assert (anti*anti.conjugate()-I).applyfunc(s.simplify)==s.zeros(32)
    assert (anti*R.conjugate()*anti.H-R.H).applyfunc(s.simplify)==s.zeros(32)
    T4=s.zeros(32);T6=s.zeros(32);geometry=a.field.cube()
    for i in range(12):
        for j in range(12):
            T4+=geometry['P4'][i,j]*Js[i]*Js[j]
            T6+=geometry['P6'][i,j]*Js[i]*Js[j]
    def rowsum(M):
        return max(sum(F(abs(s.re(x))+abs(s.im(x))) for x in M.row(i))for i in range(M.rows))
    r4,r6=rowsum(T4),rowsum(T6);cpbound=r4/4+r6/(2*root(6,False))
    g=F(1,64);v=F(1,16)+1/(12*root(6,False));q= root((13+5*root(6))/48)
    nu1=24*g*q+F(67,8)*g*g;nu0=24*g*q+8*g*g
    d1=(12*v+F(67,8))*g*g;d0=(12*v+8)*g*g
    l1=(root(5,False)/2,root(5)/2)
    l2=(root(17-4*root(3),False)/2,root(17-4*root(3,False))/2)
    deltaU=F(bridge['energy_alignment']['signed_code_vacuum_shift_over_m'][1])
    qcode=l2[0]-nu1-l1[1]
    gaps={str(c):min(l2[0]-nu1,c-deltaU)-l1[1]-d1 for c in (F(21,2),F(23,2))}
    assert all(x>0 for x in gaps.values())
    sigma0=F(bridge['fresh_vacuum']['new_pair_one_photon_current_square_over_m_squared'][1])
    sigma1=sigma0+cpbound
    # Odd sine remainder, even cosine/Coulomb part, and orthogonal photon parity.
    odd0=g*root(sigma0)+4*g**3*q**3;odd1=g*root(sigma1)+4*g**3*q**3
    even0=g*g*(12*root(3)*v+8);even1=g*g*(12*root(3)*v+F(67,8))
    r0sq=odd0*odd0+even0*even0;r1sq=odd1*odd1+even1*even1
    # PHP>=0 by bipartite-concavity of the free fixed-charge ground;
    # for a conservative bound we may instead allow PHP>=-d_q.
    lower0=-d0-r0sq/(2-nu0)
    lower1=l1[0]-d1-r1sq/qcode
    return {'source':'unchanged original Gamma bridge c21/2; c23/2 variant stated separately',
        'frozen_bridge_receipt':band.BRIDGE_SHA,
        'corner':{'actual_free_lowest_projector_rank':2,'C3_trace':str(character),'C3_square_trace':str(character2),
            'proper_rotation_lift_cubed':'-I','corner_mirror_times_T_square':'+I',
            'corner_mirror_times_T_conjugates_C3_to_inverse':True,
            'Kramers_or_exact_actual_doublet_degeneracy_claimed':False},
        'minmax':{'trial_V_upper_over_m':enc(d1),'neutral_trial_V_upper_over_m':enc(d0),
            'charged_code_gap_lower_over_m':enc(qcode-d1,False),
            'full_charged_rank2_gap_by_completion':{k:enc(v,False) for k,v in gaps.items()},
            'native_complement_included':True},
        'CP_current_bound':{'exact_T4_row_norm_upper':str(r4),'exact_T6_row_norm_upper':str(r6),
            'charged_A_squared_upper':enc(sigma1),'neutral_A_squared_upper':enc(sigma0),
            'charged_lowest_free_projector_executed':True},
        'conservative_actual_energy':{'neutral_code_shift_lower_over_m':enc(lower0,False),
            'neutral_code_shift_upper_over_m':enc(d0),'charged_reference_corrected_shift_lower_over_m':enc(lower1,False),
            'charged_reference_corrected_shift_upper_over_m':enc(l1[1]+d1),
            'charged_minus_actual_neutral_minus4_over_m':[enc(lower1-d0,False),enc(l1[1]+d1-lower0)],
            'charged_trial_residual_squared_upper':enc(r1sq),'neutral_trial_residual_squared_upper':enc(r0sq)},
        'diagnostic':{'original_full_rank2_gap':float(gaps['21/2']),'code_rank2_gap':float(qcode-d1),
            'neutral_interval':[float(lower0),float(d0)],'charged_transition_interval':[float(lower1-d0),float(l1[1]+d1-lower0)],
            'CP_current_bound':float(cpbound)}}

if __name__ == '__main__':
    print(json.dumps(compute(), indent=2, sort_keys=True))
