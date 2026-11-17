"""Actual 12-CAR one-cell spin/flavour branch and mixed colour-spin current.

This is a local candidate interaction, not a joined baryon mass solution.
The auxiliary-source eigenvalues are derived from a literal current coupling,
with inherited t and adjoint electric gap; no hyperfine coefficient is fitted.
"""
from pathlib import Path
from fractions import Fraction as F
from math import lcm,sqrt
from itertools import product
import json
import numpy as np
import sympy as sy
from staggered_chain4 import sector,colour6
from staggered_vacuum import creation,annihilation
HERE=Path(__file__).resolve().parent

def current(words,pairs):
    out={}
    for j,w in enumerate(words):
        for a,b,c in pairs:
            z,s=annihilation(b,w)
            if z is None:continue
            y,t=creation(a,z)
            if y is not None:
                if y not in out:out[y]=np.zeros(len(words),dtype=np.int64)
                out[y][j]+=c*s*t
    return np.array(list(out.values()),dtype=np.int64) if out else np.zeros((0,len(words)),dtype=np.int64)

def inside(words,pairs):
    ix={w:i for i,w in enumerate(words)};A=np.zeros((len(words),)*2,dtype=np.int64)
    for j,w in enumerate(words):
        for a,b,c in pairs:
            z,s=annihilation(b,w)
            if z is None:continue
            y,t=creation(a,z)
            if y is not None:A[ix[y],j]+=c*s*t
    return A

def main():
    s=sector(4,1);words=s['words'];R=s['R'];n=len(words)
    adjoint_zero_weight_multiplicity=len(sy.Matrix((colour6(words,range(4))-18*np.eye(n,dtype=np.int64)).tolist()).nullspace())
    assert adjoint_zero_weight_multiplicity==40
    # alpha=(spin-up flavour-up, spin-up flavour-down,
    #        spin-down flavour-up, spin-down flavour-down).
    Sp=inside(words,[(3*f+c,3*(f+2)+c,1) for f in range(2) for c in range(3)])
    Ip=inside(words,[(3*a+c,3*(a+1)+c,1) for a in [0,2] for c in range(3)])
    Zs=Sp@Sp.T-Sp.T@Sp;Zi=Ip@Ip.T-Ip.T@Ip
    assert np.array_equal(Zs@Sp-Sp@Zs,2*Sp)
    assert np.array_equal(Zi@Ip-Ip@Zi,2*Ip)
    assert not np.any(Sp@Ip-Ip@Sp)
    S4=Zs@Zs+2*(Sp.T@Sp+Sp@Sp.T)
    G24=np.zeros((n,n),dtype=np.int64)
    for a,b in product(range(3),repeat=2):
        Z=current(words,[(3*alpha+a,3*alpha+b,1 if alpha<2 else -1) for alpha in range(4)])
        P=current(words,[(3*f+a,3*(f+2)+b,1) for f in range(2)])
        M=current(words,[(3*(f+2)+a,3*f+b,1) for f in range(2)])
        G24+=3*Z.T@Z+6*(P.T@P+M.T@M)
    G24-=S4
    assert np.array_equal(G24@Sp,Sp@G24) and np.array_equal(G24@Ip,Ip@G24)
    # The complete colour-singlet triquark operator identity, not a fitted S².
    assert np.array_equal(G24@R,108*R-4*S4@R)
    branches=[]
    for S,I in product([1,3],repeat=2):
        con=np.vstack([(Zs-S*np.eye(n,dtype=np.int64))@R,(Zi-I*np.eye(n,dtype=np.int64))@R,Sp@R,Ip@R])
        null=sy.Matrix(con.tolist()).nullspace()
        if not null:continue
        cols=[]
        for v in null:
            den=lcm(*(int(x.q) for x in v));cols.append([int(x*den) for x in v])
        V=R@np.array(cols,dtype=np.int64).T
        eigen=F(9,2)-F(2,3)*F(S*(S+2),4)
        assert np.array_equal(G24@V,int(24*eigen)*V)
        branches.append(dict(twice_spin=S,twice_isospin=I,highest_weight_multiplicity=V.shape[1],
             complete_dimension=(S+1)*(I+1)*V.shape[1],mixed_current_square_eigenvalue=str(eigen)))
    assert sum(r['complete_dimension'] for r in branches)==20
    inp=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text());unit=F(inp['electron_anchor_GeV'])/F(inp['432']['declared_light_midpoint'])
    shifts=[]
    for kappa in [F(1,2),F(1),F(2)]:
        gap=3*kappa;g=F(1,8)
        energies={str(r['twice_spin']):float((float(gap)-sqrt(float(gap*gap+4*g*g*F(r['mixed_current_square_eigenvalue']))))/2) for r in branches}
        split=energies['3']-energies['1']
        assert split>0
        shifts.append(dict(kappa=str(kappa),auxiliary_adjoint_electric_gap=str(gap),current_coupling=str(g),
            exact_formula='E_minus(lambda)=(Delta-sqrt(Delta^2+4*g^2*lambda))/2',
            lower_energies_source_diagnostic=energies,delta_minus_nucleon_source_diagnostic=split,
            conditional_delta_minus_nucleon_GeV_diagnostic=split*float(unit)))
    out=dict(schema='oph.obstruction_repair.local_spin_flavour_current.v1',
        full_quark_CAR_dimension=4096,complete_three_quark_CAR_dimension=220,colour_singlet_triquark_dimension=20,
        actual_spin_and_isospin_currents=True,exact_adjoint_zero_weight_dimension=adjoint_zero_weight_multiplicity,
        exact_colour_adjoint_multiplicity=adjoint_zero_weight_multiplicity//2,all_highest_weight_branches=branches,
        exact_mixed_current_identity='sum_ai G_ai^2 = 9/2 - (2/3) S^2 on the complete colour-singlet 3-quark sector',
        current='G_ai=sum_flavour,spin,colour c† (T_a tensor sigma_i/2) c',
        generator_normalization='tr(T_a T_b)=delta_ab/2; sigma_i/2 spin',
        proposed_auxiliary_source='H_aux=3*kappa*P_adjoint_vector + t sum_ai (|ai><0|+|0><ai|)G_ai',
        auxiliary_carrier='one scalar vacuum plus an SU3 adjoint times physical-spin vector, dimension25',
        literal_auxiliary_gap_ancestry='kappa times SU3 adjoint quadratic Casimir3',
        literal_current_coupling_ancestry='inherited source t=1/8, declared new architecture; not forced by OPH',
        zero_quark_sector_auxiliary_vacuum_unchanged=True,
        source_wide_ground_vacuum_preserved=False,
        caveat='negative local baryon shift can alter the full vacuum; actual joint two-cell source must be re-grounded',
        all_physical_one_cell_N3_plus_auxiliary_colour_singlet_dimension=80,
        accessible_coupled_doublets=20,retained_dark_auxiliary_states=40,
        shifts=shifts,empirical_hyperfine_fit=False,
        auxiliary_candidate_joined_to_two_cell_vacuum=False,
        common_Dirac_or_magnetic_field_transfer=False,physical_mass_prediction=False,
        native_preparation_and_all_records_executed=False)
    (HERE/'LOCAL_SPIN_FLAVOUR_CURRENT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
