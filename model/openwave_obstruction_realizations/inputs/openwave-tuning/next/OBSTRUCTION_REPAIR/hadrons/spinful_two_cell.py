"""New two-spatial-cell, two-physical-spin SU3 source, all 4096 CAR states.

The carrier size matches the four-spinless-site experiment, but the
Hamiltonian, local gauge cut and SU2 currents below are newly constructed.
No source equivalence is asserted.
"""
from fractions import Fraction as F
from pathlib import Path
from math import lcm,floor
import json
import numpy as np
import scipy.linalg as la
import sympy as sy
from staggered_chain4 import sector,colour6,inertia_at
from staggered_vacuum import creation,annihilation
HERE=Path(__file__).resolve().parent


def action(words,orbital_pairs):
    ix={w:i for i,w in enumerate(words)};out=np.zeros((len(words),len(words)),dtype=np.int64)
    for j,w in enumerate(words):
        for a,b in orbital_pairs:
            for c in range(3):
                z,s=annihilation(3*b+c,w)
                if z is None:continue
                y,t=creation(3*a+c,z)
                if y is not None:out[ix[y],j]+=s*t
    return out


def frames(s):
    words=s['words'];R=s['R'];n=len(words);N=s['N']
    # orbital0=(left,up),1=(left,down),2=(right,up),3=(right,down).
    M=np.diag([sum((1 if x<2 else -1)*((w>>(3*x))&7).bit_count() for x in range(4)) for w in words])
    K=-action(words,[(0,2),(2,0),(1,3),(3,1)])
    E=colour6(words,[0,1])
    plus=action(words,[(0,1),(2,3)])
    minus=plus.T
    z2=np.array([sum((1 if x%2==0 else -1)*((w>>(3*x))&7).bit_count() for x in range(4)) for w in words])
    Z=np.diag(z2)
    assert np.array_equal(Z@plus-plus@Z,2*plus)
    assert np.array_equal(plus@minus-minus@plus,Z)
    for X in [M,K,E]:
        assert np.array_equal(X@plus,plus@X)
        assert np.array_equal(X@Z,Z@X)
    result=[]
    for twos in range(N%2,min(N,12-N)+1,2):
        unwanted=np.flatnonzero(z2!=twos)
        constraints=np.vstack([R[unwanted],plus@R])
        null=sy.Matrix(constraints.tolist()).nullspace()
        if not null:continue
        columns=[]
        for v in null:
            den=lcm(*(int(x.q) for x in v));columns.append([int(x*den) for x in v])
        Q=R@np.array(columns,dtype=np.int64).T
        assert not np.any(plus@Q)
        assert np.array_equal(Z@Q,twos*Q)
        G=Q.T@Q
        result.append(dict(twos=twos,R=Q,G=G,components=[Q.T@X@Q for X in [M,K,E]],rank=Q.shape[1]))
    assert sum((r['twos']+1)*r['rank'] for r in result)==s['singlet_dimension']
    return result


def bound_ground(A,G,values):
    mult=int(sum(abs(v-values[0])<1e-7 for v in values))
    den=1<<30;center=F(floor(float(values[0])*den),den);radius=F(2,den)
    for _ in range(12):
        lo,hi=center-radius,center+radius
        if inertia_at(A,G,lo)==0 and inertia_at(A,G,hi)==mult:return lo,hi,mult
        radius*=2
    raise ArithmeticError('ground band not certified')


def main():
    inp=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text())
    original_m=F(inp['432']['declared_light_midpoint']);m=F(round(original_m*(1<<36)),1<<36)
    mass_error=12*abs(m-original_m)
    sectors=[]
    for k in range(5):
        print('building spinful N',3*k,flush=True)
        s=sector(4,k);s['spin_frames']=frames(s);sectors.append(s)
        print('highest weight multiplicities',[(r['twos'],r['rank']) for r in s['spin_frames']],flush=True)
    results=[]
    for kappa in [F(1,2),F(1),F(2)]:
        records=[];numerical={}
        for s in sectors:
            for r in s['spin_frames']:
                M,K,E=r['components'];G=r['G'];n=r['rank']
                A=[[m*int(M[i,j])+F(1,8)*int(K[i,j])+kappa*F(int(E[i,j]),6) for j in range(n)] for i in range(n)]
                values,vectors=la.eigh(np.array(A,dtype=float),G.astype(float))
                numerical[s['N'],r['twos']]=(values,vectors,r)
                if s['N'] in [0,12]:lo=hi=F();mult=1
                else:lo,hi,mult=bound_ground(A,G,values);lo-=mass_error;hi+=mass_error
                records.append(dict(N=s['N'],twice_physical_spin=r['twos'],highest_weight_dimension=n,
                    complete_spin_sector_dimension=(r['twos']+1)*n,
                    lowest_band_exact=list(map(str,[lo,hi])),lowest_highest_weight_multiplicity=mult,
                    all_highest_weight_energies_diagnostic=values.tolist()))
                print('k',kappa,'N',s['N'],'2S',r['twos'],'E',float(lo),float(hi),flush=True)
        lowest=min(records,key=lambda r:F(r['lowest_band_exact'][0]))
        vlo,vhi=map(F,lowest['lowest_band_exact'])
        assert all(vhi<F(r['lowest_band_exact'][0]) for r in records if r is not lowest)
        unit=F(inp['electron_anchor_GeV'])/original_m
        baryons=[]
        for r in records:
            if r['N'] not in [3,9]:continue
            lo,hi=map(F,r['lowest_band_exact']);gap=[lo-vhi,hi-vlo]
            baryons.append(dict(N=r['N'],twice_physical_spin=r['twice_physical_spin'],
                gap_interval_exact=list(map(str,gap)),conditional_GeV_interval=list(map(float,[x*unit for x in gap]))))
        # Every exact highest-weight triquark source column defines a CAR
        # polynomial, normalized by its l1 norm to a contraction on full Fock.
        # Record all labels and all spectral outcomes, never selecting a hit.
        probes=[]
        _,vv,vframe=numerical[6,0]
        vacuum=vframe['R']@vv[:,0]
        inputwords=sectors[2]['words'];targetwords=sectors[3]['words']
        targetindex={w:i for i,w in enumerate(targetwords)}
        for three in sectors[1]['spin_frames']:
            twos=three['twos'];_,tv,tframe=numerical[9,twos]
            eigenframe=tframe['R']@tv
            for col in range(three['rank']):
                coefficients=three['R'][:,col]
                norm=sum(abs(int(x)) for x in coefficients)
                assert norm>0
                emitted=np.zeros(len(targetwords))
                for word,coefficient in zip(sectors[1]['words'],coefficients):
                    if not coefficient:continue
                    modes=[i for i in range(12) if word>>i&1]
                    assert len(modes)==3
                    for initial,amplitude in zip(inputwords,vacuum):
                        target=initial;sign=1
                        for mode in reversed(modes):
                            target,s=creation(mode,target)
                            if target is None:break
                            sign*=s
                        if target is not None:emitted[targetindex[target]]+=float(coefficient)*sign*amplitude/norm
                weights=abs(eigenframe.T@emitted)**2
                assert abs(weights.sum()-np.dot(emitted,emitted))<1e-10
                assert -1e-12<=weights.sum()<=1+1e-12
                probes.append(dict(twice_operator_spin=twos,highest_weight_polynomial_label=col,
                    exact_l1_normalization=norm,success_probability=float(weights.sum()),
                    failure_probability=float(1-weights.sum()),
                    complete_target_spectral_probabilities_diagnostic=weights.tolist(),
                    source_baryon_number_change=1,reservoir_occupancy_change_on_success=-1,
                    target_spin_selection_verified=True))
        results.append(dict(kappa=str(kappa),all_spin_sector_records=records,
            actual_vacuum_N=lowest['N'],actual_vacuum_twice_spin=lowest['twice_physical_spin'],
            baryon_sector_ground_gaps=baryons,complete_triquark_probe_records=probes))
    out=dict(schema='oph.obstruction_repair.spinful_two_cell.v1',
        carrier='two spatial cells times two physical spin components times three colours',
        full_CAR_dimension=4096,complete_SU3_singlet_dimensions=[s['singlet_dimension'] for s in sectors],
        Hamiltonian='m(N_left-N_right)-t sum_spin,color(c_left†c_right+h.c.)+kappa C_left²',
        spin_currents='S+=sum_cell,color c_(cell,up,color)†c_(cell,down,color), 2Sz=N_up-N_down',
        spin_current_algebra_exact=True,spin_H_commutators_exact=True,
        all_spin_multiplets_retained=True,t='1/8',actual_mass=str(original_m),mass_rounding_error_bound=str(mass_error),
        empirical_fit=False,results=results,
        physical_nucleon_or_neutron_proton_dictionary=False,weak_flavour_absent=True,
        native_M6_complete_program_executed=False,
        spatial_Dirac_or_continuum_transfer=False)
    (HERE/'SPINFUL_TWO_CELL_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['kappa','actual_vacuum_N','actual_vacuum_twice_spin','baryon_sector_ground_gaps']} for r in results],indent=2))


if __name__=='__main__':main()
