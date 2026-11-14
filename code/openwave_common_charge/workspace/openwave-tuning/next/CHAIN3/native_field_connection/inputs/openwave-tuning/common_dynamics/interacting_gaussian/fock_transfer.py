"""Controlled interacting photon-cap transfer by boundary flux.

Every positive oscillator is retained. The three harmonic canonical pairs
remain continuous and have identical dynamics in both compared sources.
"""
from fractions import Fraction as F
from itertools import combinations,product
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math

import branches as b
import certified_bounds as bound


def outward_dyadic(value,bits=160):
    value=F(value);den=1<<bits
    return F(-((-value.numerator*den)//value.denominator),den)


@lru_cache(None)
def axis_current_bounds(particles):
    if type(particles)is not int or particles not in (1,2):raise ValueError('one/two charges required')
    maxima=[]
    for axis in range(3):
        largest=0
        for cells in combinations(range(8),particles):
            for signs in product((-1,1),repeat=particles):
                J=sum((b.matter.route(cell,axis,sign)[2] for cell,sign in zip(cells,signs)),start=b.np.zeros(288,dtype=b.np.int64))
                largest=max(largest,int(J@J))
        maxima.append(largest)
    return tuple(maxima)


def parameters(macros,particles=1,initial_cells=1):
    c=bound.contract(macros,particles,initial_cells)
    wmin=F(c['coexact_frequency_lower_exact'])
    axis=axis_current_bounds(particles)
    gamma=[bound.sqrt_upper(F(j2,4)/wmin) for j2 in axis]
    length=macros*sum(gamma,F())
    return {'macros':macros,'particles':particles,'initial_cells':initial_cells,
            'all_branches':c['all_branches'],
            'axis_current_norm_squared_max':list(axis),
            'axis_Weyl_coefficient_norm_upper_exact':[str(x) for x in gamma],
            'cumulative_Weyl_length_upper_exact':str(length),
            'intermediate_occupation_upper_exact':str(length*length),
            'length_rule':'Sum every conditional Weyl interval, including inverse and repeated kicks. The occupation bound here is the square of that cumulative length, not a measured final or maximal occupation.'}


def contract(macros,particles,initial_cells,cap):
    if type(cap)is not int or cap<0:raise ValueError('nonnegative photon cap required')
    p=parameters(macros,particles,initial_cells)
    lam=F(p['intermediate_occupation_upper_exact']);B=p['all_branches']
    length=F(p['cumulative_Weyl_length_upper_exact'])
    if not lam:
        point=tail=error2=F()
    elif F(cap+2)<=lam:
        point=tail=F(1);error2=F(4)
    else:
        point=lam**(cap+1)/math.factorial(cap+1)
        tail=point*F(cap+2)/(F(cap+2)-lam)
        error2=min(F(4),2*B*tail+2*B*(cap+1)*length*length*point)
    dimension=math.comb(190+cap,cap)
    point_out=outward_dyadic(point);tail_out=outward_dyadic(tail);error_out=outward_dyadic(error2)
    assert point_out>=point and tail_out>=tail and error_out>=error2
    return {**p,'total_photon_cap':cap,
            'positive_series_point_upper_exact':str(point_out),
            'positive_series_tail_upper_exact':str(tail_out),
            'full_interacting_state_error_squared_upper_exact':str(error_out),
            'rounding':'Exact positive rational formulas evaluated before outward rounding to160binaryfractionbits; no floating input to the bound.',
            'unrounded_point_formula':'lambda_star^(R+1)/(R+1)!',
            'unrounded_tail_formula':'point*(R+2)/(R+2-lambda_star) when R+2>lambda_star',
            'unrounded_error_squared_formula':'min(4,2*B*tail+2*B*(R+1)*L^2*point)',
            'state_error_upper':float(bound.sqrt_upper(error_out,bits=80)),
            'coexact_Fock_dimension':str(dimension),
            'coexact_binary_wires':(dimension-1).bit_length(),
            'source':'P_R H0 P_R on all190oscillators; each routed conditional translation uses exp(i P_R L_j P_R), with exact harmonic Weyl factor and Coulomb/mass phases.',
            'native_scope':'Finite coexact quantum matrix and state-dependent transfer theorem. No emitted native circuit/reset law, physical local-edge encoder, harmonic regulator or interacting spatial continuum is included.',
            'quantum_reference':'Same branch argument tensors with arbitrary spectator references when initial coexact field is its vacuum; source clock factors identically.'}


def run():
    rows=[]
    for n,p,initial in ((1,1,8),(2,1,8),(2,2,1)):
        param=parameters(n,p,initial);lam=F(param['intermediate_occupation_upper_exact'])
        cap=4*math.ceil(lam)+128
        row=contract(n,p,initial,cap)
        assert F(row['full_interacting_state_error_squared_upper_exact'])<F(1,10**12)
        rows.append(row)
        assert contract(n,p,initial,1)['state_error_upper']==2.
    import numpy as np
    from scipy.linalg import expm
    controls=[]
    dim=65;ann=np.zeros((dim,dim),complex)
    for k in range(1,dim):ann[k-1,k]=math.sqrt(k)
    generator=.4*(ann+ann.conj().T)
    vacuum=np.eye(dim,dtype=complex)[:,0]
    for cap in (1,4,8):
        exact=expm(1j*generator)@vacuum
        finite=expm(1j*generator[:cap+1,:cap+1])@vacuum[:cap+1]
        lifted=np.pad(finite,(0,dim-cap-1))
        error=float(np.linalg.norm(exact-lifted))
        lam=F(4,25);point=lam**(cap+1)/math.factorial(cap+1)
        tail=point*F(cap+2)/(F(cap+2)-lam)
        error2=2*tail+2*(cap+1)*lam*point
        assert error*error<float(error2)+1e-27
        controls.append({'cap':cap,'numerical_state_error':error,'analytic_error_squared_upper_exact':str(error2)})
    # A closed Weyl rectangle has zero final displacement and bounded
    # intermediate occupation, while its scalar phase and path length add.
    gx=generator;gy=.4j*(ann.conj().T-ann)
    loop=expm(-1j*gy)@expm(-1j*gx)@expm(1j*gy)@expm(1j*gx)
    caploop=expm(-1j*gy[:2,:2])@expm(-1j*gx[:2,:2])@expm(1j*gy[:2,:2])@expm(1j*gx[:2,:2])
    looping=[]
    for repeats in (1,4,12):
        exact=np.linalg.matrix_power(loop,repeats)@vacuum
        finite=np.pad(np.linalg.matrix_power(caploop,repeats)@vacuum[:2],(0,dim-2))
        length=F(8,5)*repeats
        assert length*length>F(8,25)
        looping.append({'repetitions':repeats,'max_exact_intermediate_occupation':'8/25',
                        'final_displacement_exact':0,'cumulative_length_exact':str(length),
                        'retained_vacuum_scalar':[float(exact[0].real),float(exact[0].imag)],
                        'cap1_state_error':float(np.linalg.norm(exact-finite))})
    assert looping[-1]['cap1_state_error']>looping[0]['cap1_state_error']
    paths=[Path(__file__),Path(bound.__file__),Path(b.__file__),b.HERE.parent/'quadratic_field/MANIFEST.json',
           b.HERE.parent/'quadratic_field/FOCK.json',Path(b.matter.__file__)]
    pins={str(p.relative_to(b.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.interacting_gaussian.Fock_cap_transfer.v1',
            'checks':{'exact_actual_axis_current_bounds':True,'certified_positive_frequency_lower_bound':True,
                      'all_intermediate_coherent_occupations_bounded':True,
                      'boundary_flux_plus_terminal_tail':True,'all190modes_and_dimension_retained':True,
                      'cap1_transfer_is_vacuous':True,'independent_one_mode_boundary_controls':True,
                      'repeated_inverse_kicks_do_not_shrink_cumulative_length':True},
            'contracts':rows,'one_mode_controls':controls,'closed_loop_controls':looping,
            'proof':'PROOF.md, Interacting photon-cap dynamics',
            'source_sha256':pins,'official_rows_changed':[]}
    (b.HERE/'FOCK_TRANSFER.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':len(result['checks']),'caps':[r['total_photon_cap'] for r in rows],
                     'state_bounds':[r['state_error_upper'] for r in rows]}),flush=True)


if __name__=='__main__':run()
