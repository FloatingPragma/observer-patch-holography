"""Finite Gibbs targets, rational preparations and complete native flux records.

The declared source is the unchanged208-state native singlet/quantum-flux
control. Exact bounds and rational records are separated from floating
Gibbs diagnostics. No physical vacuum or spacetime limit is inferred.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
import scipy.linalg as la
import sympy as s

from model.particle_observer_outputs import field_route
from model.particle_local_gauge_state import state as local_state
from model.interacting_extension.chiral_lattice import source_data


SOURCE_INPUTS = (
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
    'model/particle_local_gauge_state/state.py',
    'model/particle_local_gauge_state/PROOF.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean',
    'Lean/Thermodynamics/FirstLawIdentity.lean',
    'Lean/Thermodynamics/LowTemperatureControl.lean',
)


def rational(value, name, positive=False):
    return local_state.rational(value, name, positive)


def integer(value, name, minimum=0):
    return local_state.integer(value, name, minimum)


def _sym(value):
    q = F(value)
    return s.Rational(q.numerator, q.denominator)


@lru_cache(maxsize=1)
def exact_control():
    native = field_route.exact_rotor(6)
    c = field_route.rotor_coefficients()
    shift = native['shift']; n = len(native['flux'])
    HE = s.kronecker_product(s.eye(16), s.diag(*native['flux'])**2/2)
    HW = (s.kronecker_product(c['onsite'], s.eye(n))
          +s.kronecker_product((c['hop']+c['hop'].H)/2, shift+shift.H))
    HD = s.kronecker_product((c['hop']-c['hop'].H)/2, shift-shift.H)
    if HE+HW+HD != native['H']:
        raise ArithmeticError('the thermal components changed the actual source')
    return {**native, 'HE': s.ImmutableMatrix(HE), 'HW': s.ImmutableMatrix(HW),
            'HD': s.ImmutableMatrix(HD), 'beta': F(1, 64), 'dimension': 208,
            'spectral_lower': -4, 'spectral_upper': 22}


@lru_cache(maxsize=1)
def connected_blocks():
    """Every source basis state belongs to exactly one exact reducing block."""
    h = exact_control()['H']; size = h.rows
    neighbors = [set() for _ in range(size)]
    for (i,j), value in h.todok().items():
        if i != j and value:
            neighbors[i].add(j); neighbors[j].add(i)
    unseen=set(range(size)); blocks=[]
    while unseen:
        seed=min(unseen); todo=[seed]; component=set()
        while todo:
            i=todo.pop()
            if i in component:continue
            component.add(i); todo.extend(neighbors[i]-component)
        unseen-=component; blocks.append(tuple(sorted(component)))
    if sorted(i for b in blocks for i in b) != list(range(size)):
        raise ArithmeticError('reducing blocks lost a source outcome')
    return tuple(blocks)


@lru_cache(maxsize=1)
def real_computation_chart():
    """Exact CAR basis change; restore every operator to source coordinates."""
    phases=tuple(s.I**(((i//13)>>2)&3).bit_count() for i in range(208))
    source=exact_control()['H']
    entries={ij:s.expand(s.conjugate(phases[ij[0]])*v*phases[ij[1]])
             for ij,v in source.todok().items()}
    if not all(v.is_Rational for v in entries.values()):
        raise ArithmeticError('the declared CAR phase chart is not real rational')
    h=s.ImmutableSparseMatrix(208,208,entries)
    if h!=h.T:raise ArithmeticError('source chart lost Hermiticity')
    return phases,h


def preparation_bound(order=12, beta='1/64'):
    m=integer(order, 'Taylor order'); b=rational(beta, 'inverse temperature', True)
    # -4 <= H <=22, hence0 <= beta(22I-H)/2 <=13 beta.
    x=13*b
    if x >= m+2:
        raise ValueError('increase Taylor order for this inverse temperature')
    tail=x**(m+1)/math.factorial(m+1)/(1-x/F(m+2))
    return {'order':m,'beta':b,'positive_amplitude_argument_upper':x,
            'amplitude_operator_error_upper':tail,
            'canonical_state_half_trace_error_upper':min(F(1),2*tail),
            'target_source_spectral_interval':[-4,22]}


def prepared_density(order=12, beta='1/64'):
    bound=preparation_bound(order,beta)
    return _prepared_density(bound['order'],bound['beta'])


@lru_cache(maxsize=16)
def _prepared_density(order, beta):
    bound=preparation_bound(order,beta); b=_sym(bound['beta'])
    phases,h=real_computation_chart(); size=h.rows
    entries={}; normalizer=s.S.Zero; block_data=[]
    for indices in connected_blocks():
        hb=h.extract(indices,indices); eye=s.eye(len(indices))
        x=b*(22*eye-hb)/2
        q=eye; term=eye
        for j in range(1,order+1):
            term=term*x/j; q=q+term
        q=s.ImmutableMatrix(q); sq=q*q
        normalizer+=s.trace(sq)
        for (i,j),value in sq.todok().items():
            if value:entries[indices[i],indices[j]]=phases[indices[i]]*value*s.conjugate(phases[indices[j]])
        block_data.append((indices,q))
    rho=s.ImmutableSparseMatrix(size,size,{ij:v/normalizer for ij,v in entries.items()})
    if s.trace(rho)!=1 or rho!=rho.H:
        raise ArithmeticError('positive polynomial state is not normalized')
    return {'rho':rho,'amplitude_blocks':tuple(block_data),'norm_squared':normalizer,
            'computational_basis_phases':phases,
            'preparation_trace_error_upper':bound['canonical_state_half_trace_error_upper'],
            'bound':bound}


def cayley_bound(duration='1/64', steps=4):
    t=rational(duration, 'duration'); j=integer(steps, 'Cayley step count', 1)
    return min(F(2),t**3*22**3/(12*j*j))


def cayley_parts(indices, duration='1/64', steps=4):
    j=integer(steps, 'Cayley step count', 1); t=rational(duration,'duration')
    return _cayley_parts(tuple(indices),t,j)


@lru_cache(maxsize=256)
def _cayley_parts(indices, duration, steps):
    j=steps; t=_sym(duration)
    h=real_computation_chart()[1].extract(indices,indices); eye=s.eye(h.rows)
    a=t/(2*j); h2=h*h
    denominator=(eye+a*a*h2).inv()
    re=(eye-a*a*h2)*denominator; im=-2*a*h*denominator
    # Real rational pairs avoid repeated symbolic complex denominator work.
    rr,ii=eye,s.zeros(h.rows); power=j
    while power:
        if power & 1:rr,ii=rr*re-ii*im,rr*im+ii*re
        power//=2
        if power:re,im=re*re-im*im,re*im+im*re
    return s.ImmutableMatrix(rr),s.ImmutableMatrix(ii)


def cayley_block(indices, duration='1/64', steps=4):
    re,im=cayley_parts(indices,duration,steps)
    phases,_=real_computation_chart()
    values={(i,j):s.expand(phases[indices[i]]*(re[i,j]+s.I*im[i,j])*s.conjugate(phases[indices[j]]))
            for i in range(len(indices)) for j in range(len(indices))}
    return s.ImmutableMatrix(len(indices),len(indices),lambda i,j:values[i,j])


@lru_cache(maxsize=4)
def record_law(order=12, duration='1/64', steps=4):
    prep=prepared_density(order); flux=exact_control()['flux']; n=len(flux)
    weights={(f,g):s.S.Zero for f in flux for g in flux}
    for indices,q in prep['amplitude_blocks']:
        ure,uim=cayley_parts(indices,duration,steps)
        for f in sorted({flux[i%n] for i in indices}):
            rows=[r for r,i in enumerate(indices) if flux[i%n]==f]
            selected=s.zeros(q.rows)
            for r in rows:selected[r,:]=q[r,:]
            after_re=ure*selected; after_im=uim*selected
            for r,i in enumerate(indices):
                amount=sum((zr*zr+zi*zi for zr,zi in zip(after_re[r,:],after_im[r,:])),s.S.Zero)
                weights[f,flux[i%n]]+=amount
    rows=[{'first_flux':f,'second_flux':g,
           'original_probability':s.factor(p/prep['norm_squared'])}
          for (f,g),p in weights.items()]
    if any(r['original_probability']<0 for r in rows) or sum(r['original_probability'] for r in rows)!=1:
        raise ArithmeticError('complete flux records lost original normalization')
    return {'first_flux_labels':list(flux),'second_flux_labels':list(flux),
            'all_joint_labels':n*n,'zero_probability_labels':sum(r['original_probability']==0 for r in rows),
            'records':rows,'source_dimension':208,'duration':rational(duration,'duration'),
            'measurement_times':[F(0),rational(duration,'duration')],
            'Cayley_steps':steps,'exact_compiled_unitarity':True,
            'Gibbs_exact_evolution_output_half_trace_error_upper':min(F(1),prep['preparation_trace_error_upper']+cayley_bound(duration,steps)),
            'all_occupations_and_reference_retained':True}


def mixed_ground_obstruction(side=2, spacing='1'):
    length=integer(side,'periodic side',2); a=rational(spacing,'spacing',True)
    source=source_data.load()
    if len(source['charges'])!=15 or sum(source['charges'])!=0:
        raise ArithmeticError('the native determinant-one matter representation changed')
    return {'periodic_side':length,'cells':length**3,'spacing':a,
            'native_internal_dimension':15,'negative_beta_spin_dimension':2,
            'full_occupied_modes_per_site':30,
            'native_U1_charge_trace':sum(source['charges']),
            'tracial_Fock_kinetic_mean':F(0),
            'alternative_kinetic_mean':-F(750*length**3,17)/a,
            'same_H_energy_decrease':F(750*length**3,17)/a,
            'same_gauge_density_electric_and_plaquette_terms':True,
            'source_regulator':'degree-zero Legendre field h0; not a smooth-field coefficient',
            'finite_mixed_state_is_ground':False,
            'physical_vacuum_continuum_impossibility_claimed':False}


@lru_cache(maxsize=1)
def thermal_diagnostics():
    """Floating diagnostics only; exact positivity/bounds are separate proofs."""
    control=exact_control(); h=np.array(control['H'].tolist(),dtype=complex)
    beta=float(control['beta']); energy,v=la.eigh(h)
    weights=np.exp(-beta*(energy-energy.min())); weights/=weights.sum()
    rho=(v*weights)@v.conj().T
    average=float(weights@energy); variance=float(weights@(energy-average)**2)
    entropy=float(-weights@np.log(weights))
    e=np.array(control['E'].tolist(),dtype=complex)
    current=np.array(control['current'].tolist(),dtype=complex)
    comm=h@e-e@h
    return {'kind':'floating diagnostics, not interval certificates',
            'full_spectral_rows':len(energy),'minimum_energy':float(energy.min()),
            'maximum_energy':float(energy.max()),'mean_energy':average,
            'energy_variance':variance,'von_Neumann_entropy':entropy,
            'beta_heat_capacity':beta*beta*variance,
            'flux_mean':float(np.trace(rho@e).real),
            'current_mean':float(np.trace(rho@current).real),
            'current_second_moment':float(np.trace(rho@current@current).real),
            'Gibbs_stationarity_residual':float(la.norm(h@rho-rho@h)),
            'current_operator_nonzero':bool(la.norm(current)>0),
            'source_current_identity_residual':float(la.norm(1j*comm-current))}


def serial(value):
    if isinstance(value,(F,s.Basic)):return str(value)
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [serial(v) for v in value]
    return value


def certificate():
    return serial({'schema':'oph.field_state_response.v1','source_inputs':list(SOURCE_INPUTS),
        'finite_state_preparation':{
            'source':'unchanged field_route.exact_rotor(6)', 'dimension':208,
            'reducing_block_dimension_counts':dict(sorted(Counter(map(len,connected_blocks())).items())),
            'all_Fock_occupation_counts':list(range(5)),
            'declared_inverse_temperature':F(1,64),
            'positive_polynomial':preparation_bound(),
            'canonical_unitary_preparation_requires_no_postselection':True,
            'joint_Gauss_and_parity_preserved':True,
            'reference_H':0,
            'source_density_stationary':True,'joint_purification_stationary':False},
        'complete_flux_records':record_law(),
        'source_thermal_identities':{
            'exact_identity_state':'target rho_beta=exp(-beta H)/Tr exp(-beta H)',
            'Gibbs_modular_H':'-log(rho_beta)=beta H+log Z I',
            'first_law':'dE=Tr(H d rho)+Tr(rho dH)',
            'source_KMS':'Tr rho A alpha_(t+i beta)(B)=Tr rho alpha_t(B) A',
            'source_cyclic_passivity':True,
            'arbitrary_reference_operations_passive':False,
            'heat_capacity':'beta^2 Var_rho(H)',
            'fixed_beta_free_energy_derivative':'dF/ds=<dH/ds>',
            'strain_entropy_derivative':'dS/ds=-beta^2 Cov_KM(H,dH/ds)',
            'polynomial_preparation_exactly_Gibbs':False,
            'polynomial_preparation_stationary_and_passive':True,
            'polynomial_preparation_complete_passivity_claimed':False},
        'mixed_spatial_ground_obstruction':mixed_ground_obstruction(),
        'diagnostic_Gibbs_values':thermal_diagnostics(),
        'physical_vacuum_or_continuum_derived':False,
        'physical_temperature_energy_or_mass_calibrations':[]})
