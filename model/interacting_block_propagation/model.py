"""Exact gapped block architecture and its moving collective excitation band.

The blocks are the actual frozen eight-cell sources. Inter-block couplings
are explicitly declared spectral dressings of their native spin currents.
No result is transferred to the uniform larger original lattice.
"""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]

import mpmath as mp
import numpy as np
import sympy as s
from model.coherent_defect_wave.chiral_joint.collective_spin_response.model import operator_certificate

BASE='model/coherent_defect_wave/chiral_joint/'
PARENTS=[BASE+name+'/results.json' for name in
         ('atomic_spin_geometry','collective_spin_response','interacting_ground_preparation')]


def decimal(value):
    value=F(value)
    return mp.mpf(value.numerator)/value.denominator


def load_parents():
    results={}
    pins={}
    for name in PARENTS:
        value=json.loads((RER/name).read_text())
        for path,expected in value['input_sha256'].items():
            actual=hashlib.sha256((RER/path).read_bytes()).hexdigest()
            assert actual==expected,('immutable parent changed',path)
            pins[path]=actual
        pins[name]=hashlib.sha256((RER/name).read_bytes()).hexdigest()
        results[name.split('/')[-2]]=value
    return results,pins


def bond_identity():
    # A fifth state represents an arbitrary untouched higher block sector.
    # The proof below acts identically on every such sector in the full block.
    d=5
    lower=[]
    for a in (1,2,3):
        b=s.zeros(d);b[0,a]=1;lower.append(b)
    pt=s.diag(0,1,1,1,0)
    q=s.diag(0,1,1,1,1)
    swap=sum((s.kronecker_product(b.H,b)+s.kronecker_product(b,b.H) for b in lower),s.zeros(d*d))
    number=s.kronecker_product(pt,s.eye(d))+s.kronecker_product(s.eye(d),pt)
    fullnumber=s.kronecker_product(q,s.eye(d))+s.kronecker_product(s.eye(d),q)
    assert swap==swap.H
    assert swap*fullnumber==fullnumber*swap
    spectra={}
    for sign in (-1,1):
        eigen=(number+sign*swap).eigenvals()
        assert all(value>=0 for value in eigen)
        spectra[str(sign)]={str(k):v for k,v in eigen.items()}
    assert all(value in (-1,0,1) for value in swap.eigenvals())
    assert all(swap[i,0]==0 for i in range(d*d))
    return {'dimension':d*d,'higher_sector_retained':True,
            'P_triplet_x_plus_P_triplet_y_plus_or_minus_exchange_positive':spectra,
            'exchange_norm':'1','vacuum_annihilated':True,'total_excited_block_number_conserved':True}


def adjacency(side,dimension=3):
    shape=(side,)*dimension
    n=side**dimension
    matrix=np.zeros((n,n),int)
    for site in np.ndindex(shape):
        i=np.ravel_multi_index(site,shape)
        for axis in range(dimension):
            for sign in (-1,1):
                target=list(site);target[axis]=(target[axis]+sign)%side
                matrix[np.ravel_multi_index(tuple(target),shape),i]+=1
    assert np.array_equal(matrix,matrix.T)
    assert np.all(matrix.sum(axis=0)==2*dimension)
    return matrix


def momentum_check(side):
    a=adjacency(side)
    sites=np.asarray(list(np.ndindex((side,)*3)))
    maximum=0.
    energies=[]
    for menu in np.ndindex((side,)*3):
        k=2*np.pi*np.asarray(menu)/side
        vector=np.exp(1j*sites@k)/np.sqrt(len(sites))
        eigen=2*sum(np.cos(k))
        maximum=max(maximum,float(np.max(abs(a@vector-eigen*vector))))
        energies.append(float(eigen))
    assert maximum<1e-12
    # e^{-i Delta t} is a common phase and cancels in position probabilities.
    tau=.25  # tau=J*t, a source-defined dimensionless clock reading
    eigvals,eigvec=np.linalg.eigh(a.astype(float))
    wave=eigvec@(np.exp(1j*tau*eigvals)*eigvec[0,:])
    probabilities=abs(wave)**2
    assert abs(sum(probabilities)-1)<1e-12
    return {'side':side,'block_count':side**3,'physical_source_cells':8*side**3,
            'complete_one_excitation_states_including_three_polarizations':3*side**3,
            'Fourier_eigenvector_max_error':maximum,
            'adjacency_spectrum_minmax':[min(energies),max(energies)],
            'clock_J_times_t':tau,'complete_position_probabilities':probabilities.tolist(),
            'position_probability_sum':float(sum(probabilities))}


def preparation_scaling(preparation):
    trial=preparation['atomic_cube_trial']
    filt=preparation['finite_energy_agreement_filter']
    p0=F(trial['exact_interacting_ground_trial_overlap_lower'])
    epsilon=F(filt['mixed_ground_excited_amplitude_upper'])
    r0=(1-p0)/p0
    contraction=r0+2*epsilon**2
    assert 0<r0<1 and 0<contraction<1
    tolerance=F(1,2**20)
    rows=[]
    for count in (1,8,64,1000,10**6,10**12):
        depth=0
        odds=r0
        while count*odds>tolerance**2:
            depth+=1;odds*=contraction
        # Independent successful subtrees prepare independent source marginals.
        # All failures remain part of the complete finite experiment.
        log_success=mp.mpf(0)
        for j in range(depth):
            odds_j=decimal(r0*contraction**j)
            log_success-=count*2**(depth-j)*mp.log1p(odds_j)
        rows.append({'blocks':count,'physical_cells':8*count,'binary_filter_depth':depth,
                     'source_trials_per_block':2**depth,
                     'one_block_excited_odds_upper':str(odds),
                     'whole_volume_trace_distance_squared_upper':str(count*odds),
                     'whole_volume_trace_distance_upper_decimal':str(mp.sqrt(decimal(count*odds))),
                     'log_all_success_probability_lower':str(log_success),
                     'all_failure_branches_retained':True})
    return {'one_block_trial_ground_overlap_lower':str(p0),
            'same_frozen_filter_M':filt['explicit_integer_M'],
            'same_frozen_filter_leakage_upper':str(epsilon),
            'initial_excited_odds_upper':str(r0),'geometric_odds_contraction':str(contraction),
            'recursion':'r_next <= r^2+2 epsilon^2 r <= (r0+2 epsilon^2) r',
            'target_global_trace_distance':str(tolerance),'rows':rows,
            'efficient_preparation_claim':False}


@lru_cache(maxsize=1)
def calculate():
    mp.mp.dps=70
    parents,pins=load_parents()
    native=operator_certificate()
    collective=parents['collective_spin_response']['same_source_interacting_records']
    line=collective['actual_interacting_energy_line']
    control=collective['complete_full_source_operator_controls']
    gaplo,gaphi=map(F,line['excitation_energy_interval_in_source_units'])
    rlo,rhi=map(F,line['first_line_probability_conditioned_on_ground_preparation'])
    nextlo=F(line['next_low_band_energy_separation_lower'])
    g=F(control['atomic_gap_lower'])
    v=F(control['exact_rational_error_bounds']['full_hopping_perturbation_norm'])
    assert collective['source']['hopping_coefficient']=='1/17179869184'
    z=6
    coupling=gaplo/(4*z*rhi)
    jlo,jhi=coupling*rlo,coupling*rhi
    uniform_gap=gaplo-z*jhi
    bandhi=gaphi+z*jhi
    higher_one_block=min(gaplo+nextlo,g-2*v)
    rest_lower=min(2*uniform_gap,higher_one_block)
    separation=rest_lower-bandhi
    assert uniform_gap==3*gaplo/4 and separation>0
    bond=bond_identity()
    momenta=[momentum_check(side) for side in (3,4,6)]
    # Taylor and walk-tail controls use no diagonalization or empirical input.
    continuation=[]
    for spacing in (F(1),F(1,2),F(1,4),F(1,8),F(1,16)):
        bound=spacing*spacing*F(3,12)  # |k_i|<=1, rescaled time=1
        maximum=0
        for k in np.linspace(-1,1,101):
            difference=abs(2*(1-np.cos(float(spacing)*k))/float(spacing)**2-k*k)
            maximum=max(maximum,difference)
        assert maximum<=float(spacing*spacing/F(12))+1e-12
        continuation.append({'spacing':str(spacing),'momentum_box':'|k_i|<=1',
                             'rescaled_time':'1','Schrodinger_generator_error_upper':str(bound),
                             'one_axis_numeric_error_max':maximum})
    walk=[]
    for radius in (4,8,16,32):
        u=mp.mpf(6)/4
        tail=2*sum(u**n/mp.factorial(n) for n in range(radius,300))
        exponential=2*mp.exp(mp.e*u-radius)
        assert tail<=exponential
        walk.append({'distance_to_removed_boundary':radius,'clock_J_times_t':'1/4',
                     'wavefunction_error_upper':str(min(mp.mpf(2),tail)),
                     'exponential_error_upper':str(min(mp.mpf(2),exponential))})
    preps=preparation_scaling(parents['interacting_ground_preparation']['preparation'])
    result={'schema':'oph.codex.block_lattice_refinement.v1',
            'native_block_operator_replayed':native,
            'immutable_parent_pins_validated':PARENTS,
            'construction':{'block':'actual eight-cell H0+2^-34 T, all native sectors retained',
                            'physical_sector':'tensor product of the complete block Gauss kernels; kinematic complement retained',
                            'block_ground':'unique frozen physical ground',
                            'lowering':'B_a=P_g O_a P_t/60',
                            'native_current':'O_a=sum_x (-1)^x_a c_x^dag(sigma_a/2 tensor I15)c_x',
                            'interblock_interaction':'-kappa sum_<xy>,a(B_xa^dag B_ya+B_ya^dag B_xa)',
                            'clock':'exp(-it H_block_lattice)',
                            'projector_dressing_is_declared_architecture':True,
                            'equal_to_uniform_larger_original_source':False,
                            'native_block_and_spectral_services_are_efficiently_compiled':False,
                            'measured_inputs':[]},
            'source_bounds':{'block_gap_interval':[str(gaplo),str(gaphi)],
                             'r_interval':[str(rlo),str(rhi)],'next_block_level_separation_lower':str(nextlo),
                             'maximum_graph_degree':z,'kappa_definition':'gap_lower/(4*z*r_upper)',
                             'kappa':str(coupling),'effective_hopping_J_interval':[str(jlo),str(jhi)],
                             'volume_uniform_gap_lower':str(uniform_gap),
                             'one_excitation_band_upper':str(bandhi),
                             'all_other_excitation_energies_lower':str(rest_lower),
                             'isolated_one_excitation_band_separation_lower':str(separation),
                             'decimals':{'block_gap_lower':str(decimal(gaplo)),'J_lower':str(decimal(jlo)),
                                         'J_upper':str(decimal(jhi)),'uniform_gap_lower':str(decimal(uniform_gap)),
                                         'band_isolation_lower':str(decimal(separation))}},
            'exact_bond_control':bond,
            'uniform_theorem':{'product_ground_unique_every_finite_volume':True,
                               'full_complement_kept':True,'zero_interblock_ground_entanglement':True,
                               'one_excitation_reducing_subspace_exact':True,
                               'dispersion':'E_a(k)=Delta-2 J sum_i cos(k_i), a=0,1,2',
                               'finite_block_momentum_is_internal_X_momentum':False,
                               'band_minimum':'Delta-6J','one_excitation_band_rank':'3*number_of_blocks',
                               'interaction_is_a_free_full_many_body_model':False},
            'finite_momentum_and_records':momenta,
            'complete_localized_excitation_record_law':{
                'pump_success_probability':'r','pump_failure_probability':'1-r',
                'position_and_retained_axis_probability':'r |U_y0(t)|^2/3',
                'exact_completeness':'(1-r)+sum_y,a r |U_y0(t)|^2/3=1',
                'all_axis_position_and_failure_outcomes_retained':True},
            'volume_boundary_error_controls':walk,'preparation_scaling':preps,
            'continuum_controls':{'fixed_source_clock_kinetic_coefficient':'J*a^2 -> 0',
                                  'nontrivial_low_energy_clock':'t=eta/(J*a^2), subtract constant Delta-6J',
                                  'limit_generator':'sum_i k_i^2',
                                  'error':'|E_a(k)-E_a(0)-J*a^2|k|^2| <= J*a^4 sum_i k_i^4/12',
                                  'fixed_band_momentum_error_controls':continuation,
                                  'maintaining_J_proportional_a^-2_violates_fixed_block_gap_control':True,
                                  'Lorentz_covariance_derived':False},
            'input_sha256':pins}
    return result
