"""Exact source moments and certified complete charge-channel survival."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from pathlib import Path
import numpy as np
from model.protected_defect_catalog import charge
from model.native_ordered_ground import source, certificate
from model.centered_particle_spectrum.certificate import sqrt_bounds, decimal_bounds

ROOT = Path(__file__).resolve().parents[2]
WITNESS = 'model/native_ordered_ground/witness.npz'
ORDER = 160


def spectral_width_certificate(op):
    h, scale = op['H'], op['scale']
    size = len(h)
    identity = np.eye(size, dtype=np.int64)
    indices = np.arange(size)
    h0 = np.zeros_like(h)
    for x, y in source.cube()[0]:
        projector = np.zeros_like(h)
        for word in range(size):
            if ((word >> x) ^ (word >> y)) & 1:
                projector[word, word] = scale//2
                projector[word ^ (1 << x) ^ (1 << y), word] = scale//2
        if not np.array_equal(projector@projector, scale*projector):
            raise ArithmeticError('each source bond must be a negative projector')
        h0 -= projector
    field = sum((identity[indices ^ (1 << x)] for x in range(8)), np.zeros_like(h))
    if not np.array_equal(h0-field*(scale//8), h):
        raise ArithmeticError('selected field and all twelve bonds must reconstruct source H')
    return {'negative_unit_projector_bonds':12, 'field_operator_norm_bound':'1',
            'full_H_spectrum_enclosure':['-13','1'], 'full_spectral_width_bound':'14',
            'H_plus_6I_operator_norm_bound':'7'}


def instrument_certificate(op):
    identity = np.eye(256, dtype=np.int64)
    x = identity[np.arange(256)^1]
    completeness = np.zeros_like(identity)
    for q in (-1, 1):
        a = op['twice_creation'][q]
        expected = op['local_Z']@(identity-q*op['signs'][0]*x)
        if not np.array_equal(a, expected):
            raise ArithmeticError('charge branch must equal local X pointer followed by Z')
        completeness += a.T@a
    if not np.array_equal(completeness, 4*identity):
        raise ArithmeticError('both creation outcomes are required')
    return {'Kq':'Z0(I-q*epsilon0*X0)/2', 'outcomes':[-1,1],
            'Kminus_adjoint_Kminus_plus_Kplus_adjoint_Kplus_equals_I':True,
            'fine_all_site_X_measurement_used_to_prepare_branch':False,
            'conditional_exterior_marginals_claimed_unchanged':False}


@lru_cache(maxsize=1)
def source_data():
    op = charge.operators()
    with np.load(ROOT/WITNESS, allow_pickle=False) as archive:
        witness = {key:archive[key] for key in archive.files}
    checked = certificate.verify(witness)
    g = witness['vectors'][:,0].astype(object)
    ground_error = F(checked['normalized_ground_vector_error'])
    norm = sum(int(x)**2 for x in g)
    width = spectral_width_certificate(op)
    instrument = instrument_certificate(op)
    hc = op['H'].astype(object)+6*op['scale']*np.eye(256,dtype=object)
    sparse = [[(j,int(row[j])) for j in np.flatnonzero(row)] for row in hc]
    result = {}
    unit = 1 << certificate.BITS
    for q in (-1,1):
        image = op['twice_creation'][q].astype(object)@g
        image_norm = int(image@image)
        weight = F(image_norm,4*norm)
        delta = 4*ground_error*sqrt_bounds(F(norm,image_norm))[1]
        vector = [int(x) for x in image]
        moments = []
        for n in range(ORDER+1):
            moments.append(F(sum(int(a)*b for a,b in zip(image,vector)), image_norm*op['scale']**n))
            if n != ORDER:
                vector = [sum(a*vector[j] for j,a in row) for row in sparse]
        coefficients = witness['vectors'].astype(object).T@image
        weights = tuple(F(int(x)**2,unit**2*image_norm) for x in coefficients)
        energies = tuple(F(int(x),unit)+6 for x in witness['energies'])
        gram = F(checked['Gram_defect'])
        if abs(sum(weights)-1) > gram:
            raise ArithmeticError('complete spectral weight exceeds certified Gram defect')
        variance = moments[2]-moments[1]**2
        var_error = 6*7**2*delta
        if variance-var_error <= 0:
            raise ArithmeticError('nonzero within-charge energy variance required')
        result[q] = {'weight':weight, 'weight_interval':(weight-2*ground_error,weight+2*ground_error),
            'vector_error':delta, 'moments':tuple(moments), 'weights':weights, 'energies':energies,
            'variance_interval':(variance-var_error,variance+var_error),
            'spectral_radius':F(checked['complete_polar_Weyl_radius']), 'gram':gram}
    if sum(result[q]['weight'] for q in result) != 1:
        raise ArithmeticError('complete original charge weights must sum to one')
    return result, checked, width, instrument


def polynomial(moments,t,order):
    if type(t) is not F or t < 0 or type(order) is not int or not 0 <= order < len(moments):
        raise ValueError('nonnegative rational time and available Taylor order required')
    re = im = F()
    for n in range(order+1):
        term = moments[n]*t**n/F(factorial(n))
        if n%2 == 0:
            re += (-1)**(n//2)*term
        else:
            im -= (-1)**((n-1)//2)*term
    return re,im


def modulus_interval(re,im,error):
    if error < 0:
        raise ValueError('nonnegative amplitude error required')
    lower,upper = sqrt_bounds(re*re+im*im)
    return max(F(),lower-error)**2, min(F(1),(upper+error)**2)


def survival(row,t,order=32):
    re,im = polynomial(row['moments'],t,order)
    remainder = (7*t)**(order+1)/factorial(order+1)
    return modulus_interval(re,im,remainder+2*row['vector_error'])


def complete_spectral_survival(row,t,order=32):
    if type(t) is not F or t<0 or not 0<=order<=ORDER:
        raise ValueError('valid exact time and series order required')
    moments = tuple(sum(w*e**n for w,e in zip(row['weights'],row['energies']))
                    for n in range(order+1))
    re,im = polynomial(moments,t,order)
    energy_bound = max(abs(e) for e in row['energies'])
    remainder = sum(row['weights'])*(energy_bound*t)**(order+1)/factorial(order+1)
    error = remainder+t*row['spectral_radius']+row['gram']*(2+row['gram'])+2*row['vector_error']
    return modulus_interval(re,im,error)


def first_threshold(row,threshold=F(99,100),bits=24):
    if type(threshold) is not F or not 0<threshold<1 or type(bits) is not int or not 1<=bits<=30:
        raise ValueError('strict probability threshold and bounded precision required')
    lo,hi = F(),F(3,14)
    if survival(row,hi)[1] >= threshold:
        raise ArithmeticError('threshold must be crossed in the certified monotonic interval')
    for _ in range(bits):
        mid = (lo+hi)/2
        lower,upper = survival(row,mid)
        if lower>threshold:lo=mid
        elif upper<threshold:hi=mid
        else:break
    if not (survival(row,lo)[0]>threshold and survival(row,hi)[1]<threshold):
        raise ArithmeticError('strict endpoint enclosure required')
    return lo,hi


def echo_error(weight_floor,event_error,weight_error):
    """Error in R/w^2, when true and estimated weights are at least the floor."""
    if not 0<weight_floor<=1 or event_error<0 or weight_error<0:
        raise ValueError('positive probability floor and nonnegative errors required')
    return event_error/weight_floor**2+2*weight_error/weight_floor**3


def build():
    rows,checked,width,instrument = source_data()
    channels = []
    for q,row in rows.items():
        lo,hi = first_threshold(row)
        times = (F(),lo,hi,F(1,5))
        observations = []
        for t in times:
            interval = survival(row,t)
            spectral = complete_spectral_survival(row,t)
            if max(interval[0],spectral[0])>min(interval[1],spectral[1]):
                raise ArithmeticError('complete spectral and exact source-moment bounds disagree')
            wl,wh = row['weight_interval']
            ret=(wl*wl*interval[0],wh*wh*interval[1])
            observations.append({'time':str(t),'conditional_survival_interval':decimal_bounds(*interval,12),
                'independent_all_level_spectral_interval':decimal_bounds(*spectral,12),
                'original_charge_and_ground_return_interval':decimal_bounds(*ret,12),
                'same_charge_and_nonreturn_interval':decimal_bounds(max(F(),wl-ret[1]),wh-ret[0],12)})
        channels.append({'internal_source_charge':q,'sector_dimension':56,
            'original_preparation_probability_interval':decimal_bounds(*row['weight_interval'],12),
            'normalized_branch_vector_error':str(row['vector_error']),
            'within_charge_energy_variance_interval':decimal_bounds(*row['variance_interval'],12),
            'retention_threshold':'99/100', 'first_threshold_time_interval':decimal_bounds(lo,hi,10),
            'exact_first_threshold_time_bracket':[str(lo),str(hi)],
            'complete_spectral_weight_proposal':[str(w) for w in row['weights']],
            'complete_spectral_weight_sum':str(sum(row['weights'])),
            'observations':observations,
            'long_time_average_survival_lower_bound':'1/56',
            'persistent_uniform_site_binary_score':{'mean_conditioned_on_charge':str(F(q,4)),
                'plus_probability_conditioned_on_charge':str(F(1,2)+F(q,8)),
                'original_charge_and_plus_probability_interval':decimal_bounds(row['weight_interval'][0]*(F(1,2)+F(q,8)),row['weight_interval'][1]*(F(1,2)+F(q,8)),12),
                'original_charge_and_minus_probability_interval':decimal_bounds(row['weight_interval'][0]*(F(1,2)-F(q,8)),row['weight_interval'][1]*(F(1,2)-F(q,8)),12),
                'original_charge_and_plus_reference_contrast_magnitude_lower':str(row['weight_interval'][0]/8)},
            'physical_particle_lifetime_or_width':None})
    rebound = {str(t):decimal_bounds(*survival(rows[-1],t,ORDER),12) for t in (F(2),F(5))}
    if F(rebound['5'][0])<=F(rebound['2'][1]):
        raise ArithmeticError('finite same-charge survival rebound must be certified')
    return {'source_spectral_certificate':checked,'source_width_certificate':width,
        'local_creation_instrument':instrument,
        'all_complete_repair_charge_branches':charge.full_source_charge(),
        'complete_source_level_count':256,
        'energy_proposal_H_plus_6I':[str(e) for e in rows[1]['energies']],
        'source_energy_shift_changes_survival':False,
        'channels':channels,
        'charge_minus_survival_rebound':rebound,
        'persistent_contrast_readout':{'setting_distribution':'uniform over all eight native sites',
            'readout':'local X eigenvalue s, recorded score epsilon_x*s',
            'complete_charge_setting_and_binary_outcome_count':32,
            'reference_ground_plus_probability':'1/2',
            'local_contrast_sum_rule':'sum_x epsilon_x (Pr_q[X_x=+]-Pr_ground[X_x=+])=q',
            'some_site_conditional_contrast_at_least':'1/8',
            'uniform_site_original_weight_contrast_strictly_greater_than':'1/50',
            'fixed_site_contrast_persistence_asserted':False,
            'propagating_particle_identified':False},
        'normalized_echo_error_example':str(echo_error(F(4,25),F(1,10**8),F(1,10**8))),
        'scope':{'same_selected_ground_Hamiltonian_and_clock':True,
            'same_finite_total_charge_retained_at_all_times':True,
            'complete_two_outcome_pure_charge_preparation_constructed':True,
            'complete_echo_ground_return_and_complement_constructed':True,
            'declared_native_coherent_service_required_for_H_evolution_and_inverse_preparation':True,
            'full_charge_echo_tape_executed':False,
            'source_moments_and_all_spectral_levels_certified':True,
            'finite_recurrence_claimed_as_physical_decay':False,
            'continuum_decay_or_resonance_proved':False,
            'physical_particle_species_or_mass_assigned':False}}
