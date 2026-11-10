"""Exact finite two-wave collision in an invariant source dark-mode sector.

The attractive contact is a consequence of the declared positive continuous
clock. The original dark Floquet phase (-1)^N is not assigned this energy law.
"""
from fractions import Fraction as F
from itertools import combinations,product
from functools import lru_cache
import sympy as s

from .model import PairFockState,MODES,mode,pair_energy_action,pair_onsite
from ..species.spatial_fock import STAR,gamma,number_measurement
from ..species.joint_readout import constraint_preparation


OCCUPATIONS=tuple(combinations(range(4),2))


def source_column(occupied):
    """Wedge of four declared dark orbitals, with both Gauss pairs retained.

    Orbital order is (left a,left b,right a,right b). Each orbital is the
    address difference at port2 and internal a=0 or b=1.
    """
    occupied=tuple(occupied)
    if tuple(sorted(set(occupied)))!=occupied or any(type(i) is not int or i not in range(4) for i in occupied):
        raise ValueError('an ordered subset of four dark orbitals is required')
    output={}
    for links in product((1,1^STAR),repeat=2):
        for addresses in product((0,1),repeat=len(occupied)):
            actual=tuple(sorted((orbital//2)*MODES+mode(2,orbital%2,0,address)
                                for orbital,address in zip(occupied,addresses,strict=True)))
            output[links,actual]=F((-1)**sum(addresses))
    return PairFockState(output,4*2**len(occupied))


def source_combination(coefficients):
    if len(coefficients)!=6:
        raise ValueError('six coefficients on the N=2 dark-sector basis required')
    output={}
    for occupied,coefficient in zip(OCCUPATIONS,coefficients,strict=True):
        for key,value in source_column(occupied).amplitudes:
            output[key]=s.expand(output.get(key,0)+coefficient*value)
    return PairFockState(output,16)


def inner(left,right):
    if left.norm_denominator!=right.norm_denominator:
        raise ValueError('one common source normalization required')
    lookup=left.as_dict()
    return s.expand(sum((s.conjugate(lookup.get(key,0))*value for key,value in right.amplitudes),s.Integer(0))
                    /left.norm_denominator)


def expected_matrix(hopping):
    t=hopping
    matrix=s.diag(2*t,4+2*t,4+2*t,4+2*t,4+2*t,2*t)
    for doublon in (0,5):
        matrix[doublon,2]=matrix[2,doublon]=-t
        matrix[doublon,3]=matrix[3,doublon]=t
    return matrix


@lru_cache(maxsize=1)
def actual_matrix():
    columns=tuple(source_column(occupation) for occupation in OCCUPATIONS)
    for i,left in enumerate(columns):
        assert all(inner(left,right)==int(i==j) for j,right in enumerate(columns))
        assert all(gamma(left,cell)==left for cell in (0,1))
        assert pair_onsite(pair_onsite(left,0),1)==left
    images=tuple(pair_energy_action(state) for state in columns)
    matrix=s.Matrix(6,6,lambda row,column:inner(columns[row],images[column]))
    assert matrix==expected_matrix(s.pi/4)
    # This is a full source intertwiner check, not just agreement of the
    # projected matrix. Every amplitude outside the six-dimensional image
    # must also vanish in the operator residual.
    for column,image in enumerate(images):
        assert image==source_combination(tuple(matrix[:,column]))
        assert all(gamma(image,cell)==image for cell in (0,1))
    return matrix


def singlet_and_doublon():
    singlet=s.Matrix([0,0,1/s.sqrt(2),-1/s.sqrt(2),0,0])
    doublon=s.Matrix([1/s.sqrt(2),0,0,0,0,1/s.sqrt(2)])
    return singlet,doublon


def transfer_probability(hopping,time):
    """One-per-cell internal singlet to either double-occupation outcome."""
    return hopping**2/(1+hopping**2)*s.sin(2*s.sqrt(1+hopping**2)*time)**2


def certificate():
    actual=actual_matrix()
    t=s.Symbol('t',positive=True);tau=s.Symbol('tau',real=True);lam=s.Symbol('lambda',real=True)
    matrix=expected_matrix(t);singlet,doublon=singlet_and_doublon()
    transform=s.Matrix.hstack(singlet,doublon)
    block=s.Matrix([[4+2*t,-2*t],[-2*t,2*t]])
    assert matrix*transform==transform*block
    characteristic=s.factor((lam*s.eye(6)-matrix).det())
    target=(lam-2*t)*(lam-4-2*t)**3*((lam-2-2*t)**2-4*(1+t**2))
    assert s.expand(characteristic-target)==0
    free=matrix+s.diag(4,0,0,0,0,4)
    free_block=s.Matrix([[4+2*t,-2*t],[-2*t,4+2*t]])
    assert free*transform==transform*free_block
    initial=source_combination(tuple(singlet))
    ambient=PairFockState({key:value for key,value in initial.amplitudes if key[0]==(1,1)},4)
    assert ambient.norm_squared()==1
    preparations=constraint_preparation(ambient)
    assert all(s.simplify(branch.norm_squared()-s.Rational(1,4))==0 for _,branch in preparations)
    radius=2*s.sqrt(1+t**2);center=2+2*t
    reduced=block-center*s.eye(2)
    assert s.simplify(reduced*reduced-radius**2*s.eye(2))==s.zeros(2)
    # Exact 2x2 exponential on the initial singlet.
    amplitudes=s.Matrix([s.cos(radius*tau)-s.I*s.sin(radius*tau)/s.sqrt(1+t**2),
                         s.I*t*s.sin(radius*tau)/s.sqrt(1+t**2)])
    assert all(s.simplify(s.I*s.diff(amplitudes[i],tau)-(reduced*amplitudes)[i])==0 for i in range(2))
    assert amplitudes.subs(tau,0)==s.Matrix([1,0])
    probability=transfer_probability(t,tau)
    assert s.simplify(s.conjugate(amplitudes[1])*amplitudes[1]-probability)==0
    assert s.simplify(sum(s.conjugate(z)*z for z in amplitudes)-1)==0
    # Retain the source link state, wave addresses and every cell occupation
    # outcome at the exact transfer maximum. The global phase is omitted;
    # the same phase multiplies every component and no readout depends on it.
    hopping=s.pi/4;den=s.sqrt(1+hopping**2)
    peak_coefficients=(-s.I*singlet+s.I*hopping*doublon)/den
    peak=source_combination(tuple(peak_coefficients))
    assert s.simplify(peak.norm_squared()-1)==0
    assert all(gamma(peak,cell)==peak for cell in (0,1))
    branches=number_measurement(peak,0)
    maximum=s.pi**2/(16+s.pi**2)
    target_probabilities=(maximum/2,1-maximum,maximum/2)
    assert [record for record,_ in branches]==[0,1,2]
    probabilities=[]
    for (record,branch),target_p in zip(branches,target_probabilities,strict=True):
        value=s.simplify(branch.norm_squared())
        assert s.simplify(value-target_p)==0
        assert all(gamma(branch,cell)==branch for cell in (0,1))
        probabilities.append({'left_wave_number':record,'probability':str(value)})
    assert s.simplify(sum(s.sympify(row['probability']) for row in probabilities)-1)==0
    # Pauli-blocked identical internal labels give no doublon response.
    polarized=s.Matrix([0,1,0,0,0,0])
    assert matrix*polarized==(4+2*t)*polarized
    # Homogeneous fixed-face link support is unchanged by the entire
    # invariant sector and by its complete occupation readout.
    assert all(set(links)<=set((1,1^STAR)) for (links,_),_ in peak.amplitudes)
    return {'source_sector':'two cells, two dark address modes per cell, two waves, positive local Gauss pairs',
            'basis_occupations':[list(x) for x in OCCUPATIONS],
            'full_source_columns_checked':6,'source_intertwiner_has_zero_leakage':True,
            'Hamiltonian_matrix':[[str(actual[i,j]) for j in range(6)] for i in range(6)],
            'hopping_coefficient':'pi/4',
            'onsite_law':'2*(n_a+n_b)-4*n_a*n_b',
            'sector_Hubbard_law':'(2+t)*N-4*(D_left+D_right)-t*sum_flavor(c_left* c_right+c_right* c_left)',
            'six_dimensional_spectrum':[
                {'energy':'2+pi/2-sqrt(16+pi**2)/2','multiplicity':1},
                {'energy':'pi/2','multiplicity':1},
                {'energy':'4+pi/2','multiplicity':3},
                {'energy':'2+pi/2+sqrt(16+pi**2)/2','multiplicity':1}],
            'singlet_doublon_block':[['4+2*t','-2*t'],['-2*t','2*t']],
            'block_energy_splitting':'4*sqrt(1+t**2)',
            'double_occupation_probability':'t**2/(1+t**2)*sin(2*sqrt(1+t**2)*tau)**2',
            'complete_occupation_PVM':['P_left_number_0','P_left_number_1','P_left_number_2'],
            'generic_probabilities':['P_double/2','1-P_double','P_double/2'],
            'first_transfer_maximum_time':'pi/sqrt(16+pi**2)',
            'maximum_double_occupation_probability':str(maximum),
            'actual_source_records_at_transfer_maximum':probabilities,
            'records_conditional_on_positive_Gauss_preparation':True,
            'Gauss_preparation_success_probability':'1/4',
            'retained_Gauss_failure_flags':[{'signs':list(signs),'probability':'1/4'}
                                          for signs,_ in preparations if signs!=(1,1)],
            'peak_original_trial_occupation_probabilities':[
                str(s.simplify(s.sympify(row['probability'])/4)) for row in probabilities],
            'polarized_internal_pair_double_occupation_probability':'0',
            'additive_onsite_countercontrol':{
                'law':'replace A_x by 2*(n_a+n_b), retain hopping and preparation',
                'double_occupation_probability':'sin(2*t*tau)**2',
                'at_t_pi_over_four_and_tau_one_free_double_probability':'1',
                'same_time_interacting_double_probability_upper_bound':str(maximum),
                'minimum_complete_detector_contrast':'16/(16+pi**2)'},
            'constraints_and_total_wave_number_preserved':True,
            'field_words_retained_at_original_cells':True,
            'intercell_regional_topological_charge_transfer':False,
            'internal_pair_singlet_identified_as_physical_spin_singlet':False,
            'contact_interaction_from_declared_positive_clock_not_original_dark_Floquet_phase':True,
            'original_two_cell_dark_onsite_tick_is_identity_on_N2':True,
            'finite_collision_occupation_transfer_executed':True,
            'asymptotic_scattering_matrix_or_cross_section':False,
            'SM_particle_or_GeV_mass_identified':False}
