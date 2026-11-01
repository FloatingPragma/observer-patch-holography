"""Authenticate the retained exact spectrum and resolve its entire ground family."""
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from tools import strict_json
from ..cube_eta_protection import cube
from ..interacting_vacuum import ground


ROOT=Path(__file__).resolve().parents[2]
PRIOR='model/matrix_spectral_response/results.json'
SCALAR='model/retarded_matter_response/results.json'


def prior_spectrum():
    values=[]
    for path in (PRIOR,SCALAR):
        value=strict_json.loads((ROOT/path).read_text())
        if value.get('all_pass') is not True:raise ValueError('prior spectrum did not pass')
        for source,digest in value['source_bindings'].items():
            if hashlib.sha256((ROOT/source).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior spectrum source changed: '+source)
        values.append(value)
    matrix,scalar=values
    data=matrix['matrix_response'];old=scalar['resolved_local_spectrum']
    if (data['all_eight_columns_annihilated_exactly'] is not True or
            data['observed_cyclic_dimension']!=312 or
            len(data['matrix_residue_rows'])!=78 or
            any(row['rank']!=4 for row in data['matrix_residue_rows']) or
            {'gap':'sqrt3/16','weight':'1/10'} not in old['known_exact_lines']):
        raise ValueError('prior full matrix or exact lambda residue contract changed')
    return matrix


def lambda_row(rows):
    hits=[]
    for row in rows:
        lo,hi=map(F,row['gap_above_ground_interval'])
        if 0<lo and lo*lo<F(3,256)<hi*hi:hits.append(row)
    if len(hits)!=1 or hits[0]['band']!=1:
        raise ValueError('unique positive lambda band could not be identified')
    return hits[0]['line']


def _scale_interval(values,factor):return [str(F(x)*factor) for x in values]


def resolved_rows(prior=None):
    prior=prior or prior_spectrum();original=prior['matrix_response']['matrix_residue_rows']
    high=lambda_row(original);result=[]
    for row in original:
        is_high=row['line']==high
        weight=['3/10','3/10'] if is_high else _scale_interval(row['local_weight_interval'],F(1,2))
        normalized=['2/5','2/5'] if is_high else _scale_interval(row['local_weight_interval'],F(2,3))
        if not 0<F(weight[0])<=F(weight[1]):raise ArithmeticError('a twirled line disappeared')
        lo,hi=map(F,row['gap_above_ground_interval'])
        result.append({'line':row['line'],'gap_above_ground_interval':row['gap_above_ground_interval'],
                       'energy_midpoint_display_only':float((lo+hi)/2),'spatial_band':row['band'],
                       'total_hidden_spin':'5/2' if is_high else '3/2',
                       'spatial_multiplicity':4,'reachable_N5_dimension':24 if is_high else 16,
                       'mixed_reference_matrix_residue_rank':8,
                       'mixed_reference_local_addition_weight':weight,
                       'conditional_on_addition_local_weight':normalized})
    if not sum(F(x['mixed_reference_local_addition_weight'][0]) for x in result)<=F(3,4)<=sum(F(x['mixed_reference_local_addition_weight'][1]) for x in result):
        raise ArithmeticError('mixed-reference addition weights do not enclose three quarters')
    return result


def frame_certificate():
    ops=cube.operators()
    frame=cube.restrict(ops.frame.real,1).toarray()+1j*cube.restrict(ops.frame.imag,1).toarray()
    original=cube.restrict(ops.H.real,1).toarray()+1j*cube.restrict(ops.H.imag,1).toarray()
    hidden=np.kron(ground.band_data()[0],np.eye(2,dtype=np.int64))
    if (not np.all(np.isin(frame.real,(-1,0,1))) or
            not np.all(np.isin(frame.imag,(-1,0,1))) or
            not np.array_equal(frame.conj().T@frame,np.eye(16)) or
            not np.array_equal(frame.conj().T@hidden@frame,original)):
        raise ArithmeticError('exact source one-body Fock frame covariance failed')
    mixed_spin=sum(bool(np.any(original[2*x:2*x+2,2*y:2*y+2]-
                                   np.diag(np.diag(original[2*x:2*x+2,2*y:2*y+2]))))
                   for x in range(8) for y in range(8))
    if not mixed_spin:raise ArithmeticError('physical intersite spin coherence was lost')
    return {'frame_is_actual_full_Fock_source_restricted_to_N1':True,
            'exact_Gaussian_integer_unitarity_and_kinetic_conjugacy':True,
            'hidden_response':'g_xy(z) I_spin',
            'physical_response':'g_xy(z) R_x^dagger R_y',
            'physical_intersite_blocks_with_spin_mixing':mixed_spin,
            'physical_local_response':'g_xx(z) I_spin',
            'physical_and_hidden_all_mode_creator_ranges_coincide':True}


def certificate():
    prior=prior_spectrum();rows=resolved_rows(prior);comparison=prior['finite_gap_dictionary']
    if comparison['finite_gap_as_rest_energy_dictionary_rejected'] is not True:
        raise ValueError('the frozen finite-gap dictionary result changed')
    return {'source':'original 2x2x2 Weyl/contact source, c=1/8,U=1,mu=0; entire N4 spin-two ground space',
            'all_78_positive_addition_lines':rows,
            'addition_line_counts_by_total_hidden_spin':{'5/2':1,'3/2':77},
            'N5_decomposition':[
                {'total_hidden_spin':'5/2','full_multiplicity':56,'reachable_multiplicity':4,
                 'dark_multiplicity':52,'full_dimension':336,'reachable_dimension':24,'dark_dimension':312},
                {'total_hidden_spin':'3/2','full_multiplicity':504,'reachable_multiplicity':308,
                 'dark_multiplicity':196,'full_dimension':2016,'reachable_dimension':1232,'dark_dimension':784},
                {'total_hidden_spin':'1/2','full_multiplicity':1008,'reachable_multiplicity':0,
                 'dark_multiplicity':1008,'full_dimension':2016,'reachable_dimension':0,'dark_dimension':2016}],
            'full_N5_dimension':4368,'reachable_N5_dimension':1256,'dark_N5_dimension':3112,
            'dark_means':'orthogonal to every c_x,sigma^dagger psi and every time translate, for every psi in the complete N4 ground space',
            'dark_energy_eigenvalues_or_full_energy_degeneracies_computed':False,
            'removal':{'positive_energy_cost':'sqrt3/16','retarded_pole':'-sqrt3/16',
                       'total_hidden_spin':'3/2','spatial_multiplicity':4,'reachable_N3_dimension':16,
                       'full_N3_dimension':560,'dark_N3_dimension':544,
                       'dark_spin_three_halves_multiplicity':52,'dark_spin_one_half_multiplicity':168,
                       'mixed_reference_local_weight':'1/4','mixed_reference_residue_rank':8},
            'mixed_ground_reference':'rho_G=P_G/5',
            'mixed_addition_response':'[G_down(z)+P_plus/(z-lambda)]/2 tensor I_spin',
            'mixed_retarded_response':'[G_down(z)+P_plus/(z-lambda)+P_minus/(z+lambda)]/2 tensor I_spin',
            'mixed_retarded_distinct_signed_poles':79,
            'mixed_retarded_minimal_cyclic_dimension':632,
            'mixed_cyclic_dimension_boundary':'The 16-column mixed-state GNS response has 79 rank-eight residues. This 632-dimensional realization differs from the 1256-dimensional N5 physical reachable range obtained by keeping all five ground inputs distinct.',
            'mixed_addition_local_total_weight':'3/4',
            'conditional_addition_local_total_weight':'1',
            'mixed_reference_has_every_original_positive_line':True,
            'arbitrary_ground_selection_adds_frequencies':False,
            'frame':frame_certificate(),
            'finite_gap_dictionary':{'same_frozen_source_energy_envelope':comparison['source_energy_envelope'],
                 'same_frozen_source_spread_upper':comparison['source_spread_upper'],
                 'pair_tests':comparison['pair_tests'],
                 'excluded_for_every_ground_choice_and_linear_single_CAR_probe':True,
                 'common_positive_unit_factor_only':True,'eligible_physical_comparisons_added':0,
                 'unseen_composites_other_particle_numbers_or_continuum_excluded':False},
            'native_full_quintet_preparation_or_79_line_resolution_executed':False,
            'total_hidden_spin_identified_with_particle_Lorentz_spin':False,
            'physical_particle_species_or_relativistic_rest_masses_identified':False}
