"""Source-bound boundary coherence, without an assumed entropy area law."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from model.particle_common_contract import contract, provenance
from model.particle_primary_geometric_cut import model as cut
from model.particle_primary_gravity_premises import model as modular

SOURCE_INPUTS = tuple('model/particle_primary_boundary_information/'+name for name in (
    '__init__.py', '__main__.py', 'model.py', 'test_model.py', 'PROOF.md', 'README.md')) + (
    'model/particle_primary_geometric_cut/PROOF.md',
    'model/particle_primary_gravity_premises/PROOF.md',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_spatial_modular_response/PROOF.md',
    'Lean/QFT/CPRestrictionNet.lean',
)


def entropy(rho):
    values = np.linalg.eigvalsh(rho)
    if min(values) < -1e-12:
        raise ValueError('positive density required')
    return float(-sum(x*np.log(x) for x in values if x > 1e-14))


def singlet_embedding(dimensions, multiplicities):
    """Exact Schur pattern; floating square roots are used only in controls."""
    if not dimensions or len(dimensions) != len(multiplicities):
        raise ValueError('matching nonempty representation data required')
    if any(type(d) is not int or d < 1 for d in dimensions):
        raise ValueError('positive integer dimensions required')
    if any(type(m) is not int or m < 1 or type(n) is not int or n < 1
           for m,n in multiplicities):
        raise ValueError('positive integer multiplicities required')
    nr = sum(d*m for d,(m,n) in zip(dimensions,multiplicities))
    no = sum(d*n for d,(m,n) in zip(dimensions,multiplicities))
    physical = sum(m*n for m,n in multiplicities)
    W = np.zeros((nr*no, physical), complex)
    offset_r = offset_o = column = 0; blocks = []
    for d,(m,n) in zip(dimensions,multiplicities):
        start = column
        for u in range(m):
            for v in range(n):
                for a in range(d):
                    W[(offset_r+u*d+a)*no+offset_o+v*d+a,column] = 1/np.sqrt(d)
                column += 1
        blocks.append((start,column,m,n,d))
        offset_r += m*d; offset_o += n*d
    return W,nr,no,blocks


def information_control(rho, dimensions, multiplicities):
    W,nr,no,blocks = singlet_embedding(dimensions,multiplicities)
    if rho.shape != (W.shape[1],W.shape[1]) or not np.allclose(np.trace(rho),1):
        raise ValueError('normalized density on the actual control support required')
    extended = W@rho@W.conj().T
    tensor = extended.reshape(nr,no,nr,no)
    inside = np.einsum('abcb->ac',tensor)
    outside = np.einsum('abad->bd',tensor)
    pinched = np.zeros_like(rho); weights=[]; edge=0.; conditional=0.
    for begin,end,m,n,d in blocks:
        block = rho[begin:end,begin:end]; p = float(np.trace(block).real)
        pinched[begin:end,begin:end] = block; weights.append(p)
        if p > 0:
            tau = block/p; tt = tau.reshape(m,n,m,n)
            tr = np.einsum('abcb->ac',tt); to = np.einsum('abad->bd',tt)
            conditional += p*(entropy(tr)+entropy(to)-entropy(tau))
            edge += p*np.log(d)
    H = -sum(p*np.log(p) for p in weights if p > 0)
    coherence = entropy(pinched)-entropy(rho)
    algebraic = H+conditional
    mutual = entropy(inside)+entropy(outside)-entropy(rho)
    if abs(mutual-algebraic-2*edge-coherence) > 1e-11:
        raise ArithmeticError('boundary information decomposition failed')
    return {'control_domain':'finite Schur representation control; not a substitute source state',
        'center_entropy':float(H), 'representation_entropy':float(edge),
        'conditional_multiplicity_mutual_information':float(conditional),
        'center_coherence':float(coherence), 'regional_algebra_mutual_information':float(algebraic),
        'ambient_source_mutual_information':float(mutual),
        'identity_residual':float(mutual-algebraic-2*edge-coherence),
        'control_state_substituted_for_primary_source':False}


@lru_cache(maxsize=1)
def native_control():
    dims = [r['internal_dimension'] for r in cut.stress.channels.native_number_control()['blocks']]
    rho = (np.eye(len(dims))+np.ones((len(dims),len(dims))))/(2*len(dims))
    return information_control(rho,dims,[(1,1)]*len(dims))


def bind_member(member):
    contract.verify_member(member)
    geometry = cut.bind_member(member); cut.verify_bound(member,geometry)
    slope = modular.selected_modular_slope(member)
    source = member['source']; a=F(source['spacing'])
    blocks = cut.stress.channels.native_number_control()['blocks']
    native_dimension = sum(r['internal_dimension'] for r in blocks)
    charges = [(sign*r['Q6'],-sign*r['Q6']) for r in blocks for sign in (-1,1)]
    if len(set(charges)) != 10 or native_dimension != 15:
        raise ArithmeticError('the ten distinct physical loop sectors changed')
    matrix_element = -F(source['magnetic_coefficient'])/(2*native_dimension*a)
    floor = geometry['result']['physical_density_floor']; exponent=floor['negative_exponent']
    if floor['base'] != 2 or exponent < 4:
        raise ArithmeticError('source floor must fit the eleven populated sectors')
    derivative_factor = slope['modular_slope_lower']
    coefficient = len(charges)*(derivative_factor*matrix_element)**2
    result = contract.serial({'member_id':member['member_id'], **member['ids'],
        'regional_cut':geometry['result']['induced_cut'],
        'source_density_floor':floor,
        'selected_modular_slope':slope,
        'blank_to_each_crossing_native_character_H_matrix_element':matrix_element,
        'ten_character_boundary_charge_pairs':charges,
        'Hamiltonian_zero_center_commutator_HS_squared_lower':2*len(charges)*matrix_element**2,
        'center_entropy_lower_factored':{'multiplier_of_c_log_2':10*exponent},
        'representation_entropy_lower_factored':{'multiplier_of_c':4,'natural_log_argument':108},
        'center_coherence_lower_factored':{'rational_coefficient':coefficient,
            'base':2,'negative_exponent':2*exponent},
        'classical_center_record_mean_source_work_lower_factored':{
            'rational_coefficient':2*len(charges)*derivative_factor*matrix_element**2,
            'base':2,'negative_exponent':exponent},
        'record_modular_work_balance':'beta_low W <= C_center + D(pinch(rho)||rho) <= beta W',
        'mutual_information_identity':'I_ambient = I_regional_algebra + 2 S_edge + C_center',
        'regional_algebra_identity':'I_regional_algebra = H(p) + sum p_lambda I(tau_lambda)',
        'center_coherence_identity':'C_center = D(rho || pinch(rho)) = S(pinch(rho))-S(rho)',
        'original_selected_state_has_strict_boundary_center_coherence':True,
        'complete_center_record_leaves_regional_entropy_unchanged':True,
        'complete_classical_center_record_disturbs_global_source':True,
        'classical_center_record_has_strictly_positive_mean_source_work':True,
        'source_mutual_information_identified_with_reference_inclusive_information':False,
        'representation_entropy_claimed_distillable_by_local_gauge_invariant_operations':False,
        'actual_regional_probability_weights_numerically_evaluated':False,
        'positive_capacity_uniform_information_bound_proved':False,
        'geometric_area_law_or_Einstein_equation_proved':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_boundary_information.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value != bind_member(member):
        raise ValueError('boundary information or source binding changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*cut.SOURCE_INPUTS,*modular.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_boundary_information_certificate.v1',
        'native_representation_control':native_control(),
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],
        'source_sha256':pins}
    if provenance.definition_hashes(paths) != pins:
        raise ValueError('source changed during boundary information replay')
    return value
