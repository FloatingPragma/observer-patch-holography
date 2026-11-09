"""Retain beta-vacuum occupation structure in the electric Schur bound.

Every bound applies to the actual finite endpoint on its full physical
carrier. The beta-filled vector enters a local operator inequality; the
selected polynomial preparation is unchanged.
"""
from fractions import Fraction as F

from model.particle_common_contract import contract, provenance
from model.particle_primary_endpoint_bounds import endpoint
from model.particle_primary_spectral_scale import scale
from model.particle_primary_variational_readout import variational
from model.particle_spatial_ground_state import ground
from model.spin_packet_readout.bounds import sqrt_upper


SOURCE_INPUTS = (
    'model/particle_primary_occupation_bound/__init__.py',
    'model/particle_primary_occupation_bound/__main__.py',
    'model/particle_primary_occupation_bound/model.py',
    'model/particle_primary_occupation_bound/test_model.py',
    'model/particle_primary_occupation_bound/PROOF.md',
    'model/particle_primary_endpoint_bounds/PROOF.md',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_primary_variational_readout/PROOF.md',
    'model/particle_common_contract/PROOF.md',
)


def occupation_loss(dimension, onsite, inverse_electric_margin):
    """Certified full-Fock lower shift for m(Dx+Dy)-w Q.

    The balanced-number vacuum has Q mean d, variance one, and its
    off-diagonal image has exactly four beta defects. Other number sectors
    retain their onsite cost. All radicals are rounded upward.
    """
    d=ground._i(dimension,'native dimension',1)
    m=ground._q(onsite,'onsite coefficient')
    w=ground._q(inverse_electric_margin,'inverse electric margin')
    if m<=0 or w<=0:raise ValueError('positive onsite and inverse electric margin required')
    q=F(2*(d+1));vacuum=w*d
    def root(defects):
        b=defects*m-w*q
        return (vacuum-b+sqrt_upper((b+vacuum)**2+4*w*w))/2
    c=root(4)
    route='four-defect inverse-order bound'
    if 2*m-w*q+c<=0:
        c=root(2)
        route='two-defect block bound'
    c=max(c,w*q-m,F())
    if c<vacuum or 2*m-w*q+c<=0:
        raise ArithmeticError('balanced-sector Schur complement is not positive')
    # These exact inequalities certify the radical rounding independently
    # of any floating-point eigensolver.
    coupling_defects=4 if route.startswith('four') else 2
    denominator=coupling_defects*m-w*q+c
    if denominator<=0 or (c-vacuum)*denominator<w*w:
        raise ArithmeticError('occupation Schur inequality failed')
    return {'dimension':d,'Haar_quartic_norm_upper':q,
        'vacuum_Haar_quartic_mean':d,'vacuum_Haar_quartic_variance':'1',
        'off_diagonal_vacuum_image_defects':4,
        'balanced_number_complement_defects_lower':2,
        'other_number_sector_defects_lower':1,
        'inverse_electric_margin':w,'occupation_loss_upper':c,
        'balanced_complement_shifted_margin_lower':2*m-w*q+c,
        'vacuum_resolvent_denominator_lower':denominator,
        'occupation_bound_route':route,
        'active_and_inactive_Dirac_rows_retained':True}


def cluster_data(member, onsite_retained_fraction='0', electric_retained_fraction='0'):
    """Independent onsite/electric allocation, before Gauss restriction."""
    contract.verify_member(member)
    r=ground._q(onsite_retained_fraction,'retained onsite fraction')
    s=ground._q(electric_retained_fraction,'retained electric fraction')
    if max(r,s)>=1:raise ValueError('retained fractions must be less than one')
    source=member['source'];degree=3 if source['side']==2 else 6
    m=F(3,degree)*(1-r);electric=F(source['electric_coefficient'])*(1-s)
    delta=60*(sqrt_upper(m*m+F(1,4))-m)
    margin=F(3,20)*electric-delta
    rows=[]
    if margin>0:
        for row in scale.vacuum_channels()['all_native_blocks']:
            native_margin=electric*F(row['trace_metric_Casimir'])/2-delta
            value=occupation_loss(row['dimension'],m,1/native_margin)
            rows.append({**value,'integer_Q6_charge':row['integer_Q6_charge'],
                'trace_metric_Casimir':F(row['trace_metric_Casimir']),
                'electric_minus_fermion_loss_lower':native_margin})
    return {'onsite_retained_fraction':r,'electric_retained_fraction':s,
        'maximum_interior_vertex_degree':degree,'onsite_per_bond_endpoint':m,
        'cluster_electric_coefficient':electric,
        'pointwise_full_Fock_shift_upper':delta,
        'all_word_nonconstant_link_margin_lower':margin,
        'nonconstant_link_block_positive':margin>0,'native_blocks':rows,
        'per_bond_occupation_loss_upper':sum((row['occupation_loss_upper'] for row in rows),F()) if margin>0 else None,
        'retained_dimensionless_defect_coefficient':3*r,
        'retained_dimensionless_electric_coefficient':F(source['electric_coefficient'])*s/2,
        'retained_H0_fraction_lower':min(r,s)}


def retained_bound(member, onsite_retained_fraction='0', electric_retained_fraction='0'):
    local=cluster_data(member,onsite_retained_fraction,electric_retained_fraction)
    source=member['source'];a=F(source['spacing']);volume=source['side']**3
    loss=lower=second=None
    if local['nonconstant_link_block_positive']:
        loss=(source['interior_hopping_bonds']*local['per_bond_occupation_loss_upper']
              +F(source['magnetic_coefficient'])*source['interior_plaquettes'])
        lower=(F(-90*volume)-loss)/a
        second=lower+6*local['retained_H0_fraction_lower']/a
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'local_bound':local,'dimensionless_anchor_energy':-90*volume,
        'dimensionless_total_loss_upper':loss,'physical_ground_lower':lower,
        'physical_ordered_second_eigenvalue_lower':second,
        'second_eigenvalue_counts_multiplicity':True,
        'full_physical_carrier':True,'primary_endpoint_lambda':'1',
        'selected_state_changed':False,'source_reference_or_clock_changed':False,
        'local_vacuum_is_used_as_selected_preparation':False,
        'endpoint_gap_or_particle_mass_proved':False})


def selected_occupation_bound(member, onsite_retained_fraction='2/5', electric_retained_fraction='0'):
    """Authenticate a ground upper and turn the operator floor into defects."""
    value=retained_bound(member,onsite_retained_fraction,electric_retained_fraction)
    upper=variational.bind_member(member)
    variational.verify_bound(member,upper)
    source=member['source'];a=F(source['spacing']);volume=source['side']**3
    energy_upper=F(upper['result']['ground_energy_upper'])
    excess=F(member['state']['target_ground_excess'])
    coefficient=F(value['local_bound']['retained_dimensionless_defect_coefficient'])
    defects=weight=None
    if value['physical_ground_lower'] is not None and coefficient>0:
        numerator=a*(energy_upper+excess)-F(value['dimensionless_anchor_energy'])
        numerator+=F(value['dimensionless_total_loss_upper'])
        if numerator<0:raise ArithmeticError('selected-state energy bounds are inconsistent')
        # D counts at most all 60 modes per cell.
        defects=min(F(60*volume),numerator/coefficient)
        weight=max(F(),1-defects/(2*volume))
    floor=retained_bound(member)
    inherited=endpoint.bind_member(member)
    return contract.serial({**value,'schema':'oph.primary_occupation_result.v1',
        'physical_ground_upper':energy_upper,
        'upper_receipt_sha256':upper['result_sha256'],
        'selected_ground_energy_excess_upper':excess,
        'selected_mean_source_energy_upper':energy_upper+excess,
        'defect_count_upper':defects,
        'matched_q6_pair_count':2*volume,
        'average_matched_q6_pair_weight_lower':weight,
        'average_matched_q6_pair_weight_simple_lower':F(11,16) if weight is not None and weight>=F(11,16) else F(),
        'weight_is_uniform_address_average':True,
        'weight_is_conditioned_on_success':False,
        'each_individual_site_weight_claimed':False,
        'zero_retention_physical_ground_lower':floor['physical_ground_lower'],
        'inherited_endpoint_physical_ground_lower':inherited['result']['physical_ground_lower'],
        'zero_retention_ground_enclosure_width':energy_upper-F(floor['physical_ground_lower']),
        'physical_mass_table_GeV':[]})


def bind_member(member, onsite_retained_fraction='2/5', electric_retained_fraction='0'):
    result=selected_occupation_bound(member,onsite_retained_fraction,electric_retained_fraction)
    return {'schema':'oph.primary_occupation_bound.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    try:
        local=value['result']['local_bound']
        expected=bind_member(member,local['onsite_retained_fraction'],local['electric_retained_fraction'])
    except (KeyError,TypeError) as error:
        raise ValueError('malformed occupation-bound receipt') from error
    if value!=expected:raise ValueError('occupation bound or its actual source changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*endpoint.SOURCE_INPUTS,*variational.SOURCE_INPUTS))
    if any(not (contract.ROOT/name).is_file() for name in paths):
        raise ValueError('an occupation-bound evidence path is absent')
    pins=provenance.definition_hashes(paths)
    zero=contract.member(0);one=contract.member(1)
    result={'schema':'oph.primary_occupation_certificate.v1',
        'selected_member_bound':bind_member(zero),
        'ground_floors':[retained_bound(m) for m in (zero,one)],
        'onsite_retention_trials':[selected_occupation_bound(zero,r) for r in ('1/5','3/10','2/5','1/2')],
        'equal_H0_retention_trials':[retained_bound(zero,r,r) for r in ('0','1/10','1/5','1/4')],
        'larger_member_default_onsite_admission':retained_bound(one,'2/5','0'),
        'uniform_refinement_defect_or_weight_bound_proved':False,
        'endpoint_gap_or_mass_proved':False,'physical_mass_table_GeV':[],
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:
        raise ValueError('occupation-bound inputs changed during replay')
    return result
