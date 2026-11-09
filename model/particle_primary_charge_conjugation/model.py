"""Exact full-source particle/hole symmetry and conjugate complete records."""
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_local_wilson_density import density
from model.particle_primary_channel_labels import channels
from model.cubic_particle_channels import group

SOURCE_INPUTS=(
    'model/particle_primary_charge_conjugation/__init__.py',
    'model/particle_primary_charge_conjugation/__main__.py',
    'model/particle_primary_charge_conjugation/model.py',
    'model/particle_primary_charge_conjugation/test_model.py',
    'model/particle_primary_charge_conjugation/PROOF.md',
    'model/particle_primary_charge_conjugation/README.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_channel_labels/PROOF.md',
    'model/particle_primary_pair_channels/PROOF.md',
    'model/particle_primary_charged_readout/PROOF.md',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_spatial_ground_state/PROOF.md',
)


@lru_cache(maxsize=1)
def exact_controls():
    gamma=density.gamma_data();beta=gamma[3];M=gamma[1]
    clean=density.clean;zero=s.zeros(4)
    if M!=M.T or M!=M.conjugate() or M*M!=s.eye(4) or M.det()!=1:
        raise ArithmeticError('actual charge-conjugation spin matrix changed')
    if clean(-M.T*beta.T*M.conjugate())!=beta:
        raise ArithmeticError('onsite term is not particle/hole invariant')
    if clean(M.T*beta.conjugate()*M.conjugate())!=-beta:
        raise ArithmeticError('inversion relation changed')
    for g in gamma[:3]:
        S=(-beta+beta*g)/2
        if clean(-M.T*S.T*M.conjugate()-S.H)!=zero:
            raise ArithmeticError('oriented Wilson term is not conjugated correctly')
    rotations=group.elements()
    for q in rotations:
        lift=channels.spin_lift(q)
        if clean(M.T*lift.conjugate()*M.conjugate()-lift)!=zero:
            raise ArithmeticError('binary cube action is not conjugation compatible')
    native=channels.native_number_control()
    blocks=native['blocks']
    charge_trace=sum(x['Q6']*x['internal_dimension'] for x in blocks)
    if charge_trace!=0:raise ArithmeticError('complete native determinant cancellation failed')
    if any(s.trace(t)!=0 for t in channels.native.native_currents()):
        raise ArithmeticError('full native gauge action has a determinant character')
    # A nonreal, nondiagonal color/weak element checks the literal exterior
    # chart; the proof uses the real-coefficient minor formula for every g.
    z=(s.Integer(3)+4*s.I)/5
    u=s.Matrix([[s.Rational(3,5),4*s.I/5],[4*s.I/5,s.Rational(3,5)]])
    v=s.Matrix([[s.Rational(5,13),s.Rational(12,13)],[-s.Rational(12,13),s.Rational(5,13)]])
    g=s.diag(z**-2*s.diag(u,s.ones(1)),z**3*(u*v))
    R=density.native_representation(g)
    if clean(density.native_representation(g.conjugate())-R.conjugate())!=s.zeros(15):
        raise ArithmeticError('native word conjugation does not match the actual chart')
    if s.simplify(s.expand(R.det()))!=1:raise ArithmeticError('native determinant is not one')
    return {'actual_spin_matrix':[[str(x) for x in M.row(i)] for i in range(4)],
        'spin_matrix_source':'gamma_data()[1]', 'oriented_Wilson_axes_checked':3,
        'binary_cube_lifts_checked':len(rotations), 'native_blocks':blocks,
        'native_integer_charge_trace':charge_trace,'full_native_generators_traceless':True,
        'exterior_coordinate_rule':'R(conjugate(g))=conjugate(R(g))',
        'native_determinant_one':True,'all_native_blocks_required_for_Gauss_phase_cancellation':True}


def bind_member(member):
    contract.verify_member(member);control=exact_controls()
    V=member['source']['side']**3;F=60*V
    if F%4:raise ArithmeticError('full CAR phase convention changed')
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'exact_source_controls':control,'CAR_modes':F,
        'native_occupation_conjugation':[
            {'Q6':b['Q6'],'internal_dimension':b['internal_dimension'],
             'anchor_occupation':2*b['internal_dimension']*V,
             'full_occupation':4*b['internal_dimension']*V,
             'map':'N -> full_occupation - N; Delta_N -> -Delta_N'}
            for b in control['native_blocks']],
        'source_transformation':'unitary full-CAR particle/hole times coordinate-conjugation pullback',
        'source_unitary_square':1,'CAR_parity_preserved':True,
        'Gauss_projector_preserved':True,'product_link_cutoff_preserved':True,
        'original_H_and_energy_zero_preserved':True,'selected_polynomial_density_preserved':True,
        'Hamiltonian_time_not_reversed':True,
        'reference_comparison':'C tensor conjugate(C); not a reference operation in the experiment',
        'opposite_occupation_sector_spectra_equal':True,
        'paired_complete_energy_record_probabilities_equal':True,
        'physical_local_integer_charge_map':'Q6_x -> -Q6_x',
        'proper_binary_cube_rotations_commute':True,
        'spatial_inversion_relation':'C I C^-1 = CAR_parity I',
        'conjugation_is_a_physical_antiparticle_identification':False,
        'Lorentz_spin_or_CPT_proved':False,'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_charge_conjugation.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('charge-conjugation source or consequence changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*channels.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_charge_conjugation_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during replay')
    return value
