"""A declared orthogonal Wilson source on the fixed native code and Gauss frame."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import argparse
import hashlib
import json
from pathlib import Path
import sys
import sympy as s
from . import algebra as a
from . import vacuum
from model.native_wilson_maxwell import field,bounds
from model.native_wilson_maxwell_ground import ground
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import field_data
from ..provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=a.ROOT
PARENT='model/native_wilson_maxwell/results.json'
PARENT_SHA='eba4963148b62f6841a038119ed5b9baf4a6f142009915b31e49d72dacd12b72'
CONE=HERE.parent/'native_wilson_cone_audit'
CONE_SHA='5a24c444fec3ca0a03f59dde8fd52614669b9d9c730dcd4887563727c4cfa62b'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(None)
def parents():
    if digest(ROOT/PARENT)!=PARENT_SHA or digest(CONE/'research_reference.json')!=CONE_SHA:
        raise ValueError('fixed native source or cone calculation differs')
    value=json.loads((ROOT/PARENT).read_text())
    for path,sha in value['input_sha256'].items():
        if digest(ROOT/path)!=sha:raise ValueError('native source pin differs: '+path)
    cone=json.loads((CONE/'research_reference.json').read_text())
    for path,sha in cone['owned_source_sha256'].items():
        if digest(CONE/'retained'/(path+'.raw'))!=sha:raise ValueError('cone source pin differs: '+path)
    return value


def coulomb_control():
    """Convex quadratic maxima on the complete code and native charge boxes."""
    G=field.cube()['G'];R=s.eye(8)-s.eye(8)[:,0]*s.ones(1,8);C=R.T*G*R
    assert C==C.T and C.is_positive_semidefinite
    corners=[]
    for q in product((-2,2),repeat=8):
        v=s.Matrix(q);value=(v.T*C*v)[0]/2
        corners.append({'charges':list(q),'Coulomb_over_m_g_squared':str(value)})
    maximum=max(F(row['Coulomb_over_m_g_squared'])for row in corners)
    assert maximum==F(116,3)
    # On total charge zero, Rq=q; ||G||=1/2 and ||q||²<=32.
    L=field.cube()['D'].T*field.cube()['D']
    assert (L*(L-2*s.eye(8))*(L-4*s.eye(8))*(L-6*s.eye(8)))==s.zeros(8)
    neutral=F(8);native=4*maximum
    assert native==F(464,3)
    return {'exact_all_256_code_box_corners':corners,
        'all_code_charge_Coulomb_norm_over_m_g_squared':'116/3',
        'neutral_code_Coulomb_norm_over_m_g_squared':'8',
        'all_native_physical_charge_Coulomb_norm_over_m_g_squared':'464/3',
        'full_code_charge_box':[-2,2],'full_native_charge_box':[-4,4],
        'static_compensator':'R=-sum_x rho_x; its energy4|R| is additional and retained',
        'finite_root_compensation_used':True}


def edge_phase_control():
    """Full four-mode edge lift after the literal rank-two bond SVD."""
    from model.native_block_dirac.algebra import annihilator
    c=[annihilator(4,j)for j in range(4)]
    z=s.Rational(3,5)+4*s.I/5
    pieces=[]
    for j in range(2):
        K=c[j+2].H*c[j]
        D=(z-1)*K+(s.conjugate(z)-1)*K.H
        # One CAR pair has spectrum0,+/-|z-1|; tensor spectators retained.
        assert D**3==s.Rational(4,5)*D
        pieces.append(D)
    assert (pieces[0]*pieces[1]-pieces[1]*pieces[0]).applyfunc(s.expand)==s.zeros(16)
    D=sum(pieces,s.zeros(16))
    assert (D*(D*D-s.Rational(4,5)*s.eye(16))*(D*D-s.Rational(16,5)*s.eye(16))).applyfunc(s.expand)==s.zeros(16)
    return {'actual_rank_two_partial_isometry_verified':True,
        'complete_four_mode_Fock_dimension':16,
        'rational_phase':'3/5+4i/5','exact_difference_spectrum':'0, +/-2/sqrt(5), +/-4/sqrt(5)',
        'general_edge_full_Fock_difference_norm':'2m|exp(iq)-1| <=2m|q|',
        'twelve_edge_Peierls_difference_norm_over_m':'18g',
        'argument':'diagonalize the retained finite Q_e for each edge bound; triangle sum does not require commuting different Q_e'}


def interaction(vac,parent):
    sigma=F(vac['new_pair_one_photon_current_norm_over_m_lower'])
    g=tau=F(1,64)
    # Neutral sector only. The older conservative Taylor bound is retained;
    # the whole-charge ground proof uses its separate116/3 coefficient.
    nu=18*g+F(59,4)*g*g
    leading=g*tau*sigma*(1-F(11,2)*tau)
    odd=F(27,16)*tau*g**3;second=tau*tau*nu*nu/2
    amplitude=a.rounded(leading-odd-second,False)
    weight=F(parent['actual_preparation']['global_code_probability_lower'])
    connected=weight*(amplitude**2-(tau*nu)**4)
    assert amplitude>0 and connected>0
    values={'g':g,'dimensionless_mT':tau,'neutral_perturbation_norm_over_m_upper':nu,
        'leading_target_amplitude_lower':leading,'odd_first_insertion_remainder_upper':odd,
        'full_second_Duhamel_remainder_upper':second,'ideal_target_amplitude_lower':amplitude,
        'actual_code_probability_lower':weight,
        'actual_original_attempt_target_probability_lower':a.rounded(weight*amplitude**2,False),
        'actual_original_attempt_target_complement_probability_upper':a.rounded(1-weight*amplitude**2),
        'actual_original_attempt_connected_covariance_lower':a.rounded(connected,False)}
    return {'exact':{k:str(v)for k,v in values.items()},
        'diagnostic_binary64':{k:float(v)for k,v in values.items()},
        'target':[1,1,1],'free_target_probability':'0',
        'input':'actual deterministic native trial, followed by the declared code-only sea rotation and the unchanged Gauss lift',
        'new_sea_count_projectors_differ_from_parent_count_projectors':True,
        'source_state_is_actual_coupled_ground':False,
        'interpretation':'finite prepared-free-vacuum dressing under the declared changed Hamiltonian',
        'old_event_probability_or_vacuum_reused':False}


def source_gap(parent,alignment):
    native=ground.native_complement(parent)
    delta_l,delta_u=map(F,native['actual_neutral_triplet_gap_interval'])
    ml,mu=map(F,parent['clock_and_capacity']['source_selected_excitation_interval'])
    assert 0<ml<mu<1 and 0<delta_l<delta_u
    shift_l,shift_u=map(F,alignment['signed_code_vacuum_shift_over_m'])
    nu=F(alignment['all_code_charge_perturbation_norm_over_m_upper'])
    c=F(alignment['robust_declared_positive_complement_coefficient'])
    code_gap=ml*F(alignment['robust_code_excitation_gap_over_m_lower'])
    complement_gap=delta_l+ml*(c-shift_u-nu)
    full_gap=min(code_gap,complement_gap)
    assert full_gap>ml/2
    _,polynomial,_=field_data()
    roots=polynomial.intervals(eps=s.Rational(1,10**100))
    assert len(roots)==6 and all(mult==1 for _,mult in roots)
    elo,ehi=map(F,roots[0][0]);v=F(360,2**34)
    # Eg lies in[8 e_l-v,8 e_u+v]. No changed sea energy is subtracted.
    eg=(8*elo-v,8*ehi+v)
    raw_lower=32*eg[0]+ml*(3+a.root(6,False)+shift_l-nu)
    raw_upper=32*eg[1]+mu*(3+a.root(6)+shift_u+nu)
    assert raw_upper<0
    return {'source_selected_mass_interval':[str(ml),str(mu)],
        'native_physical_complement':native,
        'raw_energy_baseline':'32 Eg + m(3+sqrt(6)); all native and photon zero points retained',
        'full_code_ground_energy_minus_raw_baseline_over_m':a.encode(shift_l-nu,shift_u+nu),
        'full_physical_ground_raw_energy_interval':a.encode(raw_lower,raw_upper),
        'full_physical_ground_gap_source_units_lower':str(a.rounded(full_gap,False)),
        'unique_active_physical_ground_lies_in_complete_native_code':True,
        'arbitrary_zero_energy_preparation_reference_spectators_excluded_from_uniqueness_only':True,
        'native_Gauss_violating_kinematic_states_assigned_physical_ground_interpretation':False,
        'outer_Gauss_orthogonal_zero_extension_under_cuts_physical_ground':False,
        'coupled_ground_vector_or_ground_preparation_numerically_constructed':False,
        'diagnostic_binary64':{'mass_interval':[float(ml),float(mu)],
            'raw_ground_interval':[float(raw_lower),float(raw_upper)],
            'full_gap_source_units_lower':float(full_gap),
            'native_complement_gap':[float(delta_l),float(delta_u)]}}


def resources(parent,alignment):
    upper=F(parent['clock_and_capacity']['source_selected_excitation_interval'][1])
    c=F(alignment['robust_declared_positive_complement_coefficient']);g=F(1,64)
    # Code replacement<139m; use160m. Completed photons<13m;
    # full native root-compensated Coulomb<=464/3 mg²; retain4|R|<=128.
    raw=32*(600+F(360,2**34))+upper*(173+18*g+F(464,3)*g*g+8*c)+128
    return {'full_raw_H_ext_norm_upper':str(raw),'native_constituents':32,
        'preparation_duration_source_units':'1','autonomous_wait':'1/(64m)',
        'readout_duration_source_units':'1','complete_elapsed_source_time':'2+1/(64m)',
        'preparation_unitary':'sea rotation composed with all original deterministic trial preparations; full Gauss lift',
        'preparation_control':'i Log(S)-H_ext, including explicit full-H cancellation',
        'readout_control':'i Log(exp(-i H_ext) U_count)-H_ext; retain final source evolution',
        'each_control_norm_upper':str(raw+F(22,7)),
        'net_bare_source_work_absolute_upper':str(2*raw),
        'gross_two_pulse_switching_work_upper':str(4*(raw+F(22,7))),
        'count_pointer_dimension':1161,'physical_count_labels':1160,
        'full_native_reference_field_pointer_capacity_nats_formula':parent['clock_and_capacity']['full_native_source_reference_trial_capacity_nats_formula'],
        'full_native_reference_field_pointer_capacity_nats_interval':parent['clock_and_capacity']['full_native_source_reference_trial_capacity_nats_interval'],
        'all_native_complement_photon_and_Gauss_failure_outcomes_retained':True,
        'same_elapsed_clock_on_all_failures':True,
        'controller_coefficient_workspace_in_capacity':False,
        'efficient_native_gate_compilation_proved':False}


def certificate():
    parent=parents();alignment=a.energy_alignment();vac=vacuum.calculate()
    from ..native_wilson_cone_audit import cone
    algebra=a.covariance();native_control=ground.local_completion_control()
    value={'schema':'oph.research.native_orthogonal_wilson_bridge.v1',
        'source_parent_sha256':PARENT_SHA,'cone_parent_sha256':CONE_SHA,
        'declared_changed_terms':{
            'onsite_regulator':'(degree(x)/(2a)) (Gamma-beta) in the fixed canonical field',
            'each_directed_bond_regulator':'-(Gamma-beta)/(2a), with the same Peierls phase and longitudinal orientation',
            'native_mass_selector_or_selected_native_states_changed':False,
            'canonical_field_charge_global_code_or_Gauss_map_changed':False,
            'old_code_scalar_minus_Esea_old_retained':True,'new_code_vacuum_subtraction_added':False,
            'positive_native_complement_completion':'(21m/2) sum_x(I-P_x)',
            'all_changed_terms_are_declared_architecture':True},
        'common_free_cone':cone.algebra(),'finite_static_background_control':cone.finite_background_check(),
        'fixed_charge_Gauss_and_code':algebra,
        'full_native_code_pattern_controls':native_control,'finite_full_Fock_phase_control':edge_phase_control(),
        'full_charge_Coulomb_control':coulomb_control(),'energy_alignment':alignment,
        'exact_FW_analysis':a.fw_identity(),'exact_sea_rotation_algebra':a.rotation_algebra(),'fresh_vacuum':vac,
        'full_physical_ground':source_gap(parent,alignment),
        'complete_prepared_interaction':interaction(vac,parent),'resources':resources(parent,alignment),
        'scope':{'common_quadratic_free_cone':True,'fresh_full_physical_finite_source_gap':True,
            'original_beta_Wilson_source_result_reinterpreted':False,
            'interacting_relativistic_or_spatial_continuum_limit':False,
            'all_state_native_tensor_locality':False,'measured_SM_or_GeV_mass':False,
            'finite_photon_CCR_boundary_removed':False}}
    seeds=set(parent['input_sha256'])|{PARENT,ground.GEOMETRY,
        'model/native_abelian_field_join/results.json','model/native_wilson_maxwell_ground/ground.py',
        'model/native_wilson_maxwell_ground/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/charged_ground_exchange/resolvent.py',
        'model/native_wilson_maxwell_refinement/bounds.py'}
    value['input_sha256']={p:digest(ROOT/p)for p in sorted(import_closure(seeds))}
    value['research_parent_source_sha256']={p.name.removesuffix('.raw'):digest(p)for p in sorted((CONE/'retained').iterdir())}
    value['research_parent_source_sha256']['results.json']=digest(CONE/'research_reference.json')
    value['owned_source_sha256']={p.name:digest(p)for p in sorted(HERE.iterdir())if p.suffix in ('.py','.md')}
    return json.loads(json.dumps(value,default=str))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    value=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    elif value!=json.loads(path.read_text()):raise ValueError('orthogonal source certificate differs')
    print('PASS: fixed native code/Gauss, common free cone, signed raw sea shift, full-charge gap and complete records')
