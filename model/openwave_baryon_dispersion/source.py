"""Certified momentum bands of the complete native infinite-chain source."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import numpy as np
import sympy as sy
from scipy.linalg import eigh
from model.openwave_hadron_refinement_flow.analytic import constants, full_bounds
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals, sqrt_upper
from model.openwave_capacity_stress.source import sqrt_interval, enclosure
import model.openwave_hadron_refinement_flow.parameters as parent
ROOT=Path(parent.__file__).resolve().parents[2]


def matrix(cosine,sine,proton=False):
    cosine,sine=F(cosine),F(sine)
    if cosine*cosine+sine*sine!=1:raise ValueError('exact unit-circle phase required')
    a,b,c,p,_,unit=constants(*degree_coefficients());t=p['t'];z=sy.Rational(cosine)-sy.I*sy.Rational(sine)
    A=sy.Matrix([[6*a,-3*t,-3*t,-3*t*z,-3*t*z],
                 [-3*t,c,0,-t/2,-3*t/2],[-3*t,0,b,-3*t/2,-t/2],
                 [-3*t*sy.conjugate(z),-t/2,-3*t/2,c,0],
                 [-3*t*sy.conjugate(z),-3*t/2,-t/2,0,b]])
    assert A==A.conjugate().T
    G=sy.diag(6,1,1,1,1)
    if proton:A-=sy.Rational(p['delta'])*G
    if sine==0:return A,G
    # Realification duplicates every eigenvalue and retains both quadratures.
    R,I=A.applyfunc(sy.re),A.applyfunc(sy.im)
    return R.row_join(-I).col_join(I.row_join(R)),sy.diag(G,G)


@lru_cache(None)
def phase_record(cosine,sine,repair=False,proton=False):
    A,G=matrix(cosine,sine,proton)
    bounds_source=full_bounds(*degree_coefficients())
    if proton:
        bounds_source=dict(bounds_source)
        bounds_source['Q_floor']-=constants(*degree_coefficients())[3]['delta']
        bounds_source['all_P_upper']-=constants(*degree_coefficients())[3]['delta']
    if repair:
        J, norm_upper, _ = repair_matrix(proton)
        if sine: J=sy.diag(J,J)
        A += 2*(1-sy.Rational(cosine))*J
        bounds_source=dict(bounds_source)
        bounds_source['uniform_energy_error']=5*F(1,8)**2/(bounds_source['Q_floor']-bounds_source['all_P_upper']-2*(1-F(cosine))*norm_upper)
        assert bounds_source['uniform_energy_error']>0
    bounds=eigen_intervals(A,G)
    if sine:
        assert all(a==b for a,b in zip(bounds[::2],bounds[1::2]))
        bounds=bounds[::2]
    b=bounds_source;eps=b['uniform_energy_error']
    _,_,_,p,_,unit=constants(*degree_coefficients())
    masses=[((a-eps)*unit,b*unit) for a,b in bounds]
    return {'cosine':str(cosine),'sine':str(sine),'proton':proton,
            'all_five_compact_source_intervals':[[str(x),str(y)] for x,y in bounds],
            'all_five_full_source_GeV_intervals':[[str(x),str(y)] for x,y in masses],
            'all_five_full_source_GeV_display':[enclosure(x) for x in masses],
            'complete_high_complement_source_floor':str(b['Q_floor']),
            'uniform_Schur_error_source':str(eps)}


def run():
    _,_,_,p,_,unit=constants(*degree_coefficients())
    parent_receipt=json.loads((ROOT/'model/openwave_string_positive_cone/results.json').read_text())
    length=F(parent_receipt['common_length_attachment']['a_fm'])
    constants_inputs=json.loads((ROOT/'model/openwave_transverse_gravity/inputs.json').read_text())
    hc=F(constants_inputs['hbar_c_GeV_fm']);momentum_unit=hc/length
    phases=[(F(1),F(0)),(F(3,5),F(4,5)),(F(0),F(1)),(F(-3,5),F(4,5)),(F(-1),F(0))]
    rows=[phase_record(*v) for v in phases]
    rest=tuple(map(F,rows[0]['all_five_full_source_GeV_intervals'][0]));edge=tuple(map(F,rows[-1]['all_five_full_source_GeV_intervals'][0]))
    # pi>3 suffices for a rigorous, very conservative comparison.
    relativistic_lower=sqrt_interval(rest[0]**2+9*momentum_unit**2)[0]
    assert edge[1]<relativistic_lower
    # At any phase, ||H_P(k)-H_P(0)|| <= 2 sqrt(3)t.
    eps=F(rows[0]['uniform_Schur_error_source'])
    uniform_upper=F(rows[0]['all_five_compact_source_intervals'][0][1])*unit+(2*sqrt_upper(3)*p['t'])*unit
    uniform_lower=F(rows[0]['all_five_compact_source_intervals'][0][0])*unit-eps*unit
    return {'schema':'oph.native_baryon_dispersion.v1',
            'verdict':'PASS_ALL_MOMENTUM_BAND_BOUNDS_COMMON_LENGTH_RELATIVISTIC_TEST_FAILS',
            'phase_records':rows,
            'all_real_momenta_ground_GeV_enclosure':enclosure((uniform_lower,uniform_upper)),
            'all_real_momenta_argument':'Native compact Schur theorem; diamagnetic ground lower bound; phase-difference norm upper bound.',
            'zone_edge_energy_gain_GeV':enclosure((edge[0]-rest[1],edge[1]-rest[0])),
            'unit_dictionary':{'electron_GeV':str(unit),'cell_fm':str(length),'hbar_c_GeV_fm':str(hc),'momentum_per_radian_GeV':str(momentum_unit),'length_is_supplied_from_strong_calibrations':True},
            'relativistic_comparison':{'pi_lower_used_exact':'3','required_zone_edge_energy_GeV_lower':str(relativistic_lower),
               'actual_zone_edge_energy_GeV_upper':str(edge[1]),'disjoint':True,
               'separation_GeV_lower':str(relativistic_lower-edge[1])},
            'native_proton_rest_band_shift_GeV_exact':str(-p['delta']*unit),
            'complete_observer_domain':{'low_outcomes':5,'high_complement_retained':True,
               'spectral_readout_work_exact':'0','arbitrary_external_quantum_reference':True,
               'momentum_preparation':'Normalized wave packet or fixed fibre; no uniform infinite-carrier state.'},
            'source_parameters_refitted':False,'physical_relativistic_mass_or_3D_scattering_pole_proved':False,
            'positive_correlated_propagation':repaired()}


@lru_cache(None)
def repair_matrix(proton=False):
    # One positive inverse-rest-operator rule, shared by all native channels.
    parent_receipt=json.loads((ROOT/'model/openwave_string_positive_cone/results.json').read_text())
    length=F(parent_receipt['common_length_attachment']['a_fm'])
    hc=F(json.loads((ROOT/'model/openwave_transverse_gravity/inputs.json').read_text())['hbar_c_GeV_fm'])
    _,_,_,p,_,unit=constants(*degree_coefficients())
    zeta=(hc/length)**2/(2*unit**2)
    A,G=matrix(F(1),F(0),proton)
    bounds=eigen_intervals(A,G)
    J=sy.Rational(zeta)*G*A.inv()*G
    norm=zeta/bounds[0][0]
    assert J==J.T
    return J,norm,zeta

def repaired():
    from model.openwave_hadron_refinement_flow.certified_rest import full
    from model.openwave_hadron_refinement_flow.positivity import inputs
    J,norm,zeta=repair_matrix(); base=full(*degree_coefficients())
    rest=base['rest_certificate']; q=base['full_Q_probability_upper']
    s=base['P_eigenvector_sine_upper']; ep=full_bounds(*degree_coefficients())['uniform_energy_error']
    gap=rest['gap_lower']-ep
    _,_,_,p,_,unit=constants(*degree_coefficients())
    derivative_bound=6*p['t']
    old_curvature_bound=derivative_bound+2*derivative_bound**2/gap
    mean_lo=(1-q)*zeta/rest['energy_interval'][1]-2*norm*s
    mean_hi=zeta/rest['energy_interval'][0]+2*norm*s
    curvature=(2*mean_lo-old_curvature_bound,2*mean_hi+old_curvature_bound)
    assert curvature[0]>0
    momentum_unit_squared=2*zeta*unit**2
    inertial=(momentum_unit_squared/(unit*curvature[1]),momentum_unit_squared/(unit*curvature[0]))
    points=[phase_record(c,s,True) for c,s in ((F(1),F(0)),(F(3,5),F(4,5)),(F(0),F(1)),(F(-3,5),F(4,5)),(F(-1),F(0)))]
    edge=tuple(map(F,points[-1]['all_five_full_source_GeV_intervals'][0]));zero=tuple(map(F,points[0]['all_five_full_source_GeV_intervals'][0]))
    from .native import certificate
    return {'complete_native_correlated_rule':certificate(),'inverse_rest_rule_zeta_exact':str(zeta),
       'native_positive_rule':'sum DeltaC^* [zeta M_compact(0)^(-1)] DeltaC over every native colour-singlet triquark/antitriquark channel of span at most one',
       'calibration_ancestry':'Inherited electron unit and common strong length; inverse mass is the positive native compact rest operator. No per-particle mass target enters the kinetic coefficient.',
       'two_cell_graph_coefficient_exact':'0',
       'k0_complete_N3_operator_change_exact':'0',
       'k0_cancellation_for_higher_number_sectors_claimed':False,
       'full_Fock_positive_added_operator':True,
       'one_flavour_removed_or_selected':False,
       'same_rule_for_all_spin_isospin_and_charge_conjugate_channels':True,
       'executed_repaired_fibre':'native S=I=1/2 neutron sector with all five compact states and full infinite complement',
       'kinetic_operator_norm_upper':str(norm),
       'full_complement_floor':str(full_bounds(*degree_coefficients())['Q_floor']),
       'all_momentum_P_upper':str(full_bounds(*degree_coefficients())['all_P_upper']+4*norm),
       'unmodified_curvature_absolute_bound_source':str(old_curvature_bound),
       'repaired_curvature_source_interval':[str(x)for x in curvature],
       'repaired_inertial_mass_GeV_interval':enclosure(inertial),
       'phase_records':points,'proton_repaired_phase_records':[phase_record(c,s,True,True) for c,s in ((F(1),F(0)),(F(0),F(1)),(F(-1),F(0)))],
       'shared_coefficient_universality':'The same inverse native rest matrix rule acts on proton and neutron; their kinetic coefficients are computed from their own native operator blocks, with no independently selected species constants.',
       'available_ground_bandwidth_GeV_interval':enclosure((edge[0]-zero[1],edge[1]-zero[0])),
       'exact_Lorentz_dispersion_or_scattering_pole':False,
       'other_full_Fock_spectra_transferred':False}
