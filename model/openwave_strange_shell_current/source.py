"""Certified bulk Lambda-pole current output with a separate open shell."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import sympy as sy
from model.openwave_baryon_poles.source import trial_moments,residue_bounds
from model.openwave_cluster_bandwidth import lambda_source as lam
from model.openwave_cluster_bandwidth.source import shell,angle
from model.openwave_capacity_stress.source import enclosure
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals,sqrt_upper
from .native import contraction,frame
HERE=Path(__file__).resolve().parent


def trial(record):
    assert record['denominator']==str(2**40)
    row=record['real_integer_coefficients']
    assert type(row)is list and len(row)==7
    assert all(type(x)is str and str(int(x))==x and abs(int(x))<=2**40 for x in row)
    return sy.Matrix([sy.Rational(x,record['denominator'])for x in row])


@lru_cache(None)
def lambda_pole():
    blocks,G,d=lam.compact();A=sum(blocks.values(),sy.zeros(7));levels=eigen_intervals(A,G)
    parent=lam.ground();C=F(parent['complete_relative_complement_floor_source']);b2=(6*d['p']['t'])**2
    # Uniform seven-band Schur budget also certifies ground isolation.
    eps=b2/(C-levels[-1][1]);assert C>levels[-1][1]
    record=json.loads((HERE/'lambda_trial.json').read_text());v=trial(record)
    n,mean,var=trial_moments(A,G,v)
    cert=residue_bounds(levels,C,b2,eps,mean,var)
    cert['seven_field_pole_residue_trace_interval']=cert.pop('five_field_pole_residue_trace_interval')
    cert['original_input_pole_event_probability_lower']=str(F(cert['prepared_trial_pole_weight_lower'])/7)
    return {'full_ground_GeV_interval':parent['full_Lambda_ground_GeV_interval'],
        'full_ground_GeV_display':parent['full_Lambda_ground_GeV_display'],
        'all_seven_compact_source_intervals':[list(map(str,x))for x in levels],
        'complete_relative_Q_floor_source':str(C),'uniform_seven_band_Schur_error_source':str(eps),
        'actual_trial_coordinates':record,'trial_Gram_norm_squared':str(n),
        'trial_energy_source':str(mean),'trial_variance_source':str(var),'certificate':cert}


@lru_cache(None)
def refined_shell():
    a,b=map(F,lambda_pole()['full_ground_GeV_interval'])
    # Rational endpoint witnesses are independently replayed; no root-finder
    # trajectory is a proof premise or part of the custody receipt.
    lo,hi=F(7651,40960),F(11477,61440)
    left,right=shell(lo),shell(hi)
    assert F(left['complete_pair_energy_GeV_interval'][1])<a
    assert F(right['complete_pair_energy_GeV_interval'][0])>b
    root=Path(__file__).resolve().parents[2]
    length=json.loads((root/'model/openwave_string_positive_cone/results.json').read_text())['common_length_attachment']['a_fm']
    hc=json.loads((root/'model/openwave_transverse_gravity/inputs.json').read_text())['hbar_c_GeV_fm']
    momenta=(angle(lo)[0]*F(hc)/F(length),angle(hi)[1]*F(hc)/F(length))
    return {'half_angle_tangent_interval':[str(lo),str(hi)],'opposite_momentum_GeV_interval':enclosure(momenta),
        'below_pair_GeV_interval':left['complete_pair_energy_GeV_interval'],
        'above_pair_GeV_interval':right['complete_pair_energy_GeV_interval'],
        'below_margin_GeV':str(a-F(left['complete_pair_energy_GeV_interval'][1])),
        'above_margin_GeV':str(F(right['complete_pair_energy_GeV_interval'][0])-b),
        'complete_interacting_two_cluster_scattering_state_constructed':False}


def run():
    native=contraction();_,P=frame();G=P.T*P;M=native['seven_Gram_without_Cabibbo'];cb=native['Cabibbo_squared']
    lower=M-6*G;upper=7*G-M
    # Exact block identities; no numerical positivity test defines this claim.
    assert lower==sy.diag(0,0,4,4,4,4,0)+sy.Matrix(7,7,lambda i,j:2 if (i,j)in((2,4),(4,2),(3,5),(5,3))else 0)
    assert upper==sy.diag(2,2,0,0,0,0,2)
    assert all(x>=0 for x in lower.eigenvals()) and all(x>=0 for x in upper.eigenvals())
    pole=lambda_pole();q=F(pole['certificate']['full_Q_probability_upper'])
    low=6*cb*(1-q);high=7*cb
    assert low>0
    z=F(pole['certificate']['prepared_trial_pole_weight_lower'])
    return {'schema':'oph.native_strange_bulk_pole_current.v1',
        'verdict':'PASS_NONZERO_NATIVE_BULK_LAMBDA_POLE_CURRENT_WITH_OPEN_KINEMATIC_SHELL',
        'Lambda_source_pole':pole,
        'exact_native_vertex':{'compact_uds_shapes':56,'Lambda_spin_half_isospin_zero_rank':7,
            'raw_N5_output_occupations':native['raw_N5_output_count'],'literal_CAR_words':native['literal_CAR_word_count'],
            'Cabibbo_sine_inherited':'0.22431','Cabibbo_factor_squared_exact':str(cb),
            'Gram_matrix_without_Cabibbo':[[str(x)for x in row]for row in M.tolist()],
            'source_Gram_matrix':[[str(x)for x in row]for row in G.tolist()],
            'Gram_order_exact':'6 G <= W_P* W_P/[sin(theta)^2 cos(theta)^2] <= 7 G',
            'all_56_shape_vertex_norm_squared_upper_exact':str(native['Gram_norm_squared_upper_without_Cabibbo']*cb),
            'all_compact_native_N5_outputs_retained':True,'fitted_mass_vectors_used_for_current':False,
            'full_initial_uds_sector_vertex_norm_upper':str(36*sqrt_upper(cb)),
            'full_initial_uds_sector_vertex_norm_argument':'Each term annihilates the unique strange quark. Distinct creation sites have orthogonal antiquark-position outputs; at each site at most36 norm-one CAR words act.',
            'output_flavour_counts':{'u':2,'d':2,'anti_u':1},'net_charge_change_exact':'0','net_baryon_change_exact':'0'},
        'full_pole_to_complete_compact_N5_response':{
            'amplitude_norm_squared_interval':[str(low),str(high)],
            'amplitude_norm_squared_display':enclosure((low,high)),
            'nonzero_for_every_compact_Lambda_state':True,
            'complete_full_pole_Q_contribution_discarded':False,
            'transfer':'The vertex preserves occupied spatial span, so compact and noncompact output spaces are orthogonal. The exact compact Gram lower bound and full-pole compact weight suffice; the Q output is retained separately.'},
        'complete_component_probe':{
            'success_operator':'K_s=W/36 on the complete one-u one-d one-s, no-antiquark source; zero on its full Fock complement.',
            'failure_operator':'K_f=sqrt(I-K_s* K_s), including identity on the complement.',
            'terminal_outcomes':'Complete compact N5 output, every noncompact N5 output, and the full probe failure outcome.',
            'source_pole_compact_N5_probability_interval':[str(low/1296),str(high/1296)],
            'source_pole_noncompact_N5_probability_interval':['0',str(cb*q)],
            'source_pole_probe_failure_probability_interval':[str(1-high/1296-cb*q),str(1-low/1296)],
            'original_seven_coordinate_input_trial_pole_current_compact_event_probability_lower':str(z*low/(7*1296)),
            'preparation_success_exact':'1/7','preparation_failure_exact':'6/7',
            'nonselective_preparation_work_exact':'0 on the declared P/7 input',
            'probe_work':'Tr[H(K_s rho K_s*+K_f rho K_f* - rho)] is retained; neither the probe nor its spatial-span tag is claimed energy conserving.',
            'zero_probe_work_asserted':False,
            'probabilities_are_fixed_total_momentum_fibre_quantities':True,
            'normalizable_wave_packet_transfer_executed':False,
            'operator_is_channel_resolved_native_component':True,
            'full_native_weak_Hamiltonian_replaced_by_component':False,
            'complete_autonomous_weak_dynamics_executed':False},
        'same_source_energy_shell':refined_shell(),
        'new_mass_or_current_calibration':False,
        'inherited_calibrations':'Electron GeV anchor, common strong/axial inputs, strange mass and pair coefficients, Cabibbo angle, physical speed and strong length.',
        'physical_on_shell_decay_amplitude_proved':False,'physical_Lambda_lifetime_proved':False,
        'remaining_premise':'Transfer this actual source current to the full five-body interacting scattering spectral projector at the open shell, with the complete weak Hamiltonian and its vacuum/other channels. A compact N5 response is not an on-shell p-pi decay amplitude.',
        'scope':'Positive native current transfer from the complete interacting bulk Lambda pole, with an independently enclosed same-source p-pi energy shell. The kinematic shell and nonzero compact vertex are separate statements.'}
