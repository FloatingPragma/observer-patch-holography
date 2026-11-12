"""Complete stationary charge, energy, position and native profile records."""
import sympy as s
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import outward
from .isolation import data
from .filter import design
from .local_profile import certificate as local_certificate


def certificate():
    d=data();f=design();profile=local_certificate()
    eta2=d['eta']**2
    qlo,qhi=f['preparation_success_interval'];hlo=1-f['outside_band_probability_upper']
    wlo,whi=(1-eta2)/8,s.Rational(1,8)
    q,h,w=s.symbols('q h w');p=s.symbols('p0:6')
    probabilities=[1-q,q*(1-h),q*h*(1-8*w)]
    probabilities += [q*h*w*a for _ in range(8) for a in p]
    assert len(probabilities)==51 and s.expand(sum(probabilities).subs(p[-1],1-sum(p[:-1])))==1
    rows=[{'record':'complete_preparation_failure','probability':'1-q'},
        {'record':'accepted_preparation_energy_outside_known_neutral_band','probability':'q(1-h)'},
        {'record':'inside_band_complete_atomic_profile_complement','probability':'qh(1-8w)'}]
    for cell in range(8):
        for index,record in enumerate(profile['complete_gauge_invariant_native_occupation_records']):
            plo,phi=map(s.Rational,record['probability_interval'])
            rows.append({'record':'native_neutral_position_and_occupation','cell':cell,
                'native_N_ec':record['native_N_ec'],'native_N_L':record['native_N_L'],
                'scalar_excited':record['scalar_excited'],
                'all_local_matter_plus_scalar_hypercharges':'0 before the occupation readout',
                'probability':f'q h w p{index}',
                'original_input_probability_interval':outward((qlo*hlo*wlo*plo,qhi*whi*phi),600)})
    assert len(rows)==51
    return {'initial_state':'I/D on the complete kinematic source; no atomic or interacting vacuum is substituted',
        'stationary_accepted_state':'Pi_Gauss,L7,B0,k0 M(H)^2 / Tr(Pi_Gauss,L7,B0,k0 M(H)^2)',
        'actual_source_clock':'exp(-it H), all stationary record laws hold for every real source time',
        'complete_energy_readout':['known neutral band','entire energy complement'],
        'complete_position_readout':['eight P_y: neutral local eigenprojector at y, true local grounds elsewhere, blank links',
                                     'I-sum_y P_y'],
        'occupation_readout_order':'after P_y, gauge-invariant N_ec and scalar blank/excited are resolved; their backaction is retained',
        'occupation_readout_preserves_band_energy_claimed':False,
        'global_raw_L_B3_and_total_hypercharge':[7,0,0],
        'no_external_charged_reference_needed':True,
        'conditional_band_position_each_interval':outward((wlo,whi),160),
        'conditional_band_atomic_profile_complement_upper':str(eta2),
        'conditional_band_zero_matter_current_record_lower':str(1-eta2),
        'conditional_energy_bin_success_lower':str(hlo),
        'original_input_preparation_success_interval':[str(qlo),str(qhi)],
        'original_input_success_decimal_interval':[float(qlo),float(qhi)],
        'records':rows,'complete_record_count':51,'exact_probability_sum':'1',
        'all_failures_reference_environments_and_nonblank_virtual_components_retained':True,
        'stationarity_is_closed_source_not_external_environment_lifetime':True}
