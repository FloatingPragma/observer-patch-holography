"""Exact native phase-readout boundary and a retrospective scale test."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

from model.openwave_native_process import source as native
from model.openwave_native_response.source import coherence_multiplier

ROOT=Path(__file__).resolve().parents[2]


def log_decay_interval(multiplier,terms=12):
    """Positive rational enclosure of -log(multiplier), not an energy."""
    lam=F(multiplier)
    if not 0<lam<=1 or type(terms) is not int or terms<1:
        raise ValueError('positive contraction and positive integer Taylor order required')
    x=1-lam
    lower=sum((x**k/k for k in range(1,terms+1)),F())
    upper=lower+x**(terms+1)/((terms+1)*(1-x))
    return lower,upper


def complete_phase_blind_control():
    """Every3-event branch, including padding, on the full3-port binary carrier."""
    seams=((0,1),(1,2),(0,2));branches=native.raw_local_branches()
    total_by_basis=[F()]*8;patterns=Counter();nonzero=0
    for word in product(branches,repeat=3):
        weight=F(1);mapping={i:i for i in range(8)}
        for edge,b in zip(seams,word,strict=True):
            weight*=b.weight
            mapping={i:j for i,x in mapping.items()
                     if (j:=native.apply_branch(x,edge,b)) is not None}
        assert len(set(mapping.values()))==len(mapping)
        domain=tuple(sorted(mapping))
        patterns[domain]+=1
        for i in domain:total_by_basis[i]+=weight
        if weight and domain:nonzero+=1
        # K^*K is weight times a diagonal domain projection on every branch.
        for i,j in product(range(8),repeat=2):
            gram=weight*int(i in mapping and j in mapping and mapping[i]==mapping[j])
            assert gram==(weight if i==j and i in mapping else 0)
        # |100>+|010> and |100>+i|010> give identical full-word probabilities.
        def probability(amplitudes):
            output={}
            for i,(re,im) in amplitudes.items():
                if i in mapping:
                    old=output.get(mapping[i],(0,0))
                    output[mapping[i]]=(old[0]+re,old[1]+im)
            return weight*F(sum(re*re+im*im for re,im in output.values()),2)
        plus_probability=probability({1:(1,0),2:(1,0)})
        phase_probability=probability({1:(1,0),2:(0,1)})
        assert plus_probability==phase_probability
    assert total_by_basis==[F(1)]*8
    assert sum(patterns.values())==25**3
    analyzer_plus=F((1+1)**2,4);analyzer_phase=F(1**2+1**2,4)
    assert analyzer_plus==1 and analyzer_phase==F(1,2)
    return {'complete_raw_histories':25**3,'nonzero_history_maps':nonzero,
            'complete_input_operator_columns':64,'total_effect_diagonal_exact':list(map(str,total_by_basis)),
            'domain_pattern_histories':[{'domain':list(k),'history_count':v} for k,v in sorted(patterns.items())],
            'phase_preparations':'(|port0>+|port1>)/sqrt2 versus (|port0>+i|port1>)/sqrt2',
            'every_native_record_probability_difference_exact':'0',
            'external_analyzer_positive_control':{'operation':'supplied Hadamard then port measurement, absent from native repair/hold grammar',
                'plus_probability_exact':str(analyzer_plus),'phase_probability_exact':str(analyzer_phase)}}


def source_band(level):
    g=native.graph(level);M=len(g.seams)
    classes={}
    for i in range(g.ports):
        for j in range(i+1,g.ports):
            key=(len(g.neighbours[i])+len(g.neighbours[j]),j in g.neighbours[i])
            if key not in classes:
                classes[key]=(i,j,coherence_multiplier(g,i,j))
    multipliers=[v[2] for v in classes.values()]
    minimum_decay=1-max(multipliers);maximum_decay=1-min(multipliers)
    # x <= -log(1-x) <= x/(1-x), so this is a uniform exact rate-ratio cap.
    ratio_upper=maximum_decay/(minimum_decay*(1-maximum_decay))
    tested=[]
    for i,j,lam in sorted(classes.values(),key=lambda r:r[2],reverse=True)[:3]:
        lo,hi=log_decay_interval(lam)
        tested.append({'ports':[i,j],'event_multiplier_exact':str(lam),
            'two_event_multiplier_exact':str(lam**2),'four_event_multiplier_exact':str(lam**4),
            'phase_per_attempt_interval_exact':['0','0'],
            'decay_rate_per_attempt_interval_exact':list(map(str,(lo,hi))),
            'decay_rate_per_intensity_time_interval_exact':list(map(str,(2*M*lo,2*M*hi))),
            'nonzero_stable_oscillatory_pole':False})
    return {'level':level,'ports':g.ports,'seams':M,
            'stationary_mode_phase_interval_exact':['0','0'],
            'coherence_multiplier_extrema_exact':list(map(str,(min(multipliers),max(multipliers)))),
            'all_coherence_decay_rate_ratio_upper_exact':str(ratio_upper),
            'source_to_any_unitary_channel_trace_distance_lower_exact':str(F(min(map(len,g.neighbours)),2*M)),
            'tested_modes':tested}


def run():
    parent=json.loads((ROOT/'model/openwave_common_mass_source/results.json').read_text())
    unit=F(parent['electron_unit_GeV_exact'])
    bands=[source_band(level) for level in (0,1,2)]
    upper=max(F(x['all_coherence_decay_rate_ratio_upper_exact']) for x in bands)
    inherited=[r for r in parent['same_source_free_Dirac_branches'] if r['channel'] in ('electron','muon','tau')]
    return {'schema':'oph.native_rest_energy.v1',
        'verdict':'PASS: fixed-source rest-phase test exposes the missing energy/readout transfer',
        'complete_phase_blind_control':complete_phase_blind_control(),
        'native_source_tests':bands,
        'rest_energy_test':{'frozen_source':'fair C/E, fixed uniform seam selection, global one-load carrier; no mass-dependent parameter',
            'stationary_and_nonstationary_phase_frequency_interval_per_attempt_exact':['0','0'],
            'nonzero_stable_pole_found':False,'physical_rest_mass_interval_GeV':None,
            'energy_clock_attachment_proved':False,
            'reason':'Every nonstationary eigenvalue of this unread finite channel is a real positive contraction. Its phase is zero and its amplitude decays. Native-only repair/hold record histories cannot measure source offdiagonal phases.',
            'Hamiltonian_boundary':'Phi-I is self-adjoint on Hilbert-Schmidt operator space; -i[H,.] is skew-adjoint. Equality would force both zero, whereas the actual repair has nonzero population transfer.',
            'alternative_transfer_boundary':'A Euclidean transfer construction using -log(Phi), a coherent dilation or a supplied source-word service is a distinct possible attachment. It requires an explicit physical time, algebra/state/readout and reconstruction transfer; positive real eigenvalues alone neither establish nor forbid it.',
            'record_scope':'This concerns the declared native instruction/readout grammar. A supplied coherent analyzer or another architecture changes the contract.',
            'work_and_clock':'Event ordinals, preparation, blank records and their resource counts are declared. No physical Hamiltonian, physical seconds, zero-work claim or free-reset cost is inferred.'},
        'retrospective_counterfactual':{'description':'Treating these one-load coherence relaxation rates as masses would be an additional unjustified prescription. Even that prescription cannot reproduce the inherited charged-lepton hierarchy.',
            'assignment':'slowest nonzero coherence relaxation rate is assigned the inherited electron anchor; all other coherence rates use that same conversion',
            'electron_anchor_GeV_exact':str(unit),
            'coherence_band_mass_ratio_upper_exact':str(upper),
            'counterfactual_band_GeV_upper_exact':str(unit*upper),
            'counterfactual_band_GeV_upper_decimal':float(unit*upper),
            'these_are_physical_rest_masses':False,
            'inherited_calibrated_channels':inherited,
            'muon_tau_hierarchy_reproduced':False,
            'target_exposure':'The electron, muon and tau targets and inherited calibrated source coefficients were available before this calculation; this comparison is retrospective.',
            'prospective_test':'Freeze a native coherent operation/readout transfer and an independently defined Hamiltonian/time calibration before comparing one unused species and multiple packet widths. The present grammar has no such attachment.'},
        'scope':{'native_phase_readout_boundary_proved':True,'universal_OPH_no_go':False,
            'new_GeV_prediction':False,'new_Nature_validation':False}}
