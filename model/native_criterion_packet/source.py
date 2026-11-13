"""Literal charge-criterion candidate and all31 retained authority rows.

The headline experiment executes the shared charged source, not a saved
probability receipt. Architectural premises are explicit and permitted.
Physical identification and external adjudication are separate predicates.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json
from model.native_gauss_walk import source as joint
from model.native_electric_gauss import source as electric
from model.native_observer_program.contract import descriptor_sha256

ROOT=Path(__file__).resolve().parents[2]


def equal_integer_halves(q):
    if type(q)is not int:raise ValueError('integer elementary-charge coordinate required')
    return None if q%2 else (q//2,q//2)


def cut_records(wave,region):
    region=tuple(sorted(region))
    if len(set(region))!=len(region)or any(type(v)is not int or not 0<=v<12 for v in region):
        raise ValueError('a subset of the twelve native ports is required')
    joint_weights=defaultdict(int)
    for(c,event),amplitude in wave.amplitudes.items():
        s=joint.state(c,wave.charge)
        charge=sum(electric.charge_density(s)[v]for v in region)
        flux=electric.cut_flux(s,sum(1<<v for v in region))
        if charge!=flux:raise ArithmeticError('measured integer cut flux differs from source charge')
        joint_weights[charge,flux]+=joint.norm(amplitude)
    denominator=joint.NORM*wave.denominator**2
    if sum(joint_weights.values())!=denominator:raise ArithmeticError('complete charge records lost probability')
    return {'region':list(region),'complete_joint_charge_flux_records':[
        {'charge_in_q_star':q,'outward_flux_in_q_star':f,'probability':str(F(n,denominator))}
        for(q,f),n in sorted(joint_weights.items())],
        'probability_sum':'1','half_integer_charge_probability':'0',
        'charge_and_flux_commute':True,'postselection_or_renormalization':False}


def run():
    value=joint.prepare();records=[]
    for frame in range(3):
        records.append({'frame':frame,'all_twelve_single_port_cuts':[cut_records(value,(v,))for v in range(12)],
            'multi_port_cuts':[cut_records(value,r)for r in ((0,1,2),(3,4,5,6),(7,8,9,10,11),tuple(range(12)))],
            'norm_squared':str(value.probability()),'complete_joint_arcs':len(value.amplitudes)})
        if frame<2:value=joint.step(value,electric_phase=True)
    initial=joint.prepare()
    for _ in range(2):value=joint.step(value,inverse=True,electric_phase=True)
    if not joint.equivalent(value,initial):raise ArithmeticError('headline source/link process failed its full inverse control')
    controls=[{'charge_coordinate':q,'equal_integer_halves':equal_integer_halves(q)}for q in(-3,-2,-1,0,1,2,3)]
    authority=json.loads((ROOT/'model/openwave_m10_review_gate/results.json').read_text())
    status=authority['official_status']
    if len(status['criteria'])!=31:raise ArithmeticError('complete31-row authority table required')
    return {'schema':'oph.native_charge_criterion_packet.v1',
        'verdict':'PASS_EXECUTED_INTEGER_CUT_CHARGE_CANDIDATE_WITH_PHYSICAL_GATE_OPEN',
        'common_program_sha256':descriptor_sha256(),
        'headline_command':'python -m model.native_criterion_packet --check',
        'preparation':'primary matching words0/1; q1/character1; canonical neutral string; complete sixty-event state',
        'dynamics':'two actual joint native source/tag/link steps with electric phase enabled; no intermediate charge measurement',
        'terminal_experiments':records,'full_inverse_control':True,
        'splitting_controls':controls,
        'criterion1':{'literal_test':'integer topological charge; a split e -> e/2 + e/2 impossible',
            'finite_result':'The integer Gauss flux equals regional charge as a complete commuting record; an odd integer charge cannot be two equal integer charges.',
            'topological_scope':'Relative integer-chain boundary is unchanged by adding any closed integer flux cycle. This is a discrete Gauss boundary charge, not an identified electron soliton winding.',
            'arbitrary_architecture_allowed':True,
            'unit_ancestry':'q_star is declared; no nonzero real additive charge is obtained from the Z6 character alone. Electron q=1 or q=3 are possible assignments, not established identities.',
            'even_charge_negative_control':'q=2 permits1+1; the exclusion is proved for odd q only.',
            'missing_physical_identification':['identify a persistent physical electron state in this joint construction','identify its measured charge with an odd q_star multiple','review the declared native extension or establish its in-platform implementation'],
            'local_review_recommendation':'finite charge-quantization candidate ready for review; no unconditional validated claim',
            'physical_acceptance_complete':False,'official_award_requested_or_recorded':False},
        'official_authority':{'upstream_commit':authority['upstream_commit'],
            'all31_criteria':status['criteria'],'counts':status['upstream_counts']['M10'],
            'rank':status['oph_official_rank'],'local_receipts_change_authority':False},
        'scope':{'declared_simulator_architecture_requires_uniqueness_proof':False,
            'hardware_configurability_is_a_physics_validation':False,
            'new_physical_criteria_validated':0,'external_submission_sent':False,
            'shared_full_particle_spectrum_proved':False}}
