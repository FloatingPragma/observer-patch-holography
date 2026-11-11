"""Validate eight finite-observable recommendations against their exact receipts."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
PACKAGES=['finite_observer_readout_closure','finite_current_readout_closure','finite_observable_closures']
IDS={'exec-modular-phase','exec-joint-clocks','exec-temporal-hilbert',
     'exec-axial-circulation','exec-integer-lift','exec-coherent-screen',
     'exec-incident-fields','measurement'}


def require(condition, message):
    if not condition:raise ValueError(message)


def authenticate(path):
    receipt=json.loads(path.read_text())
    bindings=receipt.get('source_bindings',receipt.get('input_sha256'))
    require(isinstance(bindings,dict) and bool(bindings),'missing evidence bindings')
    for name,digest in bindings.items():
        source=ROOT/name
        require(source.is_file(),f'missing source: {name}')
        require(hashlib.sha256(source.read_bytes()).hexdigest()==digest,f'changed source: {name}')
    return receipt


def validate_updates(snapshot,updates):
    rows=snapshot['observables'];lookup={x['id']:x for x in rows}
    require(len(lookup)==len(rows),'duplicate baseline row')
    require(len(rows)==235 and sum(x['status'] in {'theorem','exectheorem'} for x in rows)==50,
            'captured baseline must contain 235 observables and 50 completed rows')
    require(len(updates)==8 and {x['id'] for x in updates}==IDS,'unreviewed or duplicate promotion')
    for row in updates:
        require(row['id'] in lookup,'missing dashboard target')
        require(lookup[row['id']]['status']=='computation','baseline status differs')
        require(row['from_status']=='computation' and row['status']=='exectheorem','unsupported status change')
        require(row['physical_premise_promotions']==[],'physical premise promotion is outside this closure')
        require(bool(row['scope']) and bool(row['evidence']),'missing scope or evidence')
    return lookup


def build():
    snapshot=json.loads((HERE/'dashboard_inputs.json').read_text())
    receipts={name:authenticate(ROOT/'model'/name/'results.json') for name in PACKAGES}
    updates=[]
    def add(identifier,package,section,summary,scope,**extra):
        why=extra.pop('why',summary)+' '+scope
        updates.append(dict(id=identifier,from_status='computation',status='exectheorem',
            status_label='Theorem + executed',summary=summary,scope=scope,why=why,
            physical_premise_promotions=[],evidence=[f'model/{package}/PROOF.md',f'model/{package}/results.json'],
            receipt_section=section,**extra))
    r=receipts['finite_observer_readout_closure']
    require(r['scope']['physical_time_or_Universe_attachment_closed'] is False,'physical scope changed')
    p=r['periodic_Y']
    require(p['full_detector_Kraus_identity_exact'] and p['all4356_matrix_units_of_stopped_map_trace_identity_exact'],
            'periodic detector/source theorem missing')
    require(p['physical_time_identified'] is False and p['native_attempts_generate_modular_phase'] is False,'phase claim broadened')
    add('exec-modular-phase','finite_observer_readout_closure','periodic_Y',
        'The full 66-state stopped source and complete Y detector give the declared modular-intervention curve. Finite original records at caps 4, 32, 128 retain both herald signs, preparation failure and timeout; tail is bounded by (5/33)(13/15)^h.',p['recommended_scope'],
        value='omega=log(5/4); amplitudes 1/140 (C), 1/220 (E)')
    p=r['joint_event_clocks']
    require(p['full_work_operator_dimension']==16 and p['arbitrary_reference_extension_from_full_map_identity'],
            'joint clock all-input map missing')
    require(p['physical_affine_clock_map_proved'] is False,'physical clock claim broadened')
    add('exec-joint-clocks','finite_observer_readout_closure','joint_event_clocks',
        'Exact equality on all 16 work-operator coordinates binds the count-generating instruments to complete native histories and reference feedback. Mean relative count 2 and variance 16 coexist with field-plus probabilities 35/54 and 14/27.',p['recommended_scope'])
    p=r['unused_lag_correlations']
    require(p['full_source_intertwiner_KC_equals_CR'] and p['full_source_annihilator_vector_zero']
            and p['positive_definite_Gram_and_strict_positive_contraction_exact'],'full-source recurrence proof missing')
    require(p['sampled_data_certifies_exact_rank'] is False,'sampled rank improperly promoted')
    add('exec-temporal-hilbert','finite_observer_readout_closure','unused_lag_correlations',
        'The literal 66-by-7 inclusion satisfies KC=CR and C^T M C=G>0. The source annihilator proves the recurrence for every lag; independent source powers verify lags 15 through 32.',p['recommended_scope'])
    r=receipts['finite_current_readout_closure'];p=r['axial_counter']
    require(p['conditional_axial_mean_zero_all_states'] and p['independent_two_attempt_cross_covariance_zero'],
            'axial source moments missing')
    require(sum(map(sum,p['complete_two_attempt_word_counts']))==p['common_word_denominator']
            and p['complete_two_attempt_record_count']==32400,'incomplete axial record law')
    add('exec-axial-circulation','finite_current_readout_closure','axial_counter',
        'Actual native currents have conditional axial mean zero and stationary covariance (4N/99) I3. All 32400 two-attempt words are retained. Exact four-step minorization 11/90000 gives bracket mean-square error at most min(1,1440000/(11N)); the verified martingale hypotheses imply the Brownian counter limit.',
        'Stationary preparation for exact finite covariance; every initial classical binary law for the functional limit. Scaled native count is neither physical position nor laboratory time.')
    p=r['integer_face_lift']
    require(all(p['reference_checks'].values()) and p['complete_load_coin_class_contexts']==2592,'integer source replay missing')
    require(p['archive_census_reexecuted'] is False and p['winding_is_conserved_particle_hypercharge'] is False,'integer census or particle claim broadened')
    add('exec-integer-lift','finite_current_readout_closure','integer_face_lift',
        'Kernel-checked integer incidence identities give m=B n, m mod6=class and zero total on each closed carrier. All 2592 native load/coin/class contexts and the 80-carrier, 200-sweep reference attach the theorem to accepted twists, rejected branches and consistent cuts. The large-run 99.7-percent census is an imported diagnostic.',
        'Promote the exact lift and Stokes identities only. The 99.7-percent census is not a theorem or a conserved particle hypercharge.',
        name='Integer lift identities of the face register',value='m mod6=class; total lift=0',required_wording_change=True,
        why='Exact incidence and accepted-twist identities attach integer lifts to the declared native source. The imported 99.7-percent population census is a separate diagnostic; no conserved physical hypercharge or particle identity is proved.')
    r=receipts['finite_observable_closures']
    require(r['all_checks_pass'] is True and r['physical_premise_promotions']==[],'finite readout checks failed')
    expected={'exec-coherent-screen','exec-incident-fields','measurement'}
    require({x['dashboard_id'] for x in r['rows']}==expected,'finite readout targets changed')
    for row in r['rows']:
        flags=row['acceptance_flags']
        for key in ['same_declared_finite_source','complete_original_records_and_failures_retained',
                    'bounded_observable_theorem_proved','fresh_native_source_observable_executed',
                    'finite_construction_premises_closed','dashboard_finite_status_recommendation_accepted']:
            require(flags.get(key) is True,f'missing finite acceptance: {key}')
        require(flags['physical_identification_or_hardware_accuracy_claimed'] is False,'physical scope changed')
        identifier=row['dashboard_id']
        if identifier=='exec-coherent-screen':
            require(row['untouched_reference_marginal_unchanged'] and row['reference_alone_distinguishes_C_and_E'] is False,
                    'no-signalling correction missing')
            summary='The record-corrected joint screen-reference interference readout has completed P_C(+)=61/72 and P_E(+)=809/1089, difference 101/968. Both outcomes of all 21 capped words are retained. The untouched reference marginal is I2/2 for both laws.'
            extra=dict(name='Joint screen-reference interference under feedback',required_wording_change=True,
                       value='P_C(+)=61/72; P_E(+)=809/1089',
                       why='A joint screen-reference observable distinguishes the two declared feedback laws after record correction. The untouched reference marginal remains I2/2 under both laws and cannot distinguish them by itself.')
        elif identifier=='exec-incident-fields':
            summary='The actual 240-port, 28680-state native transfer executes a cap-one four-read law with 124 original records per policy, including 60 timeout records. Integer numerators sum exactly to the common denominator; no rounded probability remainder is discarded. The longer-cap density and feedback comparisons retain their separate scope.'
            extra={}
        else:
            summary='The freshly executed finite interacting source retains all 4 original terminal outcomes. Exact conversion of every binary host component and an independent rational Taylor enclosure certify complete half-trace error below 3e-7 and a nonzero unconditional feedback contrast, without branch renormalization.'
            extra={}
        add(identifier,'finite_observable_closures',f'rows/{identifier}',summary,
            'Theorem and execution for the declared finite source, preparation, apparatus and retained timeout/error law; no physical identification or hardware accuracy claim.',**extra)
    lookup=validate_updates(snapshot,updates)
    before=Counter(x['status'] for x in lookup.values())
    after=before.copy()
    for row in updates:after[row['from_status']]-=1;after[row['status']]+=1
    total=len(lookup);old=before['theorem']+before['exectheorem'];new=after['theorem']+after['exectheorem']
    paths=[HERE/'dashboard_inputs.json',HERE/'bridge.py',HERE/'PROOF.md',HERE/'test_bridge.py']
    paths +=[ROOT/'model'/name/'results.json' for name in PACKAGES]
    return dict(schema='oph.observable_score_bridge.v1',source_url=snapshot['source_url'],
        source_bindings={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        all_checks_pass=True,baseline=dict(total=total,completed=old,completed_percent=100*old/total,status_counts=dict(before)),
        proposed=dict(total=total,completed=new,completed_percent=100*new/total,status_counts=dict(after),
                      completion_percentage_point_change=100*(new-old)/total),updates=updates,
        website_modified=False,physical_completion_gates_accepted=[],new_empirical_observations=0,
        calibrated_particle_mass_promotions=[],all_eight_require_the_stated_finite_scope=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result=build();data=json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n'
    if args.write:(HERE/'results.json').write_text(data)
    if args.check:require((HERE/'results.json').read_text()==data,'score bridge is stale')
    print(json.dumps(dict(baseline=result['baseline']['completed'],proposed=result['proposed']['completed'],total=result['baseline']['total'])))


if __name__=='__main__':main()
