#!/usr/bin/env python3
"""Run checked finite compositions, adversarial joins and scoped claim gates."""
from __future__ import annotations
import os
for variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
from pathlib import Path
from fractions import Fraction as F
from dataclasses import replace, asdict
import hashlib
import importlib.util
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
CHAIN = HERE.parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
from strict_facade import StrictClock, StrictParticle, first_event, STATE_TOL
from gates import (Source, Endpoint, ActualQuantum, actual_quantum, quantum_endpoint,
                   join_quantum, field_bound, require_same_source, require_trajectory_clock, classical_witness)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def close(actual, expected, label, tolerance=3e-12):
    error = float(np.max(abs(np.asarray(actual)-np.asarray(expected))))
    if not np.isfinite(error) or error > tolerance:
        raise AssertionError(f'{label}: {error}')
    return error


def reject(name, function, output):
    try:
        function()
    except (ValueError, TypeError) as exc:
        output[name] = {'rejected':True, 'reason':str(exc)}
    else:
        raise AssertionError('mutation admitted: '+name)


def pins():
    manifest = read(HERE/'INPUT_PINS.json')
    for relative, expected in manifest['files'].items():
        if sha(ROOT/relative) != expected:
            raise ValueError('frozen input changed: '+relative)
    return manifest['files']


def quantum_run():
    print('Executing complete coarsened quantum instrument', flush=True)
    module = load('integration_particle', CHAIN/'particle/particle_chain.py')
    chain = module.ParticleChain()
    raw = chain.clock
    clock, detector = StrictClock(raw), StrictParticle(chain)
    rho, branches = raw.preparation()
    cc, pc = read(CHAIN/'clock/CLOCK_CONTRACT.json'), read(CHAIN/'particle/PARTICLE_CONTRACT.json')
    left, right = quantum_endpoint(cc,cc), quantum_endpoint(pc,cc)
    initial = rho.copy()
    split1 = clock.count_split(rho,.37,.8)
    witness = actual_quantum(raw,rho,split1)
    assert join_quantum(left,right,witness)
    state1 = sum(split1)
    q1, stops1, residual1 = detector.measure(state1)
    continued = sum(q1)
    split2 = clock.count_split(continued,.37,.8)
    q2, stops2, residual2 = detector.measure(sum(split2))
    terminal = sum(q2)
    coeff = detector.energy_coefficients(terminal)
    with np.load(CHAIN/'particle/particle_states.npz',allow_pickle=False) as archive:
        complete = archive['terminal_energy_coefficients'].sum(axis=0)
        instrument_keys=[key for key in chain.arrays if key.startswith('abort_') or key.startswith('C_Q')]
        instrument_error=max(close(chain.arrays[key],archive[key],'compiled full detector Kraus map') for key in instrument_keys)
    coefficient_error = close(coeff,complete,'coarsened complete energy poststates')
    probabilities = np.trace(coeff,axis1=1,axis2=2).real
    stops = sum(stops1.values())+sum(stops2.values())
    complete_mass_error = close(probabilities.sum()+stops,1.,'complete stop and energy law')
    frozen = read(CHAIN/'particle/RESULT.json')
    stop_error = close(stops,frozen['unconditional']['all_apparatus_stop_probability'],'all stops')
    assert np.array_equal(rho,initial)
    assert all(np.trace(state).real > 0 for state in split1)
    inspected=[rho,*[s for _,s in branches],*split1,*q1,*split2,*q2]
    with np.load(CHAIN/'clock/common_clock_states.npz',allow_pickle=False) as archive:
        inspected += [archive[row['clock_archive_key']] for row in frozen['clock_inputs']]
    off_max=0.
    for state in inspected:
        off=state.copy()
        for ids,*_ in raw.sectors.values():off[np.ix_(ids,ids)]=0
        off_max=max(off_max,float(np.max(abs(off))))
    if off_max != 0:raise ValueError('declared composition input is not exactly number-block diagonal')
    return (chain,clock,rho,left,witness), {
        'source_dimension':384, 'actual_preparation_branches':len(branches),
        'exact_number_block_inputs_inspected':len(inspected),'maximum_actual_off_number_block_entry':off_max,
        'clock_intervals':2, 'T_per_interval':.37, 'native_rate':.8,
        'count_classes':['N=0','N=1','N>=2 complete complement'],
        'first_count_probabilities':[float(np.trace(s).real) for s in split1],
        'second_count_probabilities_before_detector':[float(np.trace(s).real) for s in split2],
        'detector1_stops':stops1,'detector2_stops':stops2,
        'terminal_group_probabilities':probabilities.tolist(),
        'full_detector_Kraus_families_compared':len(instrument_keys),'full_detector_Kraus_map_max_error':instrument_error,
        'terminal_group_count':len(coeff), 'complete_probability':float(probabilities.sum()+stops),
        'frozen_576_by_348_poststate_coarsening_error':coefficient_error,
        'frozen_stop_probability_error':stop_error, 'complete_mass_error':complete_mass_error,
        'detector_completeness_errors':[residual1,residual2],
        'coarsening':'Sum all four reservoir labels, both complete count partitions and regional histories only after their linear maps; retain each apparatus stop family and all 348 terminal energy groups. Full fine record custody is pinned.',
        'proper_code_transfer':'Declared reference-preserving coherent Gauss-code load/export and routing service; no fine F or B measurement.',
        'energy_apparatus':'Declared grouped spectral Luders instrument; no native spectral compiler.'}


def guard_checks(objects):
    print('Testing admission guards against frozen comparators', flush=True)
    chain, clock, rho, _, _ = objects
    raw = chain.clock
    old = load('integration_frozen_clock', CHAIN/'audit_clock_force/frozen_common_clock.py')
    fixture = read(CHAIN/'audit_clock_force/MALFORMED_INPUTS.json')
    rejected = {}
    class Trap:
        n, sectors, ops = raw.n, raw.sectors, raw.ops
        def __getattr__(self,name):
            raise AssertionError('raw helper called before validation: '+name)
    checked = StrictClock(Trap())
    bad = {}
    for name,item in fixture['clock_sparse_fixtures'].items():
        state = np.zeros(tuple(fixture['clock_matrix_shape']),complex)
        for i,j,re,im in item['entries_row_column_real_imaginary']:
            state[i,j] = complex(float(re),float(im))
        bad[name] = state
    small_cross = np.full((384,384),-9e-14,complex)
    for ids,*_ in raw.sectors.values():small_cross[np.ix_(ids,ids)]=0
    bad['cumulative_small_offblock_negative'] = small_cross
    bad['wrong_shape'] = np.zeros((3,3))
    bad['nonhermitian'] = rho.copy(); bad['nonhermitian'][6,10] += .1j
    bad['supernormalized'] = 2*rho
    bad['infinite'] = rho.copy(); bad['infinite'][0,0] = np.inf
    for name,state in bad.items():
        for method in ('evolve','native','semigroup','count_split','trajectory'):
            arguments = {'evolve':(state,.37),'native':(state,), 'semigroup':(state,.37,.8),
                         'count_split':(state,.37,.8), 'trajectory':(state,.37,.8,1.6,[])}[method]
            reject(name+'_'+method,lambda m=method,a=arguments:getattr(checked,m)(*a),rejected)
    cases = {
        'negative_count_duration':lambda:checked.count_split(rho,-.1,.8),
        'nan_count_rate':lambda:checked.count_split(rho,.1,float('nan')),
        'negative_evolve_scale':lambda:checked.evolve(rho,.1,-1),
        'negative_semigroup_rate':lambda:checked.semigroup(rho,.1,-1),
        'infinite_idle_rate':lambda:checked.semigroup(rho,.1,.8,idle_rate=np.inf),
        'attempt_mark_above_one':lambda:checked.semigroup(rho,.1,.8,attempt_mark=2),
        'negative_move_mark':lambda:checked.semigroup(rho,.1,.8,move_mark=-1),
        'nan_idle_mark':lambda:checked.semigroup(rho,.1,.8,idle_mark=np.nan),
        'nan_trajectory_duration':lambda:checked.trajectory(rho,np.nan,.8,1.6,[]),
        'negative_trajectory_duration':lambda:checked.trajectory(rho,-.1,.8,1.6,[]),
        'zero_native_trajectory_rate':lambda:checked.trajectory(rho,.37,0,1.6,[]),
        'subnormalized_trajectory':lambda:checked.trajectory(rho/2,.37,.8,1.6,[]),
        'invalid_event_time':lambda:checked.trajectory(rho,.37,.8,1.6,[(np.nan,'idle')]),
        'unordered_events':lambda:checked.trajectory(rho,.37,.8,1.6,[(.2,'idle'),(.1,'idle')]),
        'invalid_native_mark':lambda:checked.trajectory(rho,.37,.8,1.6,[(.1,-1)]),
        'zero_probability_idle':lambda:checked.trajectory(rho,.37,.8,.8,[(.1,'idle')])}
    for name,fn in cases.items():reject(name,fn,rejected)
    cross = np.zeros((384,384),complex); cross[0,0]=cross[6,6]=cross[0,6]=cross[6,0]=.5
    reject('cross_number_semigroup',lambda:checked.semigroup(cross,.37,.8),rejected)
    sub = raw.preparation()[1][0][1].astype(complex)
    equivalence = {}
    for name,state in [('normalized',rho),('subnormalized_reservoir',sub),('zero',np.zeros_like(rho))]:
        before=state.copy()
        equivalence[name+'_evolve'] = close(clock.evolve(state,.04),old.CommonClock.evolve(raw,state,.04),'frozen evolve',0.)
        equivalence[name+'_native'] = close(clock.native(state),old.CommonClock.native(raw,state),'frozen native',0.)
        equivalence[name+'_semigroup'] = close(clock.semigroup(state,.04,.8,attempt_mark=.8,idle_rate=.8,idle_mark=.7),
            old.CommonClock.semigroup(raw,state,.04,.8,attempt_mark=.8,idle_rate=.8,idle_mark=.7),'frozen tilted semigroup',0.)
        equivalence[name+'_count_split'] = close(clock.count_split(state,.04,.8),old.CommonClock.count_split(raw,state,.04,.8),'frozen count split',0.)
        assert np.array_equal(state,before)
    equivalence['cross_number_evolve'] = close(clock.evolve(cross,.04),old.CommonClock.evolve(raw,cross,.04),'cross-sector unitary',0.)
    equivalence['cross_number_native'] = close(clock.native(cross),old.CommonClock.native(raw,cross),'cross-sector native',0.)
    a,al = clock.trajectory(rho,.04,.8,1.6,[(.02,'idle')])
    b,bl = old.CommonClock.trajectory(raw,rho,.04,.8,1.6,[(.02,'idle')])
    equivalence['trajectory_state'] = close(a,b,'frozen trajectory',0.)
    assert al == bl
    # The detector facade has the same input validation and unchanged invocation.
    qa,sa,_ = StrictParticle(chain).measure(sub)
    qb,sb,_ = chain.measure(sub)
    equivalence['particle_measure'] = close(qa,qb,'particle facade',0.)
    assert sa == sb
    for name,state in bad.items():
        reject('particle_'+name,lambda s=state:StrictParticle(chain).measure(s),rejected)
    # A positive matrix can contain tiny coherences with a large total norm.
    tiny_psd=np.eye(384,dtype=complex)/384-small_cross
    min_eigenvalue=float(np.linalg.eigvalsh(tiny_psd).min())
    assert min_eigenvalue>0
    raw_zero=old.CommonClock.semigroup(raw,tiny_psd,0.,.8)
    discarded_trace_norm=float(np.sum(abs(np.linalg.eigvalsh(tiny_psd-raw_zero))))
    assert discarded_trace_norm>STATE_TOL
    reject('PSD_tiny_cross_number_semigroup',lambda:checked.semigroup(tiny_psd,0.,.8),rejected)
    reject('PSD_tiny_cross_number_count_split',lambda:checked.count_split(tiny_psd,0.,.8),rejected)
    reject('PSD_tiny_cross_number_particle',lambda:StrictParticle(chain).measure(tiny_psd),rejected)
    diagnostics={'PSD_tiny_cross_number':{'minimum_eigenvalue':min_eigenvalue,
        'maximum_offblock_entry':9e-14,'raw_zero_time_discarded_trace_norm':discarded_trace_norm,
        'checked_compressed_apis_reject':True,'full_space_apis_remain_cross_number_capable':True}}
    return rejected,equivalence,diagnostics


def classical_run(rejected,equivalence):
    print('Executing reciprocal field, rates and native force instrument', flush=True)
    force = load('integration_force',CHAIN/'force/force_measurement.py')
    force.dependencies()
    fixed = sys.modules['field_feedback']
    reciprocal = load('integration_reciprocal',force.BRIDGE/'reciprocal_feedback.py')
    field = load('integration_controlled_field',CHAIN/'field_limit/controlled_field.py')
    old = load('integration_frozen_force',CHAIN/'audit_clock_force/frozen_force_measurement.py')
    edges, neighbours = fixed.graph()
    table = reciprocal.pair_table(neighbours)
    native = fixed.literal_native()
    witness=classical_witness(fixed,edges,native)
    reject('classical_wrong_native',lambda:classical_witness(fixed,edges,lambda *args:None),rejected)
    reject('classical_wrong_seam_order',lambda:classical_witness(fixed,list(reversed(edges)),native),rejected)
    laws, exposures = [],[]
    class Trap:
        def first_event(self,*args,**kwargs):
            raise AssertionError('raw force helper called before validation')
    for endpoint in (0,1):
        state=(endpoint,9)
        rates=reciprocal.pair_rates(force.counters(state),edges,table)
        law,g,total=first_event(force,state,edges,rates,native)
        ol,og,ot=old.first_event(state,edges,rates,native)
        equivalence['force_endpoint_'+str(endpoint)] = close(law,ol,'frozen first event',0.)
        assert g == og and total == ot
        laws.append(law);exposures.append(g)
    rates=reciprocal.pair_rates(force.counters((0,9)),edges,table)
    zero_individual=rates.copy();zero_individual[0,0]=0
    for label,valid_rates,valid_time in [('zero_horizon',rates,0.),('zero_individual_rate',zero_individual,.001)]:
        a,ag,at=first_event(force,(0,9),edges,valid_rates,native,valid_time)
        b,bg,bt=old.first_event((0,9),edges,valid_rates,native,valid_time)
        equivalence['force_'+label]=close(a,b,'frozen boundary input',0.)
        assert ag==bg and at==bt
    negative=rates.copy();negative[0,0]=-1
    nan=rates.copy();nan[0,0]=np.nan
    for name,r,t in [('negative_rate',negative,.001),('nan_rate',nan,.001),('wrong_rate_shape',rates[:1],.001),
                     ('negative_horizon',rates,-.001),('nan_horizon',rates,np.nan),('all_zero_rate_outside_raw_domain',rates*0,.001)]:
        reject('force_'+name,lambda r=r,t=t:first_event(Trap(),(0,9),edges,r,native,t),rejected)
    reject('force_wrong_state',lambda:first_event(Trap(),(96,9),edges,rates,native),rejected)
    reject('force_wrong_edges',lambda:first_event(Trap(),(0,9),edges[:2],rates,native),rejected)
    variant,variant_laws = force.arithmetic_transfer(edges,neighbours,table,field,reciprocal,native)
    # Re-enter both computed variant laws through the checked raw instrument.
    # Exact local integer energies supply every outgoing rate in the finite run.
    energies={tuple(r['state']):float(F(r['exact_energy'])) for r in variant['compiled_pairs']}
    for endpoint in (0,1):
        state=(endpoint,9); rates=np.ones((2,len(edges)))
        for side,position in enumerate(state):
            for seam,(a,b) in enumerate(edges):
                if position in (a,b):
                    dest=b if position==a else a
                    after=(dest,state[1]) if side==0 else (state[0],dest)
                    rates[side,seam]=np.exp(-2*(energies[after]-energies[state]))
        checked,_,_=first_event(force,state,edges,rates,native)
        close(checked,variant_laws[endpoint],'integer field to checked native force')
    drift=sum(float(l[:,0]@l[:,5])/g/2 for l,g in zip(laws,exposures))
    expected=read(CHAIN/'force/FORCE_CONTRACT.json')
    close(drift,expected['observable']['measured_drift'],'force drift')
    return {'source_dimension':96**2,'quantum':False,'source_witness':witness,'law_shapes':[list(l.shape) for l in laws],
        'complete_probability':[float(l[:,0].sum()) for l in laws],
        'native_holds_retained':[float(l[(l[:,8]==1)&(l[:,4]==0),0].sum()) for l in laws],
        'cap_probabilities':[float(l[l[:,8]==1,0].sum()) for l in laws],
        'drift':drift,'force':float(table[0,9]-table[1,9]),
        'controlled_field_to_rate_transfer':variant,
        'scope':'Exact integer arithmetic controller and exact-real rate coupling certificate; library exponential roundoff lacks a separate certified enclosure.'}


def transfer_checks(objects,rejected):
    print('Testing source, clock, preparation, observable and regulator joins',flush=True)
    chain,_,rho,endpoint,witness=objects
    dynamic=load('integration_dynamic_field',CHAIN/'dynamical_field/dynamical_field.py')
    dynamic.admit_actual_input(chain.clock.h,rho,F(37,100),chain.model)
    family=dynamic.Family(6)
    order=[family.index[(mask,dynamic.centered(f[3],6))] for mask,f in chain.model.states]
    h_error=close(family.h.toarray()[np.ix_(order,order)],chain.clock.h,'actual M6 Hamiltonian')
    actual_ops={tuple(meta[k] for k in ('seam','law','coin','record_value')):op.toarray()
                for meta,op in chain.clock.ops}
    instrument_error=max(close(op.toarray()[np.ix_(order,order)],actual_ops[key],'actual M6 MIX')
                         for key,op in family.instruments())
    variants={
        'wrong_source':replace(endpoint,source=replace(endpoint.source,identity='external_LOCAL12')),
        'wrong_source_array':replace(endpoint,source=replace(endpoint.source,generator='0'*64)),
        'wrong_carrier':replace(endpoint,source=replace(endpoint.source,carrier='classical_probability')),
        'wrong_clock':replace(endpoint,clock=replace(endpoint.clock,duration=F(38,100))),
        'wrong_preparation':replace(endpoint,preparation='global_reference_without_reservoir'),
        'wrong_regulator':replace(endpoint,source=replace(endpoint.source,regulator='Z42_periodic',dimension=2688)),
        'wrong_observable':replace(endpoint,observable='static_Coulomb_potential'),
        'zero_tail':replace(endpoint,complete_counts=('N=0','N=1')),
        'refusal_removed':replace(endpoint,complete_failures=False)}
    for name,other in variants.items():reject('join_'+name,lambda other=other:join_quantum(endpoint,other,witness),rejected)
    cc=read(CHAIN/'clock/CLOCK_CONTRACT.json');pc=read(CHAIN/'particle/PARTICLE_CONTRACT.json')
    for label,mutated in [('clock_id',{**pc,'clock_id':'wrong_clock'}),
                          ('count_tail',{**pc,'count_classes':pc['count_classes'][:2]}),
                          ('stop_family',{**pc,'outcome_partition':{k:v for k,v in pc['outcome_partition'].items() if k!='first_detector_stops'}})]:
        reject('contract_'+label,lambda m=mutated:quantum_endpoint(m,cc),rejected)
    reject('join_labels_only',lambda:join_quantum(endpoint,endpoint,None),rejected)
    reject('join_matching_wrong_clock',lambda:join_quantum(variants['wrong_clock'],variants['wrong_clock'],witness),rejected)
    reject('actual_wrong_clock',lambda:actual_quantum(chain.clock,rho,duration=F(38,100)),rejected)
    reject('join_matching_false_arrays',lambda:join_quantum(variants['wrong_source_array'],variants['wrong_source_array'],witness),rejected)
    reject('actual_reference_for_preparation',lambda:actual_quantum(chain.clock,chain.model.reference()[0].astype(complex)),rejected)
    reject('M6_nonvacuous_field_bound',lambda:field_bound(dynamic,6,F(1,100000)),rejected)
    bound42=field_bound(dynamic,42,F(1,100000))
    reject('M42_automatic_particle_transfer',lambda:field_bound(dynamic,42,F(1,100000),particle_modulus=6),rejected)
    reject('dynamical_zero_tail_bound',lambda:field_bound(dynamic,42,F(1,100000),claimed=F(0)),rejected)
    classical=Source('frozen-reciprocal-two-mobile-literal-C-side2-96ports','classical_probability',9216,'side2_96ports','signed_pair512_scheduler')
    local12=Source('external_LOCAL12','position_times_12_slots',1152,'size2_logical_plus_retained_padding','W=S_R')
    reject('C10_to_classical_force',lambda:require_same_source(endpoint.source,classical),rejected)
    reject('C10_to_external_LOCAL12',lambda:require_same_source(endpoint.source,local12),rejected)
    reject('static_field_trajectory_clock',lambda:require_trajectory_clock(Endpoint(classical,'signed_pair',None,'static_smoothed_3D_Coulomb',(),True)),rejected)
    # Execute a changed-regulator coherent member on every lifted reservoir column.
    reference,_,_=chain.model.reference(); _,vectors=np.linalg.eigh(reference)
    columns=np.column_stack([k@vectors[:,-1] for _,k in chain.model.source.reservoir_instrument()])
    states=[(mask,dynamic.centered(f[3],6)) for mask,f in chain.model.states]
    radius=max(max(abs(x) for x in dynamic.integer_record(mask,k)) for mask,k in states)
    assert radius==5
    model42=dynamic.Family(42); lifted=model42.lift(columns,states); final=model42.evolve(lifted)
    norm_error=close(np.sum(abs(final)**2),np.sum(abs(columns)**2),'M42 complete lifted preparation norm')
    return {'M6_actual_H_error':h_error,'M6_actual_MIX_error':instrument_error,
        'M6_complete_clock_TV_upper':str(dynamic.bound(6,5)),
        'preparation_radius_all384_labels':radius,'M42_extra_terminal_MIX_TV_upper':dynamic.outward(dynamic.bound(42,radius,terminal_events=1)),'M42_dimension':model42.n,'M42_complete_clock_TV_upper_exact':str(bound42),
        'M42_complete_clock_TV_upper':dynamic.outward(bound42),'M42_coherent_lift_norm_error':norm_error,
        'M42_preparation':'All four M6 reservoir columns isometrically lifted; no M42 ground-state substitution.',
        'particle_instrument_regulator':6,'useful_bound_regulator':42,
        'M42_particle_transfer':False,'static_3D_refinement_trajectory_transfer':False,
        'closed_quantum_links':['actual C10 native instrument to declared common Poisson/H clock',
            'same C10 clock to complete operational regional detector and actual-H energy records',
            'actual M6 source/preparation/clock into periodic integer-field family with vacuous M6 bound',
            'lifted preparation and same supplied clock into M42 integer-field comparison with explicit bound'],
        'closed_classical_links':['signed-pair native current to integer field controller',
            'controller energy to reciprocal rate law with complete marked-path TV certificate',
            'same classical native force instrument to unconditional current score'],
        'missing_transfers':['Reference-preserving M42 proper-code detector and particle instrument under its actual H/reference, with bounded comparison to integer-field observables.',
            'A common source/preparation/clock map between the C10 quantum family and the reciprocal classical force construction.',
            'A dynamic spatial refinement theorem carrying complete records and force observables into a three-dimensional Maxwell/Lorentz limit.',
            'A native compilation/resource-time theorem for coherent H services, proper-code routing and spectral apparatus, or an explicit common architecture retaining these services.',
            'Physical species, calibrated time, masses beyond supplied anchors and full literal all31 observer criteria.']}


def audit_closure():
    entries={}
    for directory in ('audit_clock_force','audit_field','audit_particle'):
        path=CHAIN/directory/'AUDIT_CONTRACT.json'
        if not path.exists():
            raise ValueError('required independent audit contract absent: '+str(path))
        document=read(path)
        if directory=='audit_clock_force':
            if document.get('valid_member_scientific_claims_pass') is not True or document.get('scientific_blockers'):
                raise ValueError('clock/force finite-member scientific audit not passed')
            expected={'clock.general_density_input_validation','force.general_rate_horizon_input_validation'}
            failures={x['claim'] for x in document.get('findings',[]) if x['status']=='FAIL'}
            if failures != expected:
                raise ValueError('unclassified clock/force audit finding')
            admission='PASS_SCOPED; two raw admission gaps assigned to checked facade'
        elif directory=='audit_field':
            if document.get('supported_claims_pass') is not True:
                raise ValueError('field supported claims audit not passed')
            accepted={'fixed_point_tracking','static_smoothed_3D_Coulomb','actual_C10_M6_baseline',
                'interacting_integer_register_limit','extra_terminal_MIX','classical_fixedpoint_energy_rate_adapter'}
            rejected={'M6_small_full_clock_limit_error','M6_integer_chart_equals_periodic_electric_energy',
                'classical_force_to_C10_join','static_charge_lift_is_dynamic_refinement','Maxwell_mechanical_force_or_physical_seconds'}
            if {c['claim'] for c in document['claims'] if c['status']=='PASS'}!=accepted or {c['claim'] for c in document['claims'] if c['status']=='FAIL'}!=rejected:
                raise ValueError('unclassified field claim audit')
            admission='PASS_SCOPED; five unsupported transfers remain refused'
        else:
            if document.get('admission')!='PASS_FINITE_DECLARED_SOURCE_STATE_CLOCK_PARTICLE_INSTRUMENT' or not all(document['checks'].values()) or not all(document['mutation_rejections'].values()):
                raise ValueError('particle finite-member audit not passed')
            admission=document['admission']
        entries[directory]={'path':str(path.relative_to(ROOT)),'sha256':sha(path),
            'admission':admission,'contract':document}
    return entries


def row_map():
    source=ROOT/'openwave-tuning/next/ROW1_NATIVE/ALL31_READY.json'
    ready=read(source)
    support={1:['quantum_operational_regional_charge','same_C10_clock_particle'],
             2:['finite_actual_H_excitation_supplied_calibration'],
             3:['declared_common_source_clock'],4:['finite_continued_excitation_only'],
             8:['complete_charged_reservoir'],20:['classical_field_rate_force_finite','static_smoothed_3D_field'],
             21:['dynamical_integer_field_family_only'],29:['dynamical_integer_field_family_only'],
             30:['finite_quantum_H_semigroup'],31:['complete_actual_H_energy_groups']}
    rows=[]
    for row in ready['rows']:
        attachment=row.get('canonical_C10_attachment',{})
        rows.append({'id':row['id'],'criterion':row['criterion'],'literal_gate':row['literal_gate'],
            'official_status_at_snapshot':row['official_status'],
            'component_evidence':support.get(row['id'],[]),
            'required_same_source_connection':attachment.get('missing_same_source_connection','No complete literal chain supplied.'),
            'full_physical_criterion_satisfied_by_integration':False,
            'verified_full_physical_review_count':row.get('complete_physical_criterion_reviews_verified',0),
            'row1_finite_literal_mechanism':row['id']==1,
            'row1_law_control_scope':'ROW1_SCOPE.json' if row['id']==1 else None,
            'row1_physical_interpretation':'partial-at-most; no uniqueness or tape-only omission award' if row['id']==1 else None,
            'claim_gate':'finite_mechanism_and_interpretation_distinct' if row['id']==1 else 'physical_criterion_unestablished'})
    assert len(rows)==31 and {r['id'] for r in rows}==set(range(1,32))
    return {'schema':'oph.chain.per_row_gate.v1','snapshot':str(source.relative_to(ROOT)),
        'snapshot_sha256':sha(source),'rows':rows,'official_rows_changed':[],
        'all31_physical_complete':False,'common_physics_complete':False,
        'rule':'Component evidence and independent component audits do not establish a full literal physical row.'}


def main():
    import snapshot_guard
    if not snapshot_guard.ACTIVE:raise RuntimeError('Use run.py for the isolated historical source execution')
    before=pins()
    audits=audit_closure()
    objects,quantum=quantum_run()
    rejected,equivalence,guard_witnesses=guard_checks(objects)
    classical=classical_run(rejected,equivalence)
    transfers=transfer_checks(objects,rejected)
    rows=row_map()
    assert pins()==before
    result={'schema':'oph.chain.checked_integration.v1','status':'PASS_SCOPED_FINITE_COMPOSITIONS',
        'quantum_composition':quantum,'classical_composition':classical,'transfer_gates':transfers,
        'guard_mutations':rejected,'mutation_count':len(rejected),'numerical_domain_witnesses':guard_witnesses,
        'frozen_valid_input_comparison_max_errors':equivalence,
        'admission_audit_closure':{'clock_matrix_preconditions':'Checked facade before raw calls; original frozen helper bytes preserved.',
            'force_rate_horizon_preconditions':'Checked facade before raw calls; original frozen helper bytes preserved.',
            'state_numerical_tolerance':STATE_TOL,'number_offblock_tolerance':0.,'no_projection_clipping_or_renormalization':True,
            'raw_helpers_untrusted_input_safe':False,'checked_facade_mutation_audit_pass':True},
        'row1_scope':read(HERE/'ROW1_SCOPE.json'),'independent_audits':audits,'external_source_map':read(HERE/'EXTERNAL_SOURCE_MAP.json'),'input_pins':before,
        'snapshot_execution':{**read(HERE/'HISTORICAL_SOURCE.json'),'runtime_isolation':snapshot_guard.receipt()},
        'source_sha256':{p.name:sha(p) for p in (HERE/'integrate.py',HERE/'strict_facade.py',HERE/'gates.py',HERE/'seal.py',HERE/'run.py',HERE/'snapshot_guard.py',HERE/'build_snapshot.py')},
        'official_rows_changed':[],'all31_physical_complete':False,'common_physics_complete':False,
        'producer_files_modified':[]}
    (HERE/'ALL31_GATE.json').write_text(json.dumps(rows,indent=2)+'\n')
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    from seal import seal
    seal()
    print(json.dumps({'status':result['status'],'rejected_mutations':len(rejected),
        'frozen_comparator_max_error':max(equivalence.values()),'quantum_energy_coarsening_error':quantum['frozen_576_by_348_poststate_coarsening_error'],
        'all31_physical_complete':False,'common_physics_complete':False},indent=2),flush=True)


if __name__=='__main__':
    main()
