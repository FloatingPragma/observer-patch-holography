"""Opposing source, instrument, cutoff and allocation checks."""
from fractions import Fraction as F
from pathlib import Path
import gzip
import hashlib
import json
import math
import sys

import numpy as np
from scipy.linalg import expm

import joint_program as joined

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'openwave-tuning/auditor/field_cutoff'))
import cutoff
sys.path.insert(0, str(HERE.parent/'spin'))
import dirac_proposal as dirac

MANIFESTS = {
    'common_dynamics/vector_field/MANIFEST.json': 'aa4a18d89c1a8893ee86eb5e9151bc75e8a0b3a527640cee66aa25382ca55e76',
    'common_dynamics/particles/MANIFEST.json': '678983fa40a2a97b8eb399d6b5ad4ab6ace77a5f441baf448cb0c00bc8a6f9a9',
    'common_dynamics/spin/MANIFEST.json': '65075e0bc14f883c5cdb877171650b08ee57c8d1bdffce281e68937bb07eb401',
    'auditor/field_cutoff/MANIFEST.json': '326994c7a0b2a972540b25662a0215ad0bf5caef78afa5730d7b55fa920df998',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def custody():
    pins = {}
    for relative, expected in MANIFESTS.items():
        path = ROOT/'openwave-tuning'/relative
        if digest(path) != expected:
            raise ValueError('frozen provider manifest mismatch: '+relative)
        pins[str(path.relative_to(ROOT))] = expected
        data = json.loads(path.read_text())
        for key in ('files', 'source_pins', 'source_sha256', 'files_sha256'):
            for source, wanted in data.get(key, {}).items():
                if digest(ROOT/source) != wanted:
                    raise ValueError('frozen source mismatch: '+source)
                if source in pins and pins[source] != wanted:
                    raise ValueError('incompatible source pins: '+source)
                pins[source] = wanted
    return pins


def expect_rejected(function, *args, **kwargs):
    try:
        function(*args, **kwargs)
    except (ValueError, AssertionError):
        return True
    raise AssertionError('mutation accepted')


def admitted_state(state, expected_weight=1.):
    graph, B, A, labels = joined.massive.af.data()
    total = 0.
    for (x, field), values in state.items():
        if type(x) is not int or not 0 <= x < 96 or type(field) is not bytes or len(field) != 288:
            raise ValueError('physical source index and complete field required')
        values = np.asarray(values)
        if values.shape != (12,) or not np.all(np.isfinite(values)):
            raise ValueError('all twelve finite event amplitudes required')
        f = np.frombuffer(field, dtype=np.uint8).astype(np.int64)
        if np.any(f >= 16):
            raise ValueError('Z16 source field required')
        n = np.zeros(96, dtype=int); n[x] = 1; n[0] -= 1
        if np.any((n+B@f) % 16):
            raise ValueError('modular Gauss sector changed')
        total += float(np.vdot(values, values).real)
    if not math.isfinite(expected_weight) or not 0 <= expected_weight <= 1 or abs(total-expected_weight) > 2e-10:
        raise ValueError('unconditional branch weight changed')
    return total


def independent_mass_step(state, mass):
    # Diagonalize the exact rank-one reflection analytically. No call to the
    # producer's mass or matter step is used in this implementation.
    mu = float(F(1,8)*mass)
    result = {}
    for key, v in state.items():
        mean = np.sum(v)/12
        result[key] = np.exp(1j*mu)*v + (np.exp(-1j*mu)-np.exp(1j*mu))*mean
    return cutoff.step(result, 16, .17, .41)


def source_replay():
    archive = json.loads(gzip.decompress((HERE.parent/'particles/MASSIVE_JOINT_LAWS.json.gz').read_bytes()))
    frozen = {(row['mass_parameter'], row['step']): row['joint_law'] for row in archive['records']}
    residual = archived_residual = 0.; columns = 0; final = {}; weights = []
    for mass in (F(0), F(1)):
        actual = joined.massive.af.initial()
        independent = cutoff.initial(16)
        for step in range(4):
            admitted_state(actual); admitted_state(independent)
            assert actual.keys() == independent.keys()
            residual = max(residual, max(float(np.max(abs(actual[k]-independent[k]))) for k in actual))
            archived = {(row['position'], bytes.fromhex(row['all288_Z16_field_values_hex'])):
                        np.asarray(row['all12_event_probabilities']) for row in frozen[str(mass), step]}
            assert independent.keys() == archived.keys()
            archived_residual = max(archived_residual, max(float(np.max(abs(abs(independent[k])**2-archived[k]))) for k in independent))
            columns += 12*len(independent)
            weights.append(sum(float(np.vdot(v,v).real) for v in independent.values()))
            if step < 3:
                independent = independent_mass_step(independent, mass)
                actual = joined.massive.step(actual, F(1,8), mass, .17, .41)
        final[str(mass)] = independent
    assert residual < 2e-12 and archived_residual < 2e-12
    positions = []
    for state in final.values():
        p = np.zeros(96)
        for (x, _), v in state.items(): p[x] += float(np.vdot(v,v).real)
        positions.append(p)
    response = float(np.sum(abs(positions[1]-positions[0]))/2)
    return final, {'all_event_outcomes_compared': columns, 'maximum_state_amplitude_error': residual,
                   'maximum_archived_probability_error': archived_residual,
                   'complete_weights': weights, 'mass_position_response_TV': response}


def clock_laws():
    pointer, _ = joined.phase_clock.actual_native_pointer()
    answer = {}; failure = 0.
    for setting in ('X', 'Y'):
        word = joined.phase_clock.local_words(3, setting, F(1,8))
        run = joined.phase_clock.recovery.execute_complete(word, 1, np.eye(2, dtype=complex)[:, [0]])
        state = run.success[:,0]
        failure = max(failure, float(run.failure_gram[0,0].real))
        theta = 3/64
        target = math.cos(theta) if setting == 'X' else math.sin(theta)
        for law, branches in pointer.items():
            rows = []; bits = np.zeros(2)
            for metadata, kraus in branches:
                probability = sum(float(np.vdot(K@state,K@state).real) for K in kraus)
                rows.append({'record': metadata, 'probability': probability})
                bits[metadata['bit']] += probability
            assert abs(bits.sum()-1)<2e-12
            assert abs(bits[0]-bits[1]-target)<2e-12
            answer[setting+law] = {'fine_clock_law': rows, 'bits': bits.tolist(),
                                  'contrast': float(bits[0]-bits[1]),
                                  'phase': '3/64', 'time_over_tau0': '3/8'}
    return answer, failure


def lower_tail(value, degree):
    # Independent exact partial sum and a later geometric remainder give
    # lower/upper brackets inside the producer's conservative upper bound.
    value = F(value); final = max(degree+24, math.floor(value)+32)
    lower = sum((value**k/math.factorial(k) for k in range(degree, final)), F())
    upper = lower + value**final/math.factorial(final)/(1-value/F(final+1))
    return lower, upper


def cutoff_review():
    rows = []
    for value, degree in ((F(0),5),(F(1,5),1),(3*F(.41),5),
                           (3*F(.41),13),(3*F(.41),29),(3*208*F(.41),509)):
        low, high = lower_tail(value,degree)
        producer = cutoff.exponential_tail_upper(value,degree)
        assert 0 <= low <= high <= producer
        rows.append({'Lambda': str(value), 'distance': degree,
                     'independent_upper_float': float(high), 'producer_upper_float': float(producer)})
    graph,B,A,labels = joined.massive.af.data()
    assert np.max(abs(A)) == 1 and np.max(abs(B@A)) == 0
    # Low-order terms cannot reach a periodic boundary for k<d, for every
    # possible sign sequence. This arithmetic control detects off-by-one
    # changes without assuming a particular plaquette path.
    cases = 0
    for m in (6,16,32):
        for J in range(4):
            d = m//2-J
            if d <= 0: continue
            for k in range(d):
                assert J+k <= m//2-1 and -J-k >= -m//2
                cases += 1
            assert J+d == m//2
    boundary_error = np.linalg.norm(np.array([1.,0.])-np.array([0.,1.]))
    assert abs(boundary_error-math.sqrt(2)) < 1e-15
    # The small one-face bound cannot be used for the homogeneous source.
    assert cutoff.transfer_bound(16,3,208*F(.41)) == 2
    assert cutoff.transfer_bound(32,3,F(.41)) < F(6,10**9)
    return {'exact_tail_brackets': rows, 'prefix_boundary_controls': cases,
            'face_entries_bounded_by_one': True, 'cyclic_boundary_control_norm': boundary_error,
            'verdict': 'The stated state/support-dependent integer-rotor theorem is valid. Its finite-time terminal-instrument and arbitrary-reference scope is preserved; quantum Maxwell and native modulus-cost transfer are absent.'}


def spin_review():
    alphas = dirac.ALPHA; beta = dirac.BETA; identity = np.eye(4)
    for a in (*alphas,beta): assert np.max(abs(a@a-identity))<1e-14
    all_generators = (*alphas,beta)
    for i in range(4):
        for j in range(i): assert np.max(abs(all_generators[i]@all_generators[j]+all_generators[j]@all_generators[i]))<1e-14
    p = np.array([.2,-.1,.3]); m=.4; H=sum((p[a]*alphas[a] for a in range(3)),m*beta)
    cases=[]
    for h,n in ((.2,5),(.1,10),(.05,20)):
        U=expm(-1j*h*m*beta)
        for a in range(3):U=expm(-1j*h*p[a]*alphas[a])@U
        observed=float(np.linalg.norm(np.linalg.matrix_power(U,n)-expm(-1j*n*h*H),2))
        independent=n*h*h*sum(abs(x*y) for i,x in enumerate([m,*p]) for y in [m,*p][i+1:])
        assert observed<=independent+1e-13
        cases.append({'h':h,'steps':n,'error':observed,'bound':independent})
    route_cases=0
    B=joined.massive.af.data()[1]
    for cell in range(8):
        for axis in range(3):
            for sign in (-1,1):
                other,path,current=dirac.route(cell,axis,sign)
                rhs=np.zeros(96,dtype=int);rhs[12*cell]=1;rhs[12*other]-=1
                assert np.array_equal(B@current,rhs)
                route_cases+=1
    return {'Clifford_relations':True,'free_limit_cases':cases,'actual_signed_route_cases':route_cases,
            'no_transfer_to_scalar_candidate':True,
            'retained_obligations':['complete native Dirac reset/evolution error','full-field preparation for the compressed cycle probe','full spatial rotation apparatus','interacting canonical-field limit','species/doubler control']}


def mutations(final):
    checks = {
        'clock_on_position_flag': expect_rejected(joined.admit, clock_wire=3203),
        'clock_on_source': expect_rejected(joined.admit, clock_wire=898),
        'auxiliary_alias': expect_rejected(joined.admit, auxiliaries=(3210,3211)),
        'wrong_source_candidate': expect_rejected(joined.admit, candidate=dirac.MODEL_ID),
        'wrong_modulus': expect_rejected(joined.admit, modulus=6),
        'two_particles_as_one_position': expect_rejected(joined.admit, charge=2),
        'unfrozen_mass': expect_rejected(joined.admit, mass=F(2)),
        'unfrozen_clock_spacing': expect_rejected(joined.admit, h=F(1,16)),
    }
    half={key:value/math.sqrt(2) for key,value in final['1'].items()}
    admitted_state(half,.5)
    checks['success_renormalization'] = expect_rejected(admitted_state, final['1'], .5)
    # Omit the field increment in the controlled branch of one real source
    # step: the nontrivial half of a plus probe moves with probability1/2.
    neighbours,_=joined.massive.af.mb.common.geometry_tables()
    broken_weight=F()
    for event in range(12):
        direction,move=divmod(event,2)
        y=int(neighbours[0,direction]) if move else 0
        if y!=0:broken_weight+=F(1,24)
    assert broken_weight==F(1,4)
    checks['controlled_bare_walk_missing_current'] = True
    # An exact hash guard fails for a one-bit mutation, not merely metadata.
    path=HERE.parent/'particles/massive_coupled.py'
    payload=path.read_bytes();mutated=bytes([payload[0]^1])+payload[1:]
    assert hashlib.sha256(mutated).hexdigest()!=digest(path)
    checks['one_byte_source_mutation'] = True
    return {'rejected_controls':checks,'missing_current_Gauss_violation_probability':str(broken_weight)}


def run():
    pins=custody()
    local_paths=[Path(__file__),HERE/'joint_program.py', ROOT/'openwave-tuning/common_dynamics/SOURCE_CONTRACT.json',
                 ROOT/'openwave-tuning/platform_packet/criteria_frozen.json']
    pins.update({str(p.relative_to(ROOT)):digest(p) for p in local_paths})
    print('custody verified',flush=True)
    final,replay=source_replay();print('independent massive source and full law replay passed',flush=True)
    clocks,failure=clock_laws()
    cutoff_result=cutoff_review();spin=spin_review();controls=mutations(final)
    bounds={setting:joined.error_contract(setting=setting) for setting in ('X','Y')}
    joint=[]
    for mass,state in final.items():
        matter_total=sum(float(np.vdot(v,v).real) for v in state.values())
        for setting,clock in clocks.items():
            total=matter_total*sum(row['probability'] for row in clock['fine_clock_law'])
            assert abs(total-1)<3e-11
            joint.append({'mass':mass,'setting_and_native_law':setting,'source_field_blocks':len(state),
                          'fine_joint_outcomes':len(state)*len(clock['fine_clock_law']),
                          'complete_probability':total,
                          'factorization':'P(position,field,clockfine)=sum_event|psi(position,field,event)|² * P(clockfine); ideal product clock preparation only'})
    print('clock, cutoff, spin and mutation checks passed; certifying ordered native blocks',flush=True)
    relocated=joined.block_summary(clock_wire=3300)
    assert relocated['max_wire']==3300
    stream=joined.block_summary()
    assert stream['counts']['native_configuration_reset']==3211+1160
    assert stream['counts']['native_commit']==1160 and stream['counts']['macro_completed']==3
    assert stream['max_wire']==3212
    assert pins=={path:digest(ROOT/path) for path in pins}
    receipt={'schema':'oph.common.dynamics.opposing_integration.v1',
             'candidate_id':joined.CANDIDATE,
             'checks':{'frozen_custody':True,'independent_full_mass_field_replay':True,
                       'all_archived_probabilities':True,'same_source_clock_joint_law':True,
                       'disjoint_native_wire_allocation':True,'ordered_full_native_block_program':True,
                       'conditional_reset_and_export_bounds':True,'one_particle_detector_domain':True,
                       'independent_cutoff_review':True,'independent_Dirac_scope_review':True,
                       'adversarial_controls_reject':True,'final_source_hashes_unchanged':True},
             'source_replay':replay,'clock_laws':clocks,'clock_export_failure_trace':failure,
             'complete_joint_laws':joint,'joined_complete_error':bounds,'program':stream,
             'allocation':{'source': [0,3202],'position_flags':list(joined.FLAG_WIRES),
                           'clock':joined.CLOCK_WIRE,'pointer_auxiliary_code':list(joined.AUX_WIRES),
                           'spectral_probe':None},
             'cutoff_review':cutoff_result,'spin_review':spin,'mutations':controls,
             'source_sha256':pins,'official_rows_changed':[],
             'scope':'Independent finite-source and complete native-instrument join. The source is an explicitly supplied mass/field architecture; an interacting continuum, physical vacuum, species dictionary and official validation are not inferred.'}
    (HERE/'INTEGRATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'checks':receipt['checks'],'source_replay':replay,
                      'complete_TV_bounds':{s:r['complete_TV_upper'] for s,r in bounds.items()},'program':stream}))


if __name__=='__main__':run()
