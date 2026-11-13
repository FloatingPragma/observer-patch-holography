"""A clock instrument conjugated by the actual electric-link walk.

The canonical source/string dressing and its complement are explicit. The
tracking word reverses the entire source/tag/link/event evolution before
decoding the source clock. It does not reset or discard cycle flux.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import lcm
from hashlib import sha256
import json
import base64
import zlib
import numpy as np

from model.native_gauss_walk import source as joint
from model.openwave_native_coherent_clock import source as clock
from model.common_source.record_code import MATCHINGS
from model.native_observer_program.contract import descriptor_sha256

PAIRS=tuple(MATCHINGS[0]+MATCHINGS[1])


def basis_wave(pair):
    return joint.Wave({(joint.canonical_config(pair,1),e):joint.ONE for e in range(60)},1,1,1)


@lru_cache(maxsize=3)
def transported_columns(n):
    if n not in (0,1,2):raise ValueError('declared complete zero through two-step frame required')
    waves=[];inverse_checks=[]
    for pair in PAIRS:
        initial=basis_wave(pair);value=initial
        for _ in range(n):value=joint.step(value,electric_phase=True)
        inverse=value
        for _ in range(n):inverse=joint.step(inverse,inverse=True,electric_phase=True)
        if not joint.equivalent(inverse,initial):
            raise ArithmeticError('full source/tag/link/event inverse did not restore its basis column')
        if inverse.probability()!=F(1,12):raise ArithmeticError('unnormalized operator basis changed')
        inverse_checks.append(True);waves.append(value)
    keys=sorted(set().union(*(w.amplitudes for w in waves)))
    common=lcm(*(w.denominator for w in waves))
    index={key:i for i,key in enumerate(keys)}
    a=np.zeros((len(keys),12),dtype=np.int64);b=np.zeros_like(a)
    for j,w in enumerate(waves):
        for key,(x,y)in w.amplitudes.items():
            a[index[key],j]=x*(common//w.denominator)
            b[index[key],j]=y*(common//w.denominator)
    bound=max(int(abs(a).max()),int(abs(b).max()))
    if 3*len(keys)*bound*bound>=2**62:
        raise ArithmeticError('exact integer Gram accumulator bound exceeded')
    denominator=60*common*common
    real=a.T@a+b.T@a+b.T@b
    imag=a.T@b-b.T@a
    if not np.array_equal(real,denominator*np.eye(12,dtype=np.int64))or np.any(imag):
        raise ArithmeticError('complete joint reference matrix-unit Gram is not identity')
    # Physical normalized source columns have event amplitude1/sqrt60.
    matrix=(a+b/2+1j*np.sqrt(3)*b/2)/np.sqrt(denominator)
    return keys,a,b,denominator,matrix,inverse_checks


def complete_joint_observation(n,rho):
    keys,a,b,denominator,matrix,_=transported_columns(n)
    if np.asarray(rho).shape!=(12,12):raise ValueError('complete admitted source density required')
    probabilities=np.einsum('ir,rs,is->i',matrix,rho,matrix.conj(),optimize=True).real
    if probabilities.min()<-2e-12 or abs(probabilities.sum()-np.trace(rho).real)>2e-10:
        raise ArithmeticError('joint source/link observation lost weight')
    source=defaultdict(float);gauss=defaultdict(float);event=np.zeros(64);energy=0.;dressing=0.
    for(key,p)in zip(keys,probabilities,strict=True):
        c,e=key;s=joint.state(c,1)
        source[c[0]]+=float(p);gauss[joint.electric.gauss(s)]+=float(p);event[e]+=p
        energy+=float(joint.electric.electric_energy(s))*p
        if c==joint.canonical_config(c[0],1):dressing+=p
    if set(gauss)!={(0,)*12}:raise ArithmeticError('the complete clock state violates Gauss law')
    payload=[[list(c[0]),c[1],list(c[2]),e,float(p)]for((c,e),p)in zip(keys,probabilities,strict=True)]
    raw=json.dumps(payload,separators=(',',':')).encode()
    return {'frame_index':n,'full_joint_nonzero_carrier_rows':len(keys),
        'complete_probability_vector_sha256':sha256(raw).hexdigest(),
        'complete_joint_table_encoding':'zlib+base64 JSON rows[pair,tag,flux30,event,diagnostic_probability]; source memory has its explicitly retained character distribution',
        'complete_joint_table_zlib_base64':base64.b64encode(zlib.compress(raw,9)).decode(),
        'source_pair_probabilities_diagnostic':[{'pair':list(pair),'probability':source[pair]}for pair in joint.electric.PAIRS],
        'all64_event_probabilities_diagnostic':event.tolist(),
        'full_probability_sum_diagnostic':float(probabilities.sum()),
        'complete_Gauss_zero_probability_diagnostic':gauss[(0,)*12],
        'electric_energy_mean_diagnostic':float(energy),
        'canonical_dressing_population_diagnostic':float(dressing),
        'cycle_flux_or_dressing_complement_discarded':False,
        'normalization_performed':False}


def frame_certificate():
    rows=[]
    for n in (0,1,2):
        keys,a,b,denominator,_,inverses=transported_columns(n)
        rows.append({'frame_index':n,'prepared_source_columns':12,'complete_matrix_units':144,
            'nonzero_complete_source_tag_link_event_rows':len(keys),
            'integer_Gram_denominator':denominator,
            'exact_complete_Gram_identity':True,'full_inverse_basis_checks':len(inverses),
            'canonical_dressing_after_inverse':'all twelve columns restored exactly; its full complement has zero amplitude on these inputs',
            'integer_amplitude_columns_sha256':sha256(a.tobytes()+b.tobytes()).hexdigest()})
    return rows


def run():
    certificates=frame_certificate()
    time,steps=F(64),4096
    bound=F(clock.rc.error_bound('C',time,steps)['trace_norm_upper_exact'])
    initial=np.outer(clock.pair_vectors(0,'X')[0],clock.pair_vectors(0,'X')[0].conj())
    first=clock.evolve(initial,'C',time,steps)
    records=[]
    for axis in ('X','Y'):
        effects=clock.effects(0,axis);probs=np.zeros((3,3));after1=np.zeros((36,36),complex)
        for i,p in enumerate(effects):
            branch=p@first@p;after1+=branch
            second=clock.evolve(branch,'C',time,steps)
            for j,q in enumerate(effects):probs[i,j]=np.trace(q@second).real
        second=clock.evolve(after1,'C',time,steps)
        # The actual inherited transfer is block diagonal in matching labels.
        for rho in (first,second):
            outside=rho.copy();outside[:12,:12]=0
            if np.any(outside):raise ArithmeticError('source clock left the admitted twelve matching rows')
        if abs(probs.sum()-1)>2e-10 or probs.min()<-2e-12:
            raise ArithmeticError('complete two-time native clock records failed')
        intervals=clock.interval_joint(F(1,256),time,axis,min(F(2),2*bound))
        if not all(float(intervals[i][j][0])-2e-11<=probs[i,j]<=float(intervals[i][j][1])+2e-11 for i in range(3)for j in range(3)):
            raise ArithmeticError('native clock bound does not cover its dressed instrument')
        records.append({'axis':axis,'law':'C','spectator':0,'integer_charge_magnitude':1,'active_memory_character':1,
            'frame_indices':[1,2],'time_per_clock_block_exact':str(time),'native_steps_per_clock_block':steps,
            'complete_two_time_native_probabilities_diagnostic':probs.tolist(),
            'complete_two_time_probability_intervals_exact':[[list(map(str,z))for z in row]for row in intervals],
            'alternative_source_field_measurements_before_quadrature':[complete_joint_observation(1,first[:12,:12]),complete_joint_observation(2,second[:12,:12])],
            'alternative_measurements_inserted_into_clock_protocol':False,
            'inherited_two_block_trace_norm_bound_exact':str(min(F(2),2*bound)),
            'full_joint_walk_or_inverse_calls':8,
            'dressing_and_primitive_implementation_errors_additional':True})
    untracked=joint.prepare()
    for _ in range(2):untracked=joint.step(untracked,electric_phase=True)
    diagonal=np.diag(np.diag(initial))
    which_path=clock.evolve(diagonal,'C',time,steps)
    plus,minus,_=clock.effects(0,'X')
    which_path_contrast=float(np.trace((plus-minus)@which_path).real)
    if abs(which_path_contrast)>2e-11:
        raise ArithmeticError('a diagonal source preparation spontaneously supplied clock coherence')
    return {'schema':'oph.native_joint_gauss_clock.v1',
        'verdict':'PASS_SAME_SOURCE_LINK_CLOCK_INSTRUMENT_WITH_EXPLICIT_TRACKING',
        'common_program_sha256':descriptor_sha256(),
        'integer_extension_descriptor':joint.electric.extension_descriptor(),
        'source':'native_Gauss_walk with q1,character1,linkwindow4,canonical shortest strings and electric phase enabled',
        'clock_service':'D Phi_C|V D* on the invariant twelve-row support V of matching words0/1, plus its complete complementary outcome, conjugated by full joint W_G^n; inverse active-memory feedback and source recovery are included.',
        'dressing_service_is_declared_architecture':True,
        'dressing_interface':'D in this executed clock maps the twelve admitted source basis pairs in V to canonical tagged strings and uniform event states. D* has the entire complementary outcome, including every other source pair, string and event mode. Its implementing unitary/processor and precision are additional.',
        'exact_transported_matrix_units':certificates,'same_joint_two_time_records':records,
        'untracked_dressing_control':joint.decode_canonical_dressing(untracked),
        'which_path_preparation_control':{'native_X_contrast_diagnostic':which_path_contrast,'source_position_coherences_removed':True,'initial_code_error_bound_applied_to_this_control':False},
        'raw_source_recovery_alone_removes_cycle_flux':False,
        'bound_hypothesis':'Initial coherent code in the span of matching words0/1, canonically dressed, with arbitrary untouched references; not arbitrary initial flux, other matching words or syndrome states.',
        'measurement_choice':'Source/tag/link/event distributions and tracked clock quadratures are alternative instruments on the same evolved marginal. Their noncommuting readouts are not assigned a simultaneous joint distribution. The repeated clock experiment includes its actual first-quadrature backaction.',
        'scope':{'same_charged_source_field_and_clock':True,
            'physical_rest_mass_identified':False,'official_OpenWave_awards':0,
            'all_source_field_modes_projected_back_to_strings':False,
            'full_inverse_then_declared_dressing_service_used':True,
            'cycle_flux_reset_used':False,'local_autonomous_clock_proved':False,
            'native_dressing_decoder_compiled_by_this_package':False,
            'shared_stationary_particle_spectrum_proved':False,
            'physical_electron_or_GeV_identification':False,
            'literal_OpenWave_row3_complete':False,'new_Nature_validation':False}}
