"""Particle/hole localization and coarse native records of the common walk."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

import mass_coin as mass
HERE=Path(__file__).resolve().parent
ROOT=mass.ROOT
sys.path[:0]=[str(HERE.parent/'manybody'),str(HERE.parent/'energy'),
              str(ROOT/'openwave-tuning/auditor/common_walk_particle')]
import manybody as mb
import energy
import instrument_contract as apparatus
from native_phase_probe import controlled_local_word
from model.certified_reference import certificate as cert


@lru_cache(None)
def reset_bound(blank,ready):
    """Uniform conditional lower bounds; no independence or exact product."""
    if any(type(n)is not int or n<0 for n in (blank,ready)):
        raise ValueError('nonnegative reset counts required')
    minima={}
    for code in (0,apparatus.pointer.READY_CODE):
        probabilities=apparatus.configuration.success_probabilities(
            4,32,full=True,mode='reset',target=apparatus.circuit.decode(code))
        minima[code]=min(map(F,probabilities))
    p0,pr=minima[0],minima[apparatus.pointer.READY_CODE]
    fail=blank*(1-p0)+ready*(1-pr)
    return {'blank_carriers':blank,'ready_carriers':ready,
            'blank_success_lower_exact':str(p0),'ready_success_lower_exact':str(pr),
            'timeout_upper_exact':str(fail),
            'semantics':'Worst-initial-state and worst-conditional success lower bounds composed by a union bound. No reset independence or exact product law is assumed. Actual first-failure probabilities use their history-dependent conditional factors; every fine failure is retained.'}


def point_preparation(vertex=0):
    if type(vertex)is not int or not 0<=vertex<96:raise ValueError('physical vertex required')
    state=np.zeros((128,16),complex)
    state[mb.common.encode_position(vertex)]=mb.compiled.controls()['uniform12_preparation'][:,0]
    return state


def compact_cycle():
    """A normalized oriented native triangle in moving event slots."""
    triangle=(0,1,5)
    neighbours,_=mb.common.geometry_tables()
    state=np.zeros((128,16),complex)
    for a,b in zip(triangle,triangle[1:]+triangle[:1]):
        d=next(d for d in range(6)if neighbours[a,d]==b)
        state[mb.common.encode_position(a),2*d+1]+=1/np.sqrt(6)
        state[mb.common.encode_position(b),2*d+1]-=1/np.sqrt(6)
    assert abs(np.vdot(state,state).real-1)<1e-14
    return state


def coarse_position_calls(source_wires=2048, undo=False):
    """Compute7 address parities, preserving event/field coherence at Q=1."""
    if type(source_wires)is not int or source_wires<2048:
        raise ValueError('fresh flag wires outside the full source required')
    calls=[(mode,source_wires+bit)for bit in range(7)for mode in range(2048)
           if (mode//16)>>bit&1]
    return list(reversed(calls))if undo else calls


def native_position_contract(field_bits=0,source_wires=2048):
    if field_bits not in (0,1152) or type(field_bits)is not int:
        raise ValueError('base matter or full288-edgeZ16 field required')
    minimum=2048+field_bits+(3 if field_bits else 0)
    if type(source_wires)is not int or source_wires<minimum:
        raise ValueError('fresh flags must follow all matter, field and source scratch wires')
    calls=coarse_position_calls(source_wires)
    assert len(calls)==7168
    # Complete one-particle columns and within-port off-diagonals: measuring
    # the seven flags identifies only the position, never its16 event labels.
    flags=[]
    for mode in range(2048):
        flag=sum(((mode//16)>>bit&1)<<bit for bit in range(7))
        assert flag==mode//16
        flags.append(flag)
    assert all(len(set(flags[16*x:16*x+16]))==1 for x in range(128))
    native=json.loads((mass.REVIEW/'SOURCE_ERROR.json').read_text())
    cx=F(native['primitive_operator_errors_exact']['cx'])
    pointer=apparatus.pointer_bound()
    reads=7+field_bits
    resets=reset_bound(7,reads)
    error=2*len(calls)*cx+reads*pointer+F(resets['timeout_upper_exact'])
    full_parities=[sum((mode//16)>>bit&1 for mode in range(2048)
                      if mb.common.decode_position(mode//16)is not None and mode%16<12)%2
                   for bit in range(7)]
    assert full_parities==[0]*7
    return {'position_flag_wires':list(range(source_wires,source_wires+7)),
            'declared_full_source_wires':source_wires,
            'compute_CNOT_calls':len(calls),'uncompute_CNOT_calls':len(calls),
            'native_pulses_excluding_resets':6*len(calls)+6*reads,
            'field_bits_measured':field_bits,'complete_native_reads':reads,
            'reset_contract':resets,'apparatus_TV_upper_exact':str(error),
            'apparatus_TV_upper':math.nextafter(float(error),math.inf),
            'instrument':'All128 coarse position outcomes P_x rho P_x on the one-particle2048-mode sector. Every event and retained field/reference coherence survives position read. Optional field-bit reads are an additional explicit dephasing intervention.',
            'outside_one_particle_sector':'The full reversible circuit computes the seven occupation parities. No many-particle position interpretation is inferred.',
            'one_hole_domain':'The fully filled1152physical-mode reference has zero address parities; omitting one mode gives its hole position, with within-port hole-event coherence retained.',
            'program':'coarse_position_calls(full_source_wires), native C/E pointer on7 flags, reverse calls. If configured, read each field bit with a fresh native pointer. Every reset/export/fine-pointer failure retained.'}


def massive_phase_contract(h,mass_parameter,steps):
    h,mass_parameter=mass.parameters(h,mass_parameter)
    if type(steps)is not int or not 0<=steps<=32:raise ValueError('bounded phase program required')
    built=mb.compiled.controls()
    counts=Counter((e['kind'],e['angle'])for e in mass.coin_word(h,mass_parameter))
    for kind,uses in Counter(k for k,t in built['source_program']).items():
        for event,number in Counter((e['kind'],e['angle'])for e in built['primitive_words'][kind]).items():
            counts[event]+=uses*number
    wrapper=F();pulses=0;table={}
    for (kind,angle),uses in counts.items():
        word=controlled_local_word(kind,angle)
        end,fixed,variables=cert.consume(word,0,apparatus.controlled_template(kind,angle))
        assert end==len(word)and not variables
        wrapper+=uses*(fixed+cert.controls_bound(word,table));pulses+=uses*len(word)
    source=mass.native_error(h,mass_parameter)
    old=json.loads((ROOT/'openwave-tuning/auditor/common_walk_particle/INSTRUMENT_CONTRACT.json').read_text())
    initial=json.loads((mass.REVIEW/'SOURCE_ERROR.json').read_text())
    reset=reset_bound(24,1)
    total=(F(initial['whole_source_preparation_error_exact'])
           +steps*(F(source['whole_step_operator_error_exact'])+wrapper)
           +F(old['phase']['extra_preparation_analysis_pointer_error_exact'])
           +F(reset['timeout_upper_exact']))
    return {'steps':steps,'controlled_massive_step_native_pulses':pulses,
            'controlled_wrapper_error_exact':str(wrapper),'complete_TV_error_exact':str(total),
            'complete_TV_error_upper':math.nextafter(float(total),math.inf),'reset_contract':reset,
            'source_word':'same frozen swap calls preceded by mass_coin.coin_word(h,m); every scalar phase controlled',
            'scope':'Base encoded one-particle source. Fock/field implementation requires control of its entire source word and fresh wires, not this detector by relabeling.'}


def run():
    paths=[Path(__file__),Path(mass.__file__),Path(mb.__file__),Path(energy.__file__),
           HERE.parent/'manybody/MANIFEST.json',HERE.parent/'manybody/HOLE_SECTOR.json',
           HERE.parent/'energy/MANIFEST.json',mass.REVIEW/'SOURCE_ERROR.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    seeds={'point_uniform_event':point_preparation(),'compact_triangle':compact_cycle()}
    dark=np.zeros((128,16),complex);dark[0,0]=1/np.sqrt(2);dark[0,2]=-1/np.sqrt(2)
    seeds['held_dark']=dark
    experiments=[]
    for name,initial in seeds.items():
        for mass_parameter in (F(0),F(1)):
            h=F(1,8);state=initial.copy();rows=[]
            for n in range(17):
                if n in (0,1,2,4,8,16):
                    record=mb.common.complete_position_records(state)
                    correlation=np.vdot(initial,state)
                    rows.append({'step':n,'source_return_correlation':[float(correlation.real),float(correlation.imag)],
                                 'all96_position_probabilities':record['valid_position_probabilities'].tolist(),
                                 'position0_probability':float(record['valid_position_probabilities'][0]),
                                 'source_padding':record['source_padding_probability'],
                                 'event_padding':record['event_padding_probability_on_valid_source'],
                                 'norm_squared':record['original_norm_squared']})
                    assert abs(record['original_norm_squared']-1)<1e-10
                    if name=='compact_triangle':assert np.linalg.norm(state-np.exp(1j*n*float(h*mass_parameter))*initial)<1e-10
                    if name=='held_dark':assert np.linalg.norm(state-(-np.exp(1j*float(h*mass_parameter)))**n*initial)<1e-10
                if n<16:state=mass.step(state,h,mass_parameter)
            experiments.append({'preparation':name,'h':str(h),'mass_parameter':str(mass_parameter),'records':rows})
    # Exact free one-particle and one-hole CAR propagators use the same U.
    # The sign of physical current reverses; supplied reference changes.
    state=point_preparation()
    for _ in range(4):state=mb.common.apply_walk(state)
    hole=point_preparation().conj()
    for _ in range(4):hole=mb.common.apply_walk(hole.conj()).conj()
    assert np.max(abs(state-hole))<1e-11
    massive_hole=point_preparation()
    for _ in range(4):massive_hole=mass.step(massive_hole.conj(),F(1,8),F(1)).conj()
    opposite_mu=point_preparation()
    for _ in range(4):opposite_mu=mass.step(opposite_mu,F(1,8),F(-1))
    assert np.max(abs(massive_hole-opposite_mu))<1e-10
    mu=1/8
    determinant=np.exp(-1j*mu)*np.exp(193j*mu)*(-np.exp(1j*mu))**768
    assert abs(determinant-np.exp(960j*mu))<2e-12
    # A finite-volume mode spacing vanishes as size increases. The actual
    # first acoustic frequency is compared on development sizes2,3,4 only.
    spacings=[]
    for L in (2,3,4):
        k=2*np.pi/L
        smallest=min(float(np.linalg.eigvalsh(mass.spatial.symbol((k,0,0)))[0]),
                     float(np.linalg.eigvalsh(mass.spatial.symbol((0,k,0)))[0]))
        spacings.append({'development_size':L,'positive_acoustic_phase':math.acos(1-smallest/12),
                         'rest_phase_at_zero_momentum':0.})
    assert spacings[-1]['positive_acoustic_phase']<spacings[0]['positive_acoustic_phase']
    # Native QPE measures the exact compact modes at0 andpi in baseW.
    qpe=[]
    for name in ('compact_triangle','held_dark'):
        output,probability=energy.phase_estimation(seeds[name],3)
        target=0 if name=='compact_triangle' else 4
        assert probability[target]>1-1e-10
        qpe.append({'preparation':name,'all8_probabilities':probability.tolist(),
                    'complete_QPE_TV_bound':energy.qpe_contract(3)['total_variation_bound']})
    result={'schema':'oph.local12.particle_records.v1',
        'checks':{'same_W_and_mass_amendment_only':True,'all_padded_position_laws_retained':True,
                  'compact_cycle_persistence_exact':True,'held_dark_mode_control':True,
                  'particle_hole_same_free_propagator':True,'finite_volume_gap_not_rest_mass':True,
                  'native_phase_measurements_of_compact_modes':True,
                  'coarse_mode_detector_preserves_event_coherence':True,
                  'massive_phase_wrapper_source_bound':True,'source_bytes_unchanged':True},
        'experiments':experiments,'finite_volume_controls':spacings,'native_QPE':qpe,
        'native_coarse_position':native_position_contract(),
        'native_joint_position_full_Z16_field':native_position_contract(288*4,3203),
        'native_massive_phase':massive_phase_contract(F(1,8),F(1),8),
        'CAR_correlation_identity':'<empty|c_f Gamma(U)^n c_g^dagger|empty>=f*U^n g. On the filled physical reference, the one-hole propagator is det(U)^n conjugate(U)^n; baseW det=1 and real.',
        'hole_boundary':'A separately supplied filled reference is used. No simultaneous opposite species, pair creation, annihilation or energy-ground preparation is inferred.',
        'massive_hole_reference':'Physical one-hole evolution is det(U_mu)*conjugate(U_mu), withdet(U_mu)=exp(i960mu) on1152physical modes. Relative to that same filled reference, the hole matrix isU_-mu. Its global determinant phase must remain in a controlled full-source experiment.',
        'persistence_boundary':'Exact compact charged modes persist in the free source. A triangle cycle haszero basequasienergy; the held-dark mode hasphasepi andzerovelocity. Neither alone establishes a mobile massive particle or Lorentz inertia. Gauge dressing/field coupling must be checked on their actual source.',
        'pole_boundary':'Finite unitary correlators have unit-circle poles and recurrences. This fact alone does not establish an infinite-volume particle pole, decay width or experimental lifetime.',
        'source_sha256':pins,'official_rows_changed':[]}
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    (HERE/'PARTICLE_RECORDS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'joint_position_field_detector_TV':result['native_joint_position_full_Z16_field']['apparatus_TV_upper'],
                      'massive_phase_TV':result['native_massive_phase']['complete_TV_error_upper']}))


if __name__=='__main__':run()
