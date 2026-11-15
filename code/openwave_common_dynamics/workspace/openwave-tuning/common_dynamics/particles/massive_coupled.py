"""The local mass coin and autonomous field on one declared source.

The mass coefficient and both field couplings are architecture inputs. This
finite experiment derives its records from the combined native program.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

import mass_coin
import particle_records
import coupled_response as response

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'vector_field'))
import autonomous as af
import composition
from model.observer.source_scalar import compiler

MODEL_ID='native-local12-mass-fermion-plaquette-v3-proposed'


def step(state,h,mass,alpha,beta):
    """One-particle restriction; no filled-hole determinant phase omitted."""
    h,mass=mass_coin.parameters(h,mass)
    mu=float(h*mass)
    R=np.ones((12,12))/6-np.eye(12)
    kick=math.cos(mu)*np.eye(12)-1j*math.sin(mu)*R
    prepared={key:kick@value for key,value in state.items()}
    return af.step(prepared,alpha,beta,charge=1)


def fock_kick_word(h,mass):
    h,mass=mass_coin.parameters(h,mass)
    mu=h*mass
    reflection=af.mb.fock_reflection_word()
    length=(len(reflection)-30)//2
    U=reflection[-length:]
    assert reflection[:length]==compiler.inverse(U)
    phase=[event for mode in range(16)for event in af.mb.occupation_phase_word(mode,float(-mu if mode==0 else mu))]
    return compiler.inverse(U)+phase+U


def error_contract(h,mass,steps):
    h,mass=mass_coin.parameters(h,mass)
    if type(steps)is not int or not 0<=steps<=3:raise ValueError('three-step development bound required')
    base=json.loads((HERE.parent/'vector_field/COMPOSITION.json').read_text())
    native=json.loads((mass_coin.REVIEW/'SOURCE_ERROR.json').read_text())
    uniform=F(base['source_preparation_error_exact'])-F(native['primitive_operator_errors_exact']['flip'])
    assert uniform>0
    mu=h*mass
    phase=F()
    for mode in range(16):
        angle=-mu if mode==0 else mu
        word=af.mb.occupation_phase_word(mode,float(angle))
        phase+=mass_coin.source_error.exact_template_bound(word,composition.phase_template(mode,angle))
    onsite=2*uniform+phase
    extra=128*onsite
    detector=particle_records.native_position_contract(1152,3203)
    reset=3203*(1-F(base['native_blank_reset_success_lower_exact']))
    total=(reset+F(base['source_preparation_error_exact'])
           +steps*(F(base['oneface_complete_V_error_exact'])+extra)
           +F(detector['apparatus_TV_upper_exact']))
    return {'steps':steps,'onsite_Gamma_mass_kick_operator_error_exact':str(onsite),
            'whole_Gamma_mass_kick_operator_error_exact':str(extra),
            'whole_step_operator_error_exact':str(extra+F(base['oneface_complete_V_error_exact'])),
            'complete_source_preparation_evolution_and_coarse_field_readout_TV_exact':str(total),
            'complete_TV_upper':math.nextafter(float(total),math.inf),
            'additional_native_pulses_per_step':128*len(fock_kick_word(h,mass)),
            'detector':detector,
            'reset_semantics':'Worst-conditional lower success bounds and a union bound. No independent reset assumption or exact product of marginal success factors is used.',
            'scope':'Exact retained source response and ideal coherent code services. Includes3203 source blank resets,7 fresh position-flag resets,1159 fresh native ready pointers, all source/export complements. No clock/probe instrument is included.'}


def run():
    paths=[Path(__file__),Path(af.__file__),Path(composition.__file__),Path(mass_coin.__file__),
           Path(particle_records.__file__),Path(response.__file__),
           HERE.parent/'vector_field/COMPOSITION.json',HERE/'MASS_COIN.json',HERE/'PARTICLE_RECORDS.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    rows=[];finals=[];complete_laws=[]
    for mass in (F(0),F(1)):
        h=F(1,8);state=af.initial();records=[];moments=[]
        for n in range(4):
            complete_laws.append({'mass_parameter':str(mass),'h':str(h),
                                  'alpha':.17,'beta':.41,'step':n,
                                  'joint_law':response.joint_law(state)})
            record=response.measure(state,1)
            assert abs(record['unconditional_norm_squared']-1)<2e-10
            record.update(step=n,time_exact=str(n*h),tau0_exact='1')
            records.append(record)
            if n<3:
                evolved=step(state,h,mass,.17,.41)
                overlap=af.phase_correlator(state,evolved)
                moments.append([float(overlap.real),float(overlap.imag)])
                state=evolved
        assert max(v[0]for v in moments)-min(v[0]for v in moments)<2e-10
        assert max(v[1]for v in moments)-min(v[1]for v in moments)<2e-10
        rows.append({'h':str(h),'mass_parameter':str(mass),'alpha':.17,'beta_face0':.41,
                     'other207_face_couplings':0.,'records':records,'conserved_V_moments':moments,
                     'native_complete_error':error_contract(h,mass,3)})
        finals.append(np.array(records[-1]['all96_position_probabilities']))
    mass_response=float(abs(finals[0]-finals[1]).sum()/2)
    assert mass_response>1e-5
    # The zero-mass interface must reproduce the actual same source.
    before=af.initial();expected=af.step(before,.17,.41);actual=step(before,F(1,8),F(0),.17,.41)
    assert set(expected)==set(actual)
    assert max(np.max(abs(expected[key]-actual[key]))for key in expected)<1e-14
    result={'schema':'oph.local12.massive_coupled_records.v1','candidate_id':MODEL_ID,
        'checks':{'same_autonomous_source_with_explicit_mass_coin':True,
                  'zero_mass_reproduces_frozen_field_source':True,
                  'all_complete_position_field_laws_and_Gauss':True,
                  'mass_intervention_changes_same_source_positions':True,
                  'same_total_Floquet_diagnostic_conserved':True,
                  'native_full_Fock_kick_and_complete_instrument_bound':True,
                  'source_bytes_unchanged':True},
        'selection':'common amendment proposal for coordinator review; no source substituted in predecessor receipts',
        'update':'V_m=U_E(alpha) U_B(beta) Gamma(S)_F Gamma(R exp(-i h*m R))',
        'native_implementation':'Prepend the emitted128 onsite Gamma(exp(-i h*m R)) words to the actual common matter-field macro. They preserve each local occupation and Gauss. The fused base one-particle implementation is separately supplied by mass_coin.',
        'experiments':rows,'mass_intervention_position_TV':mass_response,
        'complete_terminal_laws':response.write_joint_laws('MASSIVE_JOINT_LAWS.json.gz',complete_laws),
        'parameter_ancestry':'m,alpha,beta and tau0 are explicit common architecture inputs; no observed species masses, fine-structure constant or continuum equations are supplied as measured outputs.',
        'scope':'Complete finite one-face one-particle source with supplied mass coin, shared current/field and native readout/error composition. The free massive dispersion theorem does not establish an interacting pole or renormalized physical mass. Homogeneous interacting continuum, spin and observed species remain separate same-source obligations.',
        'source_sha256':pins,'official_rows_changed':[]}
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    (HERE/'MASSIVE_COUPLED.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'MASSIVE_KICK_PROGRAM.json').write_text(json.dumps({'candidate_id':MODEL_ID,
        'onsite16_mode_word':fock_kick_word(F(1,8),F(1)),
        'all128_groups':'remap localmodej to16*position_address+j',
        'other_source_words':'exact vector_field/FIELD_PROGRAM and composition source words',
        'source_sha256':pins},indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'mass_intervention_position_TV':mass_response,
                      'complete_native_TV_bounds':[r['native_complete_error']['complete_TV_upper']for r in rows]}))


if __name__=='__main__':run()
