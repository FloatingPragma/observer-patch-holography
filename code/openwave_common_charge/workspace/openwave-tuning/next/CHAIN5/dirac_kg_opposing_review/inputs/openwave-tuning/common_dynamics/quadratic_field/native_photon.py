"""Native cap-one free-photon phase, common clock and complete pointers.

This is the declared spectral-code adapter, not a local electric mesh or an
interacting-field implementation. All 190 modes and 65 padding states are
part of the eight-bit source word; its scalar vacuum phase is retained.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

import source as field
import fock

HERE=Path(__file__).resolve().parent
ROOT=field.ROOT
import native_field as nf
compiler=nf.compiler;cert=nf.cert;iv=nf.iv
sys.path.insert(0,str(HERE.parent/'clock'))
import phase_clock
from model.reference_recovery import recovery
from model.native_clock import configuration
from model.observer.source_discrimination import circuit
from model.native_interface import pointer


def diagonal_bound(word,desired):
    """Integer-outward Walsh-phase certificate for the complete 256 code."""
    assert len(desired)==256
    offset=0;error=F()
    for mask in range(256):
        label=''.join('Z' if mask>>j&1 else 'I' for j in range(8))
        offset,fixed,variables=cert.consume(word,offset,cert.pauli(label))
        assert len(variables)==1
        wanted=sum((value*(-1)**((i&mask).bit_count()) for i,value in enumerate(desired)),iv.Interval.point(0))/256
        error+=fixed+(iv.Interval.point(variables[0])-wanted).absolute_upper()
    assert offset==len(word)
    return error+cert.controls_bound(word,{})


def photon_phase(time=F(1)):
    if type(time) is not F or not 0<time<=8:raise ValueError('exact time in (0,8] required')
    spec=fock.spectral_certificate()
    frequencies=[iv.sqrt_fraction(F(3,2)*x) for x in spec['eigenvalues_exact']]
    vacuum=sum(frequencies,iv.Interval.point(0))/2
    energies=[vacuum]+[vacuum+w for w in frequencies]+[iv.Interval.point(0)]*65
    desired=[-time*x for x in energies]
    phases=[float(F(x.lo+x.hi,2*iv.SCALE)) for x in desired]
    word=compiler.compile_diagonal(phases)
    return word,diagonal_bound(word,desired),np.asarray(phases),vacuum


def word(setting='X',time=F(1)):
    if setting not in ('X','Y'):raise ValueError('configured clock X or Y required')
    diagonal,error,phases,vacuum=photon_phase(time)
    H=nf.H(0)
    template=[('I',(0,),cert._pi(-1,2)),('H',(0,),cert._pi(1,2))]
    herror=nf.root_certificate.exact_template_bound(H,template)
    clock=compiler.remap(phase_clock.local_words(1,setting,time),[8])
    # Clock preparation and analysis commute with the complete source word.
    program=H+diagonal+H+clock
    native=error+2*herror+phase_clock.exact_word_bound(1,setting,time)+8*phase_clock.pointer_bound()
    return program,native,phases,vacuum


def complete_bit_probabilities(vector,effects,bits):
    """Contract all tensor pointer effects without making 512 dense effects."""
    rho=np.outer(vector.conj(),vector).reshape([2]*(2*bits))
    effect=np.asarray(effects)
    for done in range(bits):
        rho=np.tensordot(rho,effect,axes=([0,bits-done],[1,2]))
    assert np.max(abs(rho.imag))<1e-12
    return rho.real.reshape(-1)


def run():
    paths=[Path(__file__),Path(field.__file__),Path(fock.__file__),Path(nf.__file__),
      Path(compiler.__file__),Path(cert.__file__),Path(iv.__file__),Path(recovery.__file__),
      Path(phase_clock.__file__),Path(configuration.__file__),
      HERE/'FOCK.json',HERE/'SPECTRAL_FRAME.json',
      HERE.parent/'vector_field/RESET_DOMAIN.json',HERE.parent/'vector_field/MANIFEST.json',
      ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    selected=fock.spectral_certificate()
    assert selected['frame_integers']==json.loads((HERE/'SPECTRAL_FRAME.json').read_text())['frame_integers']
    assert [str(x) for x in selected['eigenvalues_exact']]==json.loads((HERE/'FOCK.json').read_text())['spectral_certificate']['eigenvalues_exact']
    instruments,_=phase_clock.actual_native_pointer()
    p0=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(0))))
    pr=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(pointer.READY_CODE))))
    reset=9*(1-p0)+9*(1-pr)
    rows=[];programs={};worst=0.;time=F(1)
    for setting in ('X','Y'):
        program,native,phases,vacuum=word(setting,time)
        actual=recovery.execute_complete(program,9,np.eye(512,dtype=complex)[:,[0]])
        v=actual.success[:,0]
        initial=np.zeros(256,complex);initial[0]=initial[1]=1/math.sqrt(2)
        evolved=initial*np.exp(1j*phases)
        photon=compiler.execute(nf.H(0),8,evolved)
        clock=compiler.execute(phase_clock.local_words(1,setting,time),1,np.array([1,0],complex))
        expected=np.kron(clock,photon)
        discrepancy=float(np.linalg.norm(v-expected));worst=max(worst,discrepancy)
        assert discrepancy<float(native)+2e-11
        laws={}
        for law,branches in instruments.items():
            effects=[np.zeros((2,2),complex),np.zeros((2,2),complex)];fine=[]
            for metadata,kraus in branches:
                effect=sum((K.conj().T@K for K in kraus),np.zeros((2,2),complex))
                effects[metadata['bit']]+=effect
                fine.append({'record':metadata,'effect_real':effect.real.tolist(),'effect_imag':effect.imag.tolist()})
            probabilities=complete_bit_probabilities(v,effects,9)
            assert abs(float(probabilities.sum())-float(np.vdot(v,v).real))<1e-10
            assert min(probabilities)>-1e-12
            photon_plus=float(sum(p for i,p in enumerate(probabilities) if not i&1))
            clock_plus=float(probabilities[:256].sum())
            target_photon=(1+math.cos(float(phases[1]-phases[0])))/2
            target_clock=(1+(math.cos(float(time)/8) if setting=='X' else math.sin(float(time)/8)))/2
            assert abs(photon_plus-target_photon)<1e-9 and abs(clock_plus-target_clock)<1e-9
            laws[law]={'all512_joint_probabilities':probabilities.tolist(),
              'photon_plus':photon_plus,'clock_plus':clock_plus,'fine_one_bit_effects':fine,
              'complete_fine_rule':'v* (tensor in descending bit order of nine retained fine effects) v; first reset timeout and every export complement remain additional outcomes',
              'success_normalization':False}
        geometry=time*fock.spectral_certificate()['cap1_total_H_error_including_vacuum']
        complete=reset+native+geometry
        rows.append({'setting':setting,'time_exact':str(time),'native_pulses':len(program),
          'numeric_state_residual':discrepancy,'native_export_failure_probability':float(actual.failure_gram[0,0].real),
          'native_word_and_pointer_bound_exact':str(native),'rational_geometry_H_bound_exact':str(geometry),
          'complete_record_TV_bound_exact':str(complete),'complete_record_TV_bound_upper':math.nextafter(float(complete),math.inf),
          'vacuum_energy_interval':vacuum.json(),'laws':laws})
        programs[setting]=program
    # Same-source frequency mutation is separated by an actual retained bit.
    separation=(1-math.cos(float(phases[1]-phases[0])))/2
    assert separation>.05 and separation>20*max(r['complete_record_TV_bound_upper'] for r in rows)
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'NATIVE_PHOTON_PROGRAM.json').write_text(json.dumps({'words':programs,
      'source_bits':list(range(8)),'clock_bit':8,'fresh_pointer_auxiliary_bits':[9,10],
      'blank_resets':9,'ready_pointer_resets':9,'physical_code':{'0':'vacuum','1..190':'one photon in declared rational spectral frame','191..255':'retained padding with H=0'},
      'zero_point_phase':'retained in all191 physical energies; not silently subtracted',
      'complete_failure_rule':'Retain first reset timeout, every source export complement, all pointer fine records and all padding outputs; do not condition on successful source export.'},indent=2)+'\n')
    result={'schema':'oph.actual_carrier.native_quadratic_photon.v1','candidate_id':field.MODEL_ID,
      'checks':{'all190modes_and65padding_in_native_phase_word':True,'vacuum_scalar_phase_retained':True,
        'integer_outward_phase_and_geometry_certificate':True,'native9bit_source_and_shared_clock_executed':True,
        'both_native_C_and_E_pointer_laws':True,'all512_joint_outcomes_retained':True,
        'all_native_export_and_reset_failures_retained':True,'frequency_deletion_control_separated':True,
        'no_interacting_cap1_accuracy_inferred':True,'source_bytes_unchanged':True},
      'experiments':rows,'maximum_numeric_state_residual':worst,
      'initialization_failure_union_bound_exact':str(reset),'frequency_deletion_record_separation':separation,
      'scope':'A native spectral-code free-photon adapter for total photon cap1 on the actual190mode field. The three harmonic free packets remain in the continuous source and have not been discretized here. It is not a local-edge field encoder, an interacting Weyl-current word, a massive matter state or an official row validation.',
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'NATIVE_PHOTON.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'native_pulses':[r['native_pulses'] for r in rows],
      'complete_TV_bounds':[r['complete_record_TV_bound_upper'] for r in rows],'numeric_residual':worst,
      'frequency_deletion_separation':separation}),flush=True)


if __name__=='__main__':run()
