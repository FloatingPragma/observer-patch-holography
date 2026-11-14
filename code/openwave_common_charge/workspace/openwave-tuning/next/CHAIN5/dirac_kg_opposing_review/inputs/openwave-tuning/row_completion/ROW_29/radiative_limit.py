"""Native terminal radiative instruments on a controlled finite-box family."""
from fractions import Fraction as F
from pathlib import Path
from copy import deepcopy
import hashlib,json,math,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'openwave-tuning/common_dynamics/photon_size3'))
import providers,walsh,native_member
import maxwell
import numpy as np
nf=walsh.nf;iv=nf.iv;compiler=nf.compiler;cert=nf.cert
clock=native_member.phase_clock
from model.native_clock import configuration
from model.observer.source_discrimination import circuit
from model.native_interface import pointer


def phase_word(time,scale=None):
    p=2*iv.pi();p=p if scale is None else p*scale
    return walsh.compile_intervals([-time*p,-2*time*p,-2*time*p,iv.Interval.point(0)])


def apparatus(bit,time,position,setting,clock_setting,diagonal):
    position_word=[compiler.pulse('Z',-math.pi*float(position),[bit])]
    return compiler.remap(nf.H(0),[bit])+diagonal+position_word+native_member.analysis(setting,bit)+compiler.remap(clock.local_words(1,clock_setting,time),[2])


def effects_and_law(v,branches):
    effects=[np.zeros((2,2),complex),np.zeros((2,2),complex)];fine=[]
    for metadata,kraus in branches:
        e=sum((K.conj().T@K for K in kraus),np.zeros((2,2),complex));effects[metadata['bit']]+=e
        fine.append({'record':metadata,'effect_real':e.real.tolist(),'effect_imag':e.imag.tolist()})
    transition=np.array([[effects[b][x,x].real for x in (0,1)] for b in (0,1)])
    assert max(abs(effects[b][0,1]) for b in (0,1))<1e-13
    probabilities=abs(v)**2
    for bit in range(3):
        lower=np.array([x for x in range(8) if not x&(1<<bit)]);indices=np.array([lower,lower+(1<<bit)])
        probabilities[indices]=transition@probabilities[indices]
    assert abs(sum(probabilities)-np.vdot(v,v).real)<1e-11
    return probabilities,fine


def run():
    before=providers.verify();execute=native_member.independent_executor();instruments,_=clock.actual_native_pointer()
    p=2*iv.pi();plo=F(p.lo,iv.SCALE);phi=F(p.hi,iv.SCALE)
    p0=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(0))))
    pr=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(pointer.READY_CODE))))
    reset=3*(1-p0)+3*(1-pr)
    H=[('I',(0,),cert._pi(-1,2)),('H',(0,),cert._pi(1,2))]
    herror=nf.root_certificate.exact_template_bound(nf.H(0),H)
    rows=[];programs={};contrasts={}
    for time in (F(1,8),F(1,4)):
      diagonal,phase_error=phase_word(time)
      for bit in (0,1):
       for position in (F(0),F(1,4)):
        for setting in ('X','Y'):
         for clock_setting in ('X','Y'):
          name=f'T{time}_pol{bit}_x{position}_{setting}{clock_setting}'
          word=apparatus(bit,time,position,setting,clock_setting,diagonal);programs[name]=word
          v=execute(word,3,np.eye(8,dtype=complex)[:,0])[:,0]
          angle=2*math.pi*float(position-time)
          target=math.cos(angle) if setting=='X' else math.sin(angle)
          cp=float(time)/8;clock_target=math.cos(cp) if clock_setting=='X' else math.sin(cp)
          pw=[compiler.pulse('Z',-math.pi*float(position),[0])]
          pe=nf.root_certificate.exact_template_bound(pw,[('Z',(0,),cert._pi(-position))])
          ye=nf.root_certificate.exact_template_bound([compiler.pulse('Z',math.pi/4,[0])],[('Z',(0,),cert._pi(1,4))]) if setting=='Y' else F()
          native=phase_error+2*herror+pe+ye+clock.exact_word_bound(1,clock_setting,time)+2*clock.pointer_bound()
          laws={}
          for law,branches in instruments.items():
            probabilities,fine=effects_and_law(v,branches)
            q=float(sum(z*(-1 if j&(1<<bit) else 1) for j,z in enumerate(probabilities)))
            t=float(sum(z*(-1 if j&4 else 1) for j,z in enumerate(probabilities)))
            assert abs(q-target)<2*float(native)+1e-11 and abs(t-clock_target)<2*float(native)+1e-11
            laws[law]={'all8_probabilities':probabilities.tolist(),'field_contrast':q,'clock_contrast':t,'fine_one_bit_effects':fine}
            contrasts[time,bit,position,setting,clock_setting,law]=(q,t)
          rows.append({'name':name,'time':str(time),'polarization':'y' if bit==0 else 'z','position_x':str(position),'field_setting':setting,'clock_setting':clock_setting,
            'native_pulses':len(word),'word_pointer_error_exact':str(native),'reset_failure_exact':str(reset),'native_complete_TV_exact':str(native+reset),'laws':laws})
    refinements=[]
    for exponent in (34,36,38):
      h=F(1,2**exponent);r=phi*h;eta,D,zeta=maxwell.bounds(r)
      dw=F(5,4)*D/h
      # Low-source vacuum scalar and cap-one excitation both contribute.
      Herror=2*dw;field_relative=64*(zeta+dw/plo)
      assert dw<plo/4 and field_relative<F(1,100)
      bound=max(F(row['native_complete_TV_exact'])+F(row['time'])*Herror for row in rows)
      refinements.append({'L':2**exponent,'h':str(h),'allowed_mode_integer':[1,0,0],'physical_p_interval':p.json(),'periodic_box_side':1,
          'exact_macro_counts':{str(T):int(T/h) for T in (F(1,8),F(1,4))},
          'source_frequency_error_exact':str(dw),'source_cap_H_error_including_visible_vacuum_exact':str(Herror),
          'polarization_subspace_error_exact':str(zeta),'local_field_relative_error_exact':str(field_relative),
          'complete_native_source_TV_upper_exact':str(bound),'complete_native_source_TV_upper':float(bound),
          'normalized_field_mean_error_upper':float(2*bound+field_relative)})
    assert refinements[-1]['complete_native_source_TV_upper']<refinements[0]['complete_native_source_TV_upper']/15
    recovered=[]
    for time in (F(1,8),F(1,4)):
      for bit in (0,1):
        X=contrasts[time,bit,F(0),'X','X','C'][0];Y=contrasts[time,bit,F(0),'Y','X','C'][0]
        cx=contrasts[time,bit,F(0),'X','X','C'][1];cy=contrasts[time,bit,F(0),'X','Y','C'][1]
        phase=math.atan2(-Y,X);Tread=8*math.atan2(cy,cx);speed=phase/(2*math.pi*Tread)
        assert abs(speed-1)<1e-10
        recovered.append({'T':str(time),'polarization':bit,'field_phase':phase,'clock_time':Tread,'speed_in_common_units':speed})
    # Actual wrong-source word: double magnetic energy gives sqrt2 frequency.
    wrong,unused=phase_word(F(1,4),iv.sqrt_fraction(2))
    wrong_word=apparatus(0,F(1,4),F(0),'X','X',wrong)
    badv=execute(wrong_word,3,np.eye(8,dtype=complex)[:,0])[:,0]
    badp,_=effects_and_law(badv,instruments['C']);badcontrast=float(sum(v*(-1 if j&1 else 1) for j,v in enumerate(badp)))
    assert abs(badcontrast)>0.5
    # Clock phase coalescing is an exact terminal map, not a reuse of old costs.
    cw=clock.local_words(8,'Y',F(1,32));one=clock.local_words(1,'Y',F(1,4))
    assert np.linalg.norm(execute(cw,1,np.array([1.,0.]))-execute(one,1,np.array([1.,0.])))<1e-12
    programs['wrong_double_b']=wrong_word;programs['clock_split_control']=cw
    assert providers.verify()==before
    (HERE/'RADIATIVE_PROGRAM.json').write_text(json.dumps(programs,separators=(',',':'))+'\n')
    out={'schema':'oph.row29.native_radiative_refinement.v1','checks':{'actual_two_field_quadratures_and_clock':True,'all_output_and_fine_failure_contract':True,'finite_periodic_allowed_momenta':True,'controlled_source_field_and_frequency_limit':True,'both_transverse_polarizations':True,'common_native_clock_speed_one':True,'double_b_mutation_rejected':True,'terminal_clock_composition':True,'source_bytes_unchanged':True},
      'experiments':rows,'refinements':refinements,'native_recovered_speeds':recovered,'mutant_double_b_X_contrast':badcontrast,
      'source_domain':'Exact two-low-mode radiative observable algebra of the noncompact v4 finite periodic source, total observed occupation at most1, arbitrary normalized unobserved/source reference marginal. Source and native preparation are the specified spectral-code isometry.',
      'observables':'A_e=1/3. At volume1, A amplitude=sqrt(A_e/(2p)), E and B amplitude=sqrt(A_e*p/2). A=e_pol*X, E=-e_pol*Y, B=-(e_x cross e_pol)*Y with these amplitudes. These are first moments of the selected radiative Fourier component in its vacuum/one-photon domain; they equal full-field means when spectator first moments vanish.',
      'terminal_contract':'Only terminal field and clock observations. Free source macro powers are coalesced by exact group law; processor pulses and their retained export records are those emitted here, not transferred dense-lattice or per-macro costs.',
      'complete_fine_law':'All8 basis outcomes and tensor products of3 one-bit fine effects on the unnormalized native success vector, plus first export failures and3 blank/3 pointer-ready reset failures. No success normalization.',
      'phase_reference':'Visible two-mode E0=Tr(Omega)/2 retained and bounded. Unobserved coexact and harmonic factors are traced; their scalar phase cancels by autonomy. No coherent padded reference or controlled full-field evolution transfer.',
      'native_precision_boundary':'The exact continuum source limit is h->0. Fixed native word/pointer precision adds the stated nonzero instrument floor; no claim that it vanishes at fixed control precision.',
      'program_sha256':providers.sha(HERE/'RADIATIVE_PROGRAM.json'),'source_sha256':before}
    (HERE/'RADIATIVE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':out['checks'],'TV':[r['complete_native_source_TV_upper'] for r in refinements],'speed':[r['speed_in_common_units'] for r in recovered]}),flush=True)
if __name__=='__main__':run()
