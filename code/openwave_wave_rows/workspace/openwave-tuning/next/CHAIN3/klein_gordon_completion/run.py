"""Complete native spectroscopy development and massive refinement checks."""
from fractions import Fraction as F
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import math
import sys
import numpy as np
from scipy.linalg import expm
import native_spectroscopy as ns

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
a,b=ns.a,ns.b


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data): (HERE/name).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
def peak(x): return float(np.max(abs(x)))
def complex_array(x): return {'real':np.asarray(x).real.tolist(),'imag':np.asarray(x).imag.tolist()}


def pins():
    result={}
    packages={
       'openwave-tuning/next/CHAIN3/larger_walk_connection':'ed7803941a479a257067202d2f758fbb4468a32e318ba13eadaf25f52e4af35c',
       'openwave-tuning/next/CHAIN2/particle_field_limit_audit':'fba3a996050b09398aebc44ebc0d3fa8000911581ce06789a3bcbf17d3f5ac08',
       'openwave-tuning/next/CHAIN2/local12_limit_audit':'997916c3b216f3006b9c88b8ef39419a213df5075e8159f47298d87ad295e2ad',
       'openwave-tuning/common_dynamics/particles':'678983fa40a2a97b8eb399d6b5ad4ab6ace77a5f441baf448cb0c00bc8a6f9a9',
       'openwave-tuning/common_dynamics/continuum':'3067f7bdea43ca5137d6b115f7e1bd92b4d753ba99684ed198ba2da170240adf',
    }
    for relative,h in packages.items():
        p=ROOT/relative/'MANIFEST.json'
        assert sha(p)==h,relative
        result[str(p.relative_to(ROOT))]=h
    result.update(ns.walk.verify_sources())
    for mod in (ns,ns.walk,ns.native_word,a,b,ns.walk.mass_source,b.continuum):
        p=Path(mod.__file__).resolve();result[str(p.relative_to(ROOT))]=sha(p)
    return result


def full_native(word,bits):
    result=b.complete_run(word,bits,np.eye(1<<bits),96)
    b.require_complete(result)
    gram=peak(result['success'].conj().T@result['success']+sum(result['failure_grams'])-np.eye(1<<bits))
    assert gram<2e-10
    return result,gram


def primitive_checks():
    rows=[]
    for kind in ('X','Y','Z','I','H','J'):
        count=2 if kind=='J' else 1
        angle=a.angle(0,0,F(1,7))
        word=ns.controlled_event(kind,tuple(range(count)),angle,count)
        result,gram=full_native(word,count+1)
        local=a.execute_full_response(a.certified_pulse(kind,tuple(range(count)),angle,96))[0][:(1<<count)]
        expected=np.eye(1<<(count+1),dtype=complex);expected[1<<count:,1<<count:]=local
        error=peak(result['success']-expected);assert error<2e-12
        rows.append({'kind':kind,'complete_columns':1<<(count+1),'pulses':len(word),
                     'operator_residual':error,'complete_gram_residual':gram,
                     'all_export_complements_retained':True,'exact_operator_bound':str(a.template_bound(word,96))})
    for mass in (F(0),F(1,2),F(1)):
        word=list(ns.controlled_word(ns.native_word.coin_template(F(1,4),mass),4))
        result,gram=full_native(word,5)
        expected=np.eye(32,dtype=complex);expected[16:,16:]=ns.walk.mass_source.coin_matrix(F(1,4),mass)
        error=peak(result['success']-expected);assert error<2e-11
        rows.append({'kind':'full_controlled_coin','mass_exact':str(mass),'complete_columns':32,'pulses':len(word),
                     'operator_residual':error,'complete_gram_residual':gram,
                     'all_export_complements_retained':True,'exact_operator_bound':str(a.template_bound(word,96))})
    return rows


def bloch_source(k,mu):
    """Direct actual seam generator; no spectral/dispersion input."""
    S=np.eye(144,dtype=complex)
    for direction,matching in enumerate(a.common.INTERNAL_MATCHINGS):
        e=2*direction+1
        for x,y in matching:
            i,j=12*x+e,12*y+e;S[i,i]=S[j,j]=0;S[i,j]=S[j,i]=1
    for axis,row in enumerate(a.spatial.links()):
        for x,y in row:
            i,j=12*x+11,12*y+11;S[i,i]=S[j,j]=0
            S[i,j],S[j,i]=np.exp(1j*k[axis]),np.exp(-1j*k[axis])
    R=np.kron(np.eye(12),np.ones((12,12))/6-np.eye(12))
    return S@(math.cos(mu)*R-1j*math.sin(mu)*np.eye(144)),S,R


def readout(rho,quadrature,law,clock_steps=0,size=2):
    """Actual M6 analysis plus native C/E pointer, all fine branches."""
    word=[]
    if clock_steps:
        word += [('Z',(0,),a.angle(0,0,-F(1,size*16)))]*clock_steps
    if quadrature=='Y': word += [('Z',(0,),a.angle(F(1,4)))]
    word += b.H(0)
    analysis=b.complete_run(word,1,np.eye(2),96)
    alive=analysis['success']@rho@analysis['success'].conj().T
    rows=[]
    for i,g in enumerate(analysis['failure_grams']):
        rows.append({'outcome':'analysis_export_failure','index':i,'probability':float(np.trace(g@rho).real)})
    families,pointer=b.pointer_blocks(96)
    for i,g in enumerate(pointer['failure_grams']):
        rows.append({'outcome':'pointer_export_failure','index':i,'probability':float(np.trace(g@alive).real)})
    for meta,ks in families[law]:
        for environment,K in enumerate(ks):
            rows.append({'outcome':'native_pointer','record':meta,'environment':environment,
                         'probability':float(np.trace(K@alive@K.conj().T).real)})
    assert abs(sum(row['probability'] for row in rows)-np.trace(rho).real)<2e-11
    assert min(row['probability'] for row in rows)>-1e-14
    return rows


def contrast(rows):
    return sum((1-2*r['record']['bit'])*r['probability'] for r in rows if r['outcome']=='native_pointer')


def finite_experiments():
    rows=[];fine=[]
    for size in (2,4):
        _,d=b.member(size)
        for mode in ((1,0,0),(1,1,0),(1,1,1)):
            word=ns.preparation_template(size,mode)
            initial=np.zeros(1<<d,complex);initial[0]=1
            prep=b.complete_run(word,d,initial,96);b.require_complete(prep)
            psi=prep['success'][:,0]
            expected=np.zeros((size**3,16,16),complex)
            for cell in range(size**3):
                x,y,z=cell//(size*size),(cell//size)%size,cell%size
                expected[cell,:12,:12]=np.exp(2j*np.pi*np.dot(mode,(x,y,z))/size)/(12*size**1.5)
            prep_error=float(np.linalg.norm(psi-expected.reshape(-1)))
            assert prep_error<2e-11
            prep_norm=float(np.vdot(psi,psi).real)
            preparation_failures=[float(np.trace(g).real) for g in prep['failure_grams']]
            assert abs(prep_norm+sum(preparation_failures)-1)<2e-11
            for mass in (F(0),F(1,2),F(1)):
                source=ns.walk.Walk(size,h=F(1,size),mass=mass)
                state=psi.reshape(source.geometry.shape).copy()
                W,_,_=bloch_source(2*np.pi*np.asarray(mode)/size,float(mass)/size)
                e=np.ones(144)/12
                fibre=e.copy()
                for step in range(1,5):
                    state=source.step(state);fibre=W@fibre
                    c=np.vdot(psi,state.reshape(-1))
                    c_ideal=np.vdot(e,fibre)
                    assert abs(c-c_ideal)<2e-11
                    rho=np.array([[prep_norm,c.conjugate()],[c,prep_norm]])/2
                    records={}
                    for quadrature in ('X','Y'):
                        for law in ('C','E'):
                            rr=readout(rho,quadrature,law)
                            clock=readout(np.ones((2,2))/2,'X',law,step,size)
                            observed=contrast(rr)
                            target=c.real if quadrature=='X' else c.imag
                            assert abs(observed-target)<2e-11
                            assert abs(contrast(clock)-math.cos(step/(8*size)))<2e-11
                            # These are the factorized full native terminal
                            # laws at the certified compressed source centre.
                            joint=[{'matter':i,'clock':j,'probability':r['probability']*s['probability']}
                                   for i,r in enumerate(rr) for j,s in enumerate(clock)]
                            total=sum(x['probability'] for x in joint)+sum(preparation_failures)
                            assert abs(total-1)<3e-11
                            records[quadrature+law]={'matter':rr,'clock':clock,'joint':joint,
                                'preparation_export_failures':preparation_failures,
                                'normalization_residual':abs(total-1)}
                    fine.append({'size':size,'mode':list(mode),'mass_exact':str(mass),'steps':step,'laws':records})
                rows.append({'size':size,'mode':list(mode),'mass_exact':str(mass),
                    'native_preparation_error':prep_error,'native_preparation_pulses':len(word),
                    'steps_executed':4,'full_padded_dimension':psi.size,
                    'final_return_amplitude':[float(c.real),float(c.imag)],
                    'source_evaluation':'actual finite source array; controlled native word compression with complete operator bound',
                    'no_continuum_claim_for_small_member':True})
        print('Finite native preparation/source/pointer member passed',size,flush=True)
    save('COMPLETE_RECORDS.json',fine)
    return rows


def exact_continuum_algebra():
    data=b.continuum.rational_bloch_data()
    denominator=data['exact_projector_denominator']
    Pn=np.rint(data['P']*denominator).astype(np.int64)
    Rn=np.rint(data['R']*6).astype(np.int64)
    Hn=[np.rint(h*1152).astype(np.int64) for h in data['H']]
    one=np.ones(144,dtype=np.int64)
    assert np.array_equal(Pn@Pn,denominator*Pn)
    assert np.array_equal(Rn@Pn,Pn@Rn)
    # RP = 2 ee* - P, with e=one/12.
    assert np.array_equal(Rn@Pn,denominator//12*np.ones((144,144),dtype=np.int64)-6*Pn)
    for h in Hn: assert not np.any(Rn@h+h@Rn)
    g=[h@one for h in Hn]
    for x in range(3):
        for y in range(3):
            assert int(g[x]@g[y])==(3*1152**2 if x==y else 0)
    return {'rank_P0':int(np.trace(Pn)//denominator),'exact_denominator':denominator,
            'mass_reflection_preserves_P0':True,'mass_anticommutes_with_three_gradients':True,
            'gradient_Gram_exact':'delta_ab/48','propagating_branches':2,'stationary_slow_branches':24,
            'all_source_outputs_retained':True}


def continuum_experiments():
    data=b.continuum.rational_bloch_data();e=np.ones(144)/12
    R=data['R'];P=data['P'];rows=[];mutations={}
    A=np.kron(np.eye(12),np.ones((12,1))/math.sqrt(12))
    for mode in ((1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,1,1),(1,2,0)):
        p=2*np.pi*np.array(mode);K=F(7 if max(mode)==1 else 13)
        for mass in (F(0),F(1,2),F(1)):
            for size in (32,128,512,2048):
                h=F(1,size);n=size
                W,S,_=bloch_source(float(h)*p,float(h*mass))
                target=b.continuum.bloch_walk(float(h)*p)@expm(-1j*float(h*mass)*R)
                source_residual=peak(W-target);assert source_residual<3e-14
                T=A.conj().T@S@A
                values,vectors=np.linalg.eigh(T)
                t=float(values[-1]);f=A@vectors[:,-1]
                spectral_identity=float(np.linalg.norm(W@W@f-2*t*math.cos(float(h*mass))*W@f+f))
                assert spectral_identity<2e-13
                omega=math.acos(max(-1.,min(1.,t*math.cos(float(h*mass)))))
                kg=float(mass)**2+float(p@p)/48
                deviation=abs((omega/float(h))**2-kg)
                law=ns.walk.mass_source.continuum_bound(h,mass,K,n)
                assert deviation<=float(F(law['squared_frequency_error_upper_exact']))+1e-8
                actual=np.linalg.matrix_power(W,n)@e
                gen=1j*sum(p[i]*data['H'][i] for i in range(3))-1j*float(mass)*(R@P)
                limit=expm(gen)@e
                state_error=float(np.linalg.norm(actual-limit))
                assert state_error<=float(F(law['uncapped_vector_bound_exact']))+2e-10
                c=np.vdot(e,actual)
                energy=math.sqrt(kg)
                expected=math.cos(energy)-1j*float(mass)*math.sin(energy)/energy
                assert abs(np.vdot(e,limit)-expected)<1e-12
                rows.append({'mode':list(mode),'mass_exact':str(mass),'size':size,
                    'omega_per_source_macro':omega,'energy_squared':(omega/float(h))**2,
                    'KG_target':kg,'squared_dispersion_deviation':deviation,
                    'squared_dispersion_bound':float(F(law['squared_frequency_error_upper_exact'])),
                    'source_matrix_residual':source_residual,'spectral_polynomial_residual':spectral_identity,
                    'full_state_continuum_error':state_error,'full_state_bound':float(F(law['uncapped_vector_bound_exact'])),
                    'actual_return':[float(c.real),float(c.imag)],'KG_return':[expected.real,expected.imag],
                    'all_optical_components_kept':True})
                if mode==(1,0,0) and mass==1 and size==2048:
                    correct=expected
                    mutations['double_time_return_difference']=abs(c-(math.cos(2*energy)-1j*float(mass)*math.sin(2*energy)/energy))
                    mutations['massless_substitution_return_difference']=abs(c-math.cos(math.sqrt(float(p@p)/48)))
                    mutations['source_mass_sign_difference']=abs(c-correct.conjugate())
                    assert min(mutations.values())>.1
            print('Controlled KG direction/mass passed',mode,str(mass),flush=True)
    # Flat cycle and sharp-record processes are retained counterexamples to
    # overbroad assertions that every source component obeys the KG branch.
    P=np.asarray(data['P']);scalar=np.outer(e,e)
    transverse=P-scalar
    assert int(round(np.trace(transverse)))==25
    return rows,mutations


def formal_family():
    rows=[]
    for size in (32,128,512,2048,8192):
        row=ns.certificate(size,size,F(1),(1,0,0),bandwidth=F(7))
        h=F(1,size)
        assert F(row['delta_preparation_exact'])<=h*h
        assert F(row['delta_controlled_step_exact'])<=h*h
        assert F(row['delta_clock_and_readout_exact'])<=h*h
        rows.append(row)
    assert all(F(x['native_vs_continuum_uncapped_exact'])>F(y['native_vs_continuum_uncapped_exact']) for x,y in zip(rows,rows[1:]))
    assert F(rows[-1]['native_vs_continuum_uncapped_exact'])<F(1,20)
    return rows


def main():
    before=pins()
    primitives=primitive_checks();print('Controlled native primitives and full coins passed',flush=True)
    algebra=exact_continuum_algebra();print('Exact full-sector massive generator passed',flush=True)
    finite=finite_experiments()
    continuum,mutations=continuum_experiments()
    family=formal_family();print('Complete native refinement family passed',flush=True)
    # Mutated native controls and missing complement accounting must reject.
    event=a.certified_pulse('Z',(0,),a.angle(F(1,4)),96);event['controls_exact'][0]='1'
    try:a.verify_pulse(event)
    except ArithmeticError:mutations['native_control_mutation_rejected']=True
    else:raise AssertionError('changed controls accepted')
    result={'schema':'oph.row30.native_kg_completion.v1','candidate_id':'native-local12-mass-coin-rational-spectroscopy-v1',
      'criterion':'emergent ω² = k² + m² for matter waves','native_primitive_checks':primitives,
      'exact_continuum_algebra':algebra,'finite_native_instruments':finite,'continuum_spectroscopy':continuum,
      'formal_complete_native_family':family,'controls':mutations,'full_criterion_candidate':True,
      'full_criterion_met':False,'two_independent_full_reviews':False,'official_rows_changed':[],
      'declared_architecture_not_uniqueness_proof':True,'mass_parameter_is_architecture_input':True,
      'native_full_source_streams_executed':False,'source_compression_error_certified':True,
      'all_sources_unchanged':pins()==before,'source_sha256':before,
      'scope':'Declared coherent native processor architecture; scalar prepared matter-wave sector, full source outputs retained, all failure outcomes. No interacting field-limit or empirical mass claim.'}
    assert result['all_sources_unchanged']
    save('RESULT.json',result)
    imported={}
    for mod in tuple(sys.modules.values()):
        path=getattr(mod,'__file__',None)
        if path:
            p=Path(path).resolve()
            if p.suffix=='.pyc' and p.with_suffix('.py').exists():p=p.with_suffix('.py')
            if p.is_file() and p.is_relative_to(ROOT) and not p.is_relative_to(HERE):
                imported[str(p.relative_to(ROOT))]=sha(p)
    save('IMPORTS.json',{'files_sha256':imported})
    print(json.dumps({'full_criterion_candidate':True,'finite_instruments':len(finite),
                      'continuum_cases':len(continuum),'native_family_TV_upper':float(F(family[-1]['native_vs_continuum_uncapped_exact']))}),flush=True)


if __name__=='__main__':main()
