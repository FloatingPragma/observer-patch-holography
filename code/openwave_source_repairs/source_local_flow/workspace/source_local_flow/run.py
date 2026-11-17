#!/usr/bin/env python3
"""Prospective source-local words, complete error and unfiltered phase reader."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,time,sys,traceback
import numpy as np
from scipy.linalg import expm
import local_flow as lf
import native_catalog as nc
import rounding
from finalize import finish
from native_event import native
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,obj):p.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
def log(x):print(json.dumps(x,sort_keys=True),flush=True)

def source_census(edges,bases):
    counts=np.zeros((240,240),dtype=np.int64);calls=0
    for p in range(240):
        q=[0]*240;q[p]=1
        for a,b in edges:
            for c in (0,1):
                out,record=native(q,(a,b),c);j=record[1]
                equivalent,k=lf.endpoint(q,a,b,c);assert list(out)==equivalent and j==k
                counts[out.index(1),p]+=1;calls+=1
    L=1260*np.eye(240,dtype=np.int64)-counts
    assert np.array_equal(L,L.T) and np.array_equal(L.sum(axis=0),np.zeros(240))
    single=[];pair=[];pair_calls=0
    for v in bases:
        energy=F(int(v@(L@v)),int(v@v));res=L@v-float(energy)*v
        assert np.array_equal(res,np.zeros(240));single.append(str(energy))
    for v in bases[:2]:
        for w in bases[2:]:
            initial=np.outer(v,w);mean=np.zeros((240,240),dtype=np.int64);vacuum={}
            for p,n in zip(*np.nonzero(initial)):
                weight=int(initial[p,n])
                for e,(a,b)in enumerate(edges):
                    for coin in (0,1):
                        initial_source=[0]*240;initial_source[p]=1;initial_source[n]=-1
                        actual,record=native(initial_source,(a,b),coin);j=record[1];pair_calls+=1
                        out,k=lf.pair_endpoint(int(p),int(n),a,b,coin)
                        expected=[0]*240
                        if out is not None:expected[out[0]]=1;expected[out[1]]=-1
                        assert tuple(expected)==actual and j==k
                        if out is None:vacuum[(e,coin,j)]=vacuum.get((e,coin,j),0)+weight
                        else:mean[out]+=weight
            generator=1260*initial-mean
            energy=F(int(np.sum(initial*generator)),int(np.sum(initial*initial)))
            assert np.array_equal(generator,float(energy)*initial) and not any(vacuum.values())
            pair.append({'generator':str(energy),'full_charged_residual':0,'vacuum_fine_amplitudes':list(vacuum.items())})
    for a,b in edges:
        for coin in (0,1):
            out,record=native([0]*240,(a,b),coin)
            assert out==(0,)*240 and record==(coin,0)
    return L,{'literal_native_single_calls':calls,'literal_native_pair_mode_calls':pair_calls,
              'literal_native_vacuum_calls':2*len(edges),
              'copied_pinned_native_kernel_sha256':sha(HERE/'inputs/native_event.py'),
              'all_formula_comparisons_pass':True,'vacuum_zero_columns':2*len(edges),
              'single_generators':single,'pair_generators':pair,'all240_modes_retained':True}

def certificate(initial,program,steps,pair,energy,time):
    one=lf.step(initial,program,pair);target=lf.phase(-energy*time/steps)*initial
    residual=float(np.linalg.norm((one-target).reshape(-1,initial.shape[-1]),'fro'))
    # Per matching: disjoint two-point mixing, twice for pairs, plus a
    # diagonal contact overwrite. Directed sine/cosine and complex floating
    # arithmetic each receive a conservative 512u relative Frobenius guard.
    accounting=rounding.guard(program,steps,pair,energy,time)
    guard=float(F(accounting['one_step_state_arithmetic_upper']))
    total=steps*(residual+guard)
    return {'one_step_Frobenius_residual_centre':residual,'one_step_rounding_upper':guard,
            'directed_rounding':accounting,
            'repetitions':steps,'selected_subspace_operator_error_upper':total}

def read(initial,final):
    cols=initial.shape[-1];z=np.vdot(initial,final)/cols
    # Complete ideal-reference record centres; physical native pointer,
    # processor and loader outcomes are retained by the comparison radius.
    return {'fringe_centre':[float(z.real),float(z.imag)],
            'X_centres':[(1+float(z.real))/2,(1-float(z.real))/2],
            'Y_centres':[(1+float(z.imag))/2,(1-float(z.imag))/2],
            'final_norm_centre':float(np.vdot(final,final).real/cols)}

def decode(short,long,epsilon):
    # No expected energy appears here. Complete law error bounds both
    # pointer means and retains failed outcomes in the denominator.
    radius=2*math.sqrt(2)*epsilon/(1-epsilon)
    assert radius<1
    phase_error=radius/(1-radius)
    z=complex(*short['fringe_centre']);e=-math.atan2(z.imag,z.real)/.1
    coarse=[e-phase_error/.1,e+phase_error/.1]
    raw=-math.atan2(long['fringe_centre'][1],long['fringe_centre'][0])
    winding_lo=math.ceil((4*coarse[0]-raw-phase_error)/(2*math.pi))
    winding_hi=math.floor((4*coarse[1]-raw+phase_error)/(2*math.pi))
    assert winding_lo==winding_hi
    fine=(raw+2*math.pi*winding_lo)/4
    return {'short_interval':coarse,'winding_from_short':winding_lo,
            'energy_interval':[fine-phase_error/4,fine+phase_error/4]}

def main():
    started=time.monotonic();reg=json.loads((HERE/'PREREG.json').read_text())
    for name,h in reg['files'].items():assert sha(HERE/name)==h,name
    config=reg['configuration'];edges=[tuple(e)for e in json.loads((HERE/'inputs/graph.json').read_text())['seams']]
    cert=json.loads((HERE/'inputs/certificate.json').read_text());basis=np.array(cert['integer_modes'],dtype=np.int64)
    L,census=source_census(edges,basis);local=lf.local_pair_census();routes=lf.route_audit(edges)
    local_calls=0
    for p,n,coin,current,destination in local['fine_records']:
        q=[0]*4;q[p]=1;q[n]=-1;out,record=native(q,(0,1),coin)
        expected=[0]*4
        if destination is not None:expected[destination[0]]=1;expected[destination[1]]=-1
        assert out==tuple(expected)and record==(coin,current);local_calls+=1
    census['literal_native_local_type_calls']=local_calls
    source0=np.column_stack((basis[0]/2,(2*basis[1]-basis[0])/math.sqrt(12)))
    source1=np.column_stack((basis[2]/2,(2*basis[3]-basis[2])/math.sqrt(12)))
    pair=np.stack([np.outer(source0[:,i],source1[:,j])for i in range(2)for j in range(2)],axis=-1)
    assert np.max(np.abs(pair.reshape(-1,4).conj().T@pair.reshape(-1,4)-np.eye(4)))<1e-14
    energies={'single':float(F(census['single_generators'][0])),
              'pair':float(F(census['pair_generators'][0]['generator']))}
    programs={};descriptions={};hashes={}
    for t in config['times']:
        program,desc,h=lf.description(edges,t,config['steps']);programs[str(t)]=program;descriptions[str(t)]=desc;hashes[str(t)]=h
    catalog,matrices=nc.execute([h for p in programs.values()for _,h in p])
    b=F(catalog['complete_per_pulse_allowance']);words={};native_error={}
    for name,is_pair in [('single',False),('pair',True)]:
        words[name]=lf.word_counts(edges,config['steps'],is_pair,catalog['ccx_pulses'])
        native_error[name]=nc.complete(words[name]['total_pulses'],b)
    log({'source_census':census,'native_catalog':catalog['executed_native_M6_matrices'],
         'per_pulse':float(b),'represented_pulses':{k:v['total_pulses']for k,v in words.items()}})
    records={};certificates={};arrays={}
    for t in config['times']:
        records[str(t)]={};certificates[str(t)]={};program=programs[str(t)]
        for name,initial,is_pair in [('single',source0,False),('pair',pair,True)]:
            c=certificate(initial,program,config['steps'],is_pair,energies[name],t)
            log({'time':t,'object':name,'one_step_certificate':c})
            final=lf.propagate(initial,program,config['steps'],is_pair,
                     lambda n:log({'time':t,'object':name,'completed_macrosteps':n}))
            records[str(t)][name]=read(initial,final);certificates[str(t)][name]=c
            arrays[f'{name}_{str(t).replace(".","p")}']=final
    # Same generic source programme and same reader, with one actual edge
    # disabled. It changes the operator on the source, not a target energy.
    changed_edges=[e for e in edges if set(e)!={4,5}];assert len(changed_edges)==629
    changed,changed_desc,changed_hash=lf.description(changed_edges,4,config['steps'])
    changed_final=lf.propagate(source0,changed,config['steps'],False)
    changed_record=read(source0,changed_final)
    changed_distance=abs(complex(*changed_record['fringe_centre'])-complex(*records['4']['single']['fringe_centre']))
    # A contact input probes the nonfactorizing killed charged block.
    contact=np.zeros((4,4,1),complex);contact[0,1,0]=1
    actual=contact.copy();lf.matching(actual,np.array([[0,1]]),.25,True)
    single_edge=np.eye(4,dtype=complex);lf.matching(single_edge,np.array([[0,1]]),.25,False)
    unrestricted=single_edge@contact[:,:,0]@single_edge.T
    contact_difference=float(np.linalg.norm(actual[:,:,0]-unrestricted))
    # Native local census/exponential equality includes all twelve charged
    # local occupancy types, not only the selected compact sector.
    local_columns=np.zeros((4,4,12),complex)
    for i,(a,c)in enumerate(local['states']):local_columns[a,c,i]=1
    local_out=local_columns.copy();lf.matching(local_out,np.array([[0,1]]),.25,True)
    local_matrix=np.array([local_out[a,c]for a,c in local['states']])
    local_residual=float(np.linalg.norm(local_matrix-expm(-.25j*np.array(local['generator'])),'fro'))
    assert local_residual<1e-13
    inherited=F(cert['common_full_preparation_and_apparatus_allowance'])
    reset=F(35,2147483648);pointer=F(200000000001,20000000000000000000000)
    budget={};decoded={}
    for name in ('single','pair'):
        numerical=max(certificates[str(t)][name]['selected_subspace_operator_error_upper']for t in config['times'])
        comparison=F(str(numerical));apparatus=words[name]['initial_and_terminal_scratch_bits']*(reset+2*pointer)+2*pointer+4*b
        stored=max(2*F(certificates[str(t)][name]['directed_rounding']['stored_full_propagation_arithmetic_upper'])+
                   F(certificates[str(t)][name]['directed_rounding']['final_inner_product_readout_upper'])for t in config['times'])
        s=inherited+native_error[name]+comparison+stored+apparatus
        total=s/(1-2*s);budget[name]={'inherited':str(inherited),'source_word':str(native_error[name]),
                 'selected_source_flow':str(comparison),'stored_propagation_and_readout':str(stored),
                 'reset_read':str(apparatus),'complete':str(total)}
        decoded[name]=decode(records['0.1'][name],records['4'][name],float(total))
    one=decoded['single']['energy_interval'];two=decoded['pair']['energy_interval']
    ratio=[two[0]/one[1],two[1]/one[0]];kev=[511*x for x in ratio]
    lmat=np.array(local['generator'])
    gates={'source_local_hermitian_and_complete':bool(np.array_equal(lmat,lmat.T)and np.array_equal(lmat@lmat,2*lmat)
                                                   and len(local['fine_records'])==24 and len(local['vacuum_records'])==4),
           'local_full12_column_match':local_residual<1e-13,
           'both_complete_errors_below_001':all(F(v['complete'])<F(1,1000)for v in budget.values()),
           'pair_interval_contains2':ratio[0]<2<ratio[1],
           'pair_interval_excludes_old1024_45085':kev[1]<1024.45085,
           'changed_seam_same_reader_separates':changed_distance>.05,
           'contact_distinguishes_unrestricted_tensor':contact_difference>.1,
           'vacuum_zero_by_source_action':census['vacuum_zero_columns']==1260,
           'full_coupled_reflection_checked':routes['complete_coupled_pair_reflection_columns']==2*240*239*630}
    result={'candidate':reg['candidate'],'prereg_sha256':sha(HERE/'PREREG.json'),
      'source_census':census,'local_charged_census':local,'route_audit':routes,'native_catalog':catalog,'words':words,
      'source_program_hashes':hashes,'source_programs':descriptions,'certificates':certificates,
      'complete_unfiltered_reference_record_centres':records,'complete_budgets':budget,
      'decoder_intervals':decoded,'pair_to_single_interval':ratio,'pair_keV_at511_anchor':kev,
      'changed_seam':{'program_sha256':changed_hash,'record':changed_record,'fringe_separation':changed_distance},
      'contact_tensor_difference':contact_difference,'local_full12_column_residual':local_residual,
      'outcomes':['native source position and reference bit','initial collision/invalid/vacuum type',
                  'inherited preparation/route/load failure','processor first failure at full prefix',
                  'dirty scratch','source padding export failure','reset timeout','pointer outcomes'],
      'execution_semantics':'Fresh M6 primitive catalogue and literal source census; exact clean-address routing compression plus matching product propagation, with complete native prefix error. No new full raw multibillion-pulse replay claimed.',
      'gates':gates,'all_gates_pass':all(gates.values()),'elapsed_seconds':time.monotonic()-started}
    out=HERE/'attempts'/config['attempt'];out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/'SOURCE_OUTPUTS.npz',**arrays)
    np.savez_compressed(out/'NATIVE_M6.npz',**{str(i):x for i,x in enumerate(matrices.values())})
    result['matrix_catalog_keys']=list(matrices);exit_code=finish(out/'RECEIPT.json',result)
    log({'gates':gates,'pair_keV':kev,'seconds':result['elapsed_seconds']})
    raise SystemExit(exit_code)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        write(HERE/'LAST_EXECUTION_FAILURE.json',{'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()})
        raise
