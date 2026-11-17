#!/usr/bin/env python3
"""Independent common-source review. No producer module is imported."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,math,sys,time
import numpy as np
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent;P=HERE.parent/'selected_compact_source_join';I=P/'inputs'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gate_apply(state,matrix,targets):
    # Tensor-axis implementation independent of producer's bit-index gather.
    axes=[12-bit for bit in reversed(targets)]
    rest=[i for i in range(13) if i not in axes];order=axes+rest+[13]
    tensor=state.reshape((2,)*13+(state.shape[1],)).transpose(order)
    changed=(matrix@tensor.reshape(len(matrix),-1)).reshape(tensor.shape)
    return changed.transpose(np.argsort(order)).reshape(state.shape)
def spectral_power(q,initial,dephase=False):
    E=np.eye(12)-q.conj().T@q
    S=np.kron(q,q.conj())+np.outer((np.eye(12)/12).ravel(),E.T.ravel())
    if dephase:S[np.arange(144)%13!=0]=0
    x=np.zeros(144);x[13*initial]=1
    return(np.linalg.matrix_power(S,512)@x).reshape(12,12)
def word_bound(n):
    b=F(1,50000000000);x=n*b;e=x/(1-x);f=n*b*b/(1-x)**2+e*e
    return(e*(2+e)+f)/2
started=time.perf_counter();r=read(P/'RECEIPT.json');assert r['all_gates_pass']
reg=read(P/'PREREG.json')
for n,h in reg['sha256'].items():assert sha(P/n)==h,n
provenance=read(P/'PROVENANCE.json');manifest_maps={}
for name,c in provenance['components'].items():
    m=P/c['copied_manifest'];assert sha(m)==c['manifest_sha256']
    manifest_maps[name]={n:h for h,n in(line.split('  ',1) for line in m.read_text().splitlines())}
for n,binding in provenance['inputs'].items():
    assert sha(P/n)==binding['sha256']==manifest_maps[binding['component']][binding['source_relative_path']]

# One-step signed action: derive destinations by occupancy cases, without
# reusing the producer's dictionary balancing formula.
loader=json.loads(gzip.decompress((I/'loader/RECORDS.json.gz').read_bytes()))
edges=read(I/'loader/graph.json')['seams'];count=0
for code,e,c,j,out in loader['source_join']['loaded_branches']:
    p,n=divmod(code,256);a,b=edges[e];hp=p in(a,b);hn=n in(a,b)
    if hp and hn:want=0;current=1 if p==a else -1
    elif hp:
        pp=a if c else b;want=256*pp+n;current=int(p==a)-int(pp==a)
    elif hn:
        nn=b if c else a;want=256*p+nn;current=-int(n==a)+int(nn==a)
    else:want=code;current=0
    assert(want,current)==(out,j);count+=1
assert count==181440

maps=np.load(P/'INTERFACE_MAPS.npz');app=np.load(I/'energy/apparatus/maps.npz')
local={'ccx':app['native_routine_ccx'],'flip':app['native_routine_flip']}
pointer=np.load(I/'flow/attempts/stage2_002/evidence/reference_preparation.npz')['pointer_effects_C']
errors=[];broad=F()
reset=F(35,2**31);ptr=F(200000000001,20000000000000000000000)
for k,m in enumerate(r['new_broad_loader_programs']):
    cell=m['cell'];labels=m['labels'];state=np.zeros((8192,len(labels)),complex)
    for j,label in enumerate(labels):state[label,j]=1
    truth=np.arange(512)
    for kind,t in m['operations']:
        state=gate_apply(state,local[kind],t)
        if kind=='flip':truth^=1<<t[0]
        else:truth^=(((truth>>t[0])&1)&((truth>>t[1])&1))<<t[2]
    want=np.array([i^(256 if cell<=(i&255)<cell+12 else 0) for i in range(512)])
    assert np.array_equal(truth,want)
    clean=state[:512].reshape(2,256,len(labels))
    grams=np.array([state[512*s:512*(s+1)].conj().T@state[512*s:512*(s+1)] for s in range(16)])
    effects=np.array([sum(pointer[out,i,j]*clean[i].conj().T@clean[j] for i in range(2) for j in range(2)) for out in range(2)])
    residual=max(float(np.max(abs(clean-maps['loader%d_clean_columns'%k]))),float(np.max(abs(grams-maps['loader%d_scratch_word_grams'%k]))),float(np.max(abs(effects-maps['loader%d_effects'%k]))))
    errors.append(residual);assert residual<1e-12
    N=m['represented_native_pulses'];assert N==sum(40 if name=='ccx' else 2 for name,t in m['operations'])
    b=F(1,50000000000)
    assert F(m['processor_first_failure_probability_upper'])==N*b*b/(1-N*b)**2
    assert m['blank_reset_bits_charged']==10 and m['reset_pointer_reads_charged']==10 and m['flag_and_scratch_pointer_reads_charged']==5
    broad+=word_bound(N)+10*reset+15*ptr
assert broad==F(r['complete_broad_loader_error_upper'])

certs=[read(I/'energy'/name/'RECEIPT.json')['scientific'] for name in('first_cell_certified','second_cell')]
Qs=[np.load(I/'energy'/name/'maps.npz')['certified_Q_proposal'] for name in('first_cell_certified','second_cell')]
assert np.array_equal(Qs[0],np.load(I/'spatial/inputs/first_cell_maps.npz')['certified_Q_proposal'])
states=[spectral_power(Qs[0],0),spectral_power(Qs[1],2),spectral_power(Qs[0],0,True)]
for state,name in zip(states,['first_final_density','second_final_density','dephasing_final_density']):
    assert np.max(abs(state-maps[name]))<1e-11
field=np.load(I/'spatial/FIELD_A_RECORDS.npz')
assert np.max(abs(states[0]-field['native_density_history'][-1]))<1e-11
assert np.max(abs(states[2]-field['dephasing_density_history'][-1]))<1e-11

# Reconstruct the actual broad-loaded positive-position marginal. Region1 is
# a spectator acceptance factor; neither cell is filtered by the dark sector.
largest_record_error=0.
for name,first in [('native',states[0]),('dephasing',states[2])]:
    conditional=[]
    for k,rho in enumerate((first,states[1])):
        blocks=maps['loader%d_clean_columns'%k][:,:,:12]
        values=[]
        for out in range(2):
            values.append(sum(pointer[out,i,j]*blocks[j]@rho@blocks[i].conj().T for i in range(2) for j in range(2)))
        conditional.append(values)
    density=conditional[0][1]*np.trace(conditional[1][1]).real
    law=[]
    for address in range(256):
        E=np.ones((1,1),complex)
        for bit in reversed(range(8)):E=np.kron(E,pointer[(address>>bit)&1])
        law.append(float(np.sum(E.T*density).real))
    err=float(np.max(abs(np.array(law)-r['spatial_full_outcome_centres'][name]['all256_position_centres'])))
    largest_record_error=max(largest_record_error,err);assert err<1e-11

flow=F(read(I/'energy/SHARP_FLOW_CERTIFICATE.json')['512_complete_U4_halfdiamond_upper'])
raw=[flow+F(c['512epoch_apparatus_error_upper']) for c in certs]
gain=[x/(1-2*x) for x in raw];bad=[F(c['512epoch_ideal_bad_weight_upper']) for c in certs]
route=F(read(I/'route/RECEIPT.json')['per_trial_timeout_upper']);terminal=8*(reset+ptr)
S=sum(gain)+sum(bad)+route+broad+terminal;ferr=S/(1-2*S)
assert 2*ferr==F(r['spatial_halfturn_TV_upper']) and 2*ferr<F(1,250)
fa=read(I/'spatial/FIELD_A_RESULT.json');comparison=F(fa['comparison_512_error_upper'])
delta=F(certs[1]['certified_subblock_error_upper']);kap=4*delta+2*delta*delta;u=F(1,2**53);gamma=128*u/(1-128*u)
second=512*kap/(1-512*kap)+512*12*4096*gamma/(1-512*kap)+F(1,10**12)
x=raw[0]+512*terminal;C=x/(1-2*x)+gain[1]+route+broad+terminal
control=F(fa['dephasing_coordinate_halfturn_TV'])-2*(C/(1-2*C)+comparison+second)
assert control==F(r['active_dephasing_TV_lower']) and control>F(1,100)
energy=read(I/'energy/routed_join/RECEIPT.json')['scientific']
S=F(energy['linear_physical_TP_composition_bound'])+broad
assert S/(1-2*S)==F(r['energy_complete_error_upper'])<F(1,1000)
assert r['pair_to_single_interval'][0]<2<r['pair_to_single_interval'][1]<r['old_fixed_count_ratio_retained']
assert r['broad_loader_before_assay_choice'] and r['spatial_assay_before_dark_sector_diagnosis']
result={'same_group_only':True,'producer_modules_imported':[],
 'signed_action_columns_recomputed':count,'independent_loader_tensor_errors':errors,
 'independent_full256_positive_position_error':largest_record_error,
 'same_preparation_operator_and_density':True,'full_cell_loader_no_dark_filter':True,
 'scratch_outcomes':16,'processor_first_failure_bounds_explicit':True,
 'broad_loader_error_upper':str(broad),'spatial_halfturn_TV_upper':str(2*ferr),
 'active_dephasing_TV_lower':str(control),'energy_complete_error_upper':r['energy_complete_error_upper'],
 'pair_keV_interval':r['pair_keV_at_single511_anchor_interval'],
 'bindings':{name:sha(P/name) for name in ['PREREG.json','interface.py','RECEIPT.json','INTERFACE_MAPS.npz','AMENDMENT_V2.md','PROVENANCE.json']},
 'seconds':time.perf_counter()-started}
target=HERE/'RESULT.json'
if target.exists():
    old=read(target);assert {k:v for k,v in result.items() if k!='seconds'}=={k:v for k,v in old.items() if k!='seconds'}
else:
    with target.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps({k:result[k] for k in ['signed_action_columns_recomputed','independent_loader_tensor_errors','independent_full256_positive_position_error','pair_keV_interval','seconds']}))
