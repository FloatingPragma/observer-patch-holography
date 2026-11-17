#!/usr/bin/env python3
"""Independent final bookkeeping, detector and conditioning review."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,time,sys
import numpy as np
sys.set_int_max_str_digits(0) # All rational input lengths are bounded by frozen receipts.
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2]
P=HERE.parents[1]/'work/2026-11-15_codex_source_repair_v2'
def load(name):return json.loads((P/name).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def phase_bound(e):
    radius=2*math.sqrt(2)*float(e)/(1-2*float(e))
    return radius/(1-radius)
started=time.perf_counter()
c=load('coupled_energy/RECEIPT.json')['scientific'];r=load('routed_join/RECEIPT.json')['scientific']
assert r['coupled_receipt_sha256']==sha(P/'coupled_energy/RECEIPT.json')
source=np.load(P/'poisson_energy/source_sector_maps.npz');centers=np.load(P/'coupled_energy/maps.npz')
maximum_pump_error=0.;upstream=[];flags=[]
for name,folder,offset,initial in [('first','first_cell_certified',0,0),('second','second_cell',12,2)]:
    maps=np.load(P/folder/'maps.npz');Q=maps['certified_Q_proposal'];E=np.eye(12)-Q.conj().T@Q
    S=np.kron(Q,Q.conj())+np.outer((np.eye(12)/12).ravel(),E.T.ravel())
    v=np.zeros(144);v[13*initial]=1
    rho=(np.linalg.matrix_power(S,512)@v).reshape(12,12)
    maximum_pump_error=max(maximum_pump_error,float(np.max(abs(rho-centers['source_pump_'+name]))))
    assert maximum_pump_error<1e-11
    f=[float(np.trace(effect@rho).real) for effect in source[name+'_effects']]
    assert np.max(abs(np.asarray(f)-centers['source_flags_'+name]))<1e-11
    flags.append(f);upstream.append(float(np.trace(source[name+'_processor_failure']@rho).real))
assert abs(upstream[0]-c['upstream_sector_processor_failure_centres']['first'])<1e-20
assert abs(upstream[1]-c['upstream_sector_processor_failure_centres']['second'])<1e-20
fp=upstream[0]+sum(flags[0])*upstream[1]
assert abs(fp-c['upstream_sector_processor_failure_centres']['pair_first_failure'])<1e-20
for timekey,group in c['complete_source_coupled_records'].items():
    for name,record in group.items():
        meta=c['native_flag_controlled_words'][timekey][name]
        f=np.array(flags[0] if name=='single' else np.kron(flags[1],flags[0]))
        for axis in ('X','Y'):
            law=f@np.array(meta['complete_per_input_flag_output_centres'][axis])
            target=record[axis];assert np.max(abs(law-target['diagnostic_probability_centres']))<1e-11
            assert target['upstream_sector_processor_failure_centre']>=0
            assert abs(1-sum(law)-target['upstream_sector_processor_failure_centre']-
                       target['downstream_processor_failure_centre']-target['arithmetic_normalization_residual'])<1e-11

# The source flag is discarded after its actual pointer record. Its post-read
# source state therefore depends only on the flag POVM effect; this gives a
# different implementation of the producer's local-pointer Kraus contraction.
base=OW/'shared/compact_logarithmic_flow/attempts/stage2_002'
old=json.loads((base/'RECEIPT.json').read_text())
archive=base/'evidence/reference_preparation.npz'
assert sha(archive)==old['artifacts']['evidence/reference_preparation.npz']
effects=np.load(archive)['pointer_effects_C']
witness=centers['negative_source_sector_columns'].reshape(2,256,3)
expected=c['negative_unresolved_native_address_records'];address_error=0.
for j in range(3):
    A=witness[:,:,j];rho=A.T@effects[0].T@A.conj();law=[]
    for address in range(256):
        M=np.ones((1,1),complex)
        for bit in reversed(range(8)):M=np.kron(M,effects[(address>>bit)&1])
        law.append(float(np.sum(M.T*rho).real))
    address_error=max(address_error,float(np.max(abs(np.array(law)-expected[j]['complete_flag0_address_centres']))))
    assert address_error<1e-12
    category=['charged_unresolved','physical_vacuum','invalid_padding'][j]
    assert expected[j]['measured_class_centres'][category]>1-1e-7

reset=F(35,2**31);pointer=F(200000000001,20000000000000000000000);b=F(1,50000000000)
for row in c['additional_coupled_costs']:
    N=row['pulses'];x=N*b;e=x/(1-x);failure=N*b*b/(1-x)**2
    bound=(e*(2+e)+failure)/2
    assert bound==F(row['pulse_bound'])
    assert F(row['reset_bound'])==row['reset_bits']*reset
    assert F(row['pointer_bound'])==row['pointer_reads']*pointer
    assert F(row['complete_bound'])==bound+row['reset_bits']*reset+row['pointer_reads']*pointer
assert sorted(set(row['pulses'] for row in c['additional_coupled_costs']))==[69,151]
assert F(c['unresolved_address_read_error_upper'])==8*(reset+pointer)
for linear,gain in zip(c['linear512_preparation_word_and_apparatus_bounds'],c['gain_corrected512_preparation_word_and_apparatus_bounds']):
    linear=F(linear);assert F(gain)==linear/(1-2*linear)
prep=list(map(F,c['complete_source_preparation_error_upper']));reader=F(c['complete_coupled_sector_reader_error_upper'])
route=1-(1-F(1,2**32))**8;assert route==F(r['route_timeout_upper'])
S=sum(prep)+reader+route;assert S==F(r['linear_physical_TP_composition_bound'])
joint=S/(1-2*S);assert joint==F(r['complete_join_error_upper']) and joint<F(1,1000)
singleS=prep[0]+reader+route;single=singleS/(1-2*singleS)
# Exact product-of-gains composition is tighter than the filed geometric guard.
product=((1+2*prep[0])*(1+2*prep[1])*(1+2*reader)-1)/2+route
assert joint-product>F(1,10**7)
intervals={};center_bias=0.
for name,error in [('single',single),('pair',joint)]:
    bound=phase_bound(error);record=r['handoffs'][0]['complete_energy_output_centres']
    short=complex(*record['1/10'][name]['complex_accepted_fringe'])
    long=complex(*record['4'][name]['complex_accepted_fringe'])
    center=-10*math.atan2(short.imag,short.real);shortinterval=[center-10*bound,center+10*bound]
    basephase=-math.atan2(long.imag,long.real)/4
    ks=[k for k in range(32) if shortinterval[0]-bound/4<=basephase+k*math.pi/2<=shortinterval[1]+bound/4]
    assert len(ks)==1 # Resolves aliases from the short measurement, not6 or12.
    value=basephase+ks[0]*math.pi/2;intervals[name]=[value-bound/4,value+bound/4]
    expected=r['energy_intervals'][name]
    assert expected['long_winding_from_short_interval']==ks[0]
    assert np.max(abs(np.array(expected['long_interval'])-intervals[name]))<1e-12
    assert abs(expected['acceptance_floor']-(1-float(error)))<1e-15
    # Direct central-array comparison to the exact ideal-source law. The
    # 1e-7 outward guard exceeds all observed arithmetic/trig roundoff and is
    # smaller than the filed global gain margin, so no negative center is
    # interpreted as a negative physical probability.
    energy=6 if name=='single' else 12
    for t in ('1/10','4'):
        m=c['native_flag_controlled_words'][t][name];p0,p1=m['success_outputs']
        for axis in ('X','Y'):
            ideal=np.zeros(2**m['bits']);value=(math.cos if axis=='X' else math.sin)(-float(F(t))*energy)
            ideal[p0]=(1+value)/2;ideal[p1]=(1-value)/2
            item=c['complete_source_coupled_records'][t][name][axis]
            bias=sum(abs(np.array(item['diagnostic_probability_centres'])-ideal))/2
            bias+=(abs(item['upstream_sector_processor_failure_centre'])+abs(item['downstream_processor_failure_centre']))/2
            center_bias=max(center_bias,float(bias));assert bias<F(1,10**7)
lo,hi=intervals['single'];pl,ph=intervals['pair'];ratio=[pl/hi,ph/lo]
assert np.max(abs(np.array(ratio)-r['pair_to_single_interval']))<1e-12
assert ratio[0]<2<ratio[1] and 511*ratio[1]<1024.45084870234
assert len(r['handoffs'])==12
for hand in r['handoffs']:
    assert hand['actual_routed_pair_code']==14
    assert hand['actual_source']==[1]+[0]*13+[-1]+[0]*225
    assert hand['native_frame']
result={'same_group_only':True,'producer_modules_imported':[],
 'source_pump_superoperator_error':maximum_pump_error,'full_negative_address_reader_error':address_error,
 'upstream_sector_failure_retained':True,'exact_cost_and_gain_arithmetic':True,
 'complete_join_error_upper':str(joint),'product_gain_error_upper':str(product),
 'gain_margin_over_product_bound':str(joint-product),'maximum_signed_center_bias':center_bias,
 'principal_short_and_unique_long_decoder':True,'energy_intervals':intervals,'pair_ratio_interval':ratio,
 'pair_keV_interval':[511*x for x in ratio],'old_excess_excluded':True,
 'accepted_preparation_loading_boundary_required':True,
 'bindings':{str(p.relative_to(P)):sha(p) for p in [P/'run_coupled_energy.py',P/'coupled_energy/RECEIPT.json',P/'coupled_energy/maps.npz',P/'join_routed_source.py',P/'routed_join/RECEIPT.json',P/'routed_join/maps.npz',P/'PROOF.md']},
 'master_bindings':{name:sha(OW/name) for name in ['CAUTIONS.md','SUBMISSION_CHECKLIST.md']},
 'seconds':time.perf_counter()-started}
target=HERE/'FINAL_RESULT.json'
if target.exists():
    prior=json.loads(target.read_text())
    assert {k:v for k,v in result.items() if k!='seconds'}=={k:v for k,v in prior.items() if k!='seconds'}
else:
    with target.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps({k:result[k] for k in ['full_negative_address_reader_error','principal_short_and_unique_long_decoder','pair_keV_interval','old_excess_excluded','seconds']}))
