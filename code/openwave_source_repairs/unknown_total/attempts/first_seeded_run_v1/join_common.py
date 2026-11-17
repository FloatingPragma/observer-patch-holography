"""Consume actual unknown-total handoffs with the sealed common pump maps."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
import numpy as np
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def pump(Q,start):
    rho=np.zeros((12,12),complex);rho[start,start]=1
    escape=np.eye(12)-Q.conj().T@Q
    for _ in range(512):rho=Q@rho@Q.conj().T+np.trace(escape@rho)*np.eye(12)/12
    return rho

def calculate():
    origin=read(HERE/'CONSUMER_ORIGINS.json')
    assert sha(HERE/'inputs/common_MANIFEST.sha256')==origin['common_manifest_sha256']
    manifest={n:h for h,n in (line.split('  ',1)for line in(HERE/'inputs/common_MANIFEST.sha256').read_text().splitlines())}
    for name,record in origin['copied'].items():assert sha(HERE/'inputs'/name)==record['sha256']==manifest[record['common_relative']]
    h=read(HERE/'HANDOFFS.json');prior=read(HERE/'inputs/common_RECEIPT.json');own=read(HERE/'RECEIPT.json')
    assert h['producer_records_sha256']==own['records_sha256']
    seeds=[]
    for trial in h['handoffs']:
        q=trial['observed_source'];assert q==trial['actual_source']
        pos=[i for i,x in enumerate(q)if x==1];neg=[i for i,x in enumerate(q)if x==-1]
        assert len(pos)==len(neg)==1 and sum(abs(x)for x in q)==2
        assert pos==[0] and neg==[14] and trial['positive_frame']==[0,1]
        seeds.append((pos[0],neg[0]-12))
    assert len(seeds)==24 and set(seeds)=={(0,2)}
    m=np.load(HERE/'inputs/common_maps.npz')
    maps=[np.load(HERE/'inputs'/name)['certified_Q_proposal']for name in ['first_maps.npz','second_maps.npz']]
    rho=[pump(Q,start)for Q,start in zip(maps,seeds[0])]
    difference=max(float(np.max(abs(rho[0]-m['first_final_density']))),float(np.max(abs(rho[1]-m['second_final_density']))))
    assert difference<1e-12
    flags=[];dirty=[];scratch=[]
    for i,state in enumerate(rho):
        flags.append([float(np.trace(e[:12,:12]@state).real)for e in m[f'loader{i}_effects']])
        dirty.append(float(np.trace(m[f'loader{i}_dirty_gram'][:12,:12]@state).real))
        scratch.append([float(np.trace(g[:12,:12]@state).real)for g in m[f'loader{i}_scratch_word_grams']])
    joint_flags=np.outer(*flags)
    residual=float(1-joint_flags.sum()-dirty[0]-sum(flags[0])*dirty[1])
    assert abs(residual)<1e-7
    eta=F(own['preparation_inconclusive_upper'])
    energy=F(prior['energy_complete_error_upper'])+eta
    spatial=F(prior['spatial_halfturn_TV_upper'])+2*eta
    control=F(prior['active_dephasing_TV_lower'])-2*eta
    # The original admitted-state fringe centers are unchanged. For eps<.001,
    # derivative of the phase-radius/(4*timeunit) long-energy half-width is <1.
    # Four eta plus explicit 1e-12 outward float slack is conservative.
    widening=4*float(eta)+1e-12
    intervals={k:[v['long'][0]-widening,v['long'][1]+widening]for k,v in prior['energy_intervals'].items()}
    ratio=[intervals['pair'][0]/intervals['single'][1],intervals['pair'][1]/intervals['single'][0]]
    result=dict(candidate=own['candidate'],common_candidate='A-old240-selected-compact-source-v1',
        common_manifest_sha256=origin['common_manifest_sha256'],actual_handoffs_consumed=24,
        unique_actual_basis_inputs=1,compressed_pump_epochs_executed=1024,
        terminal_density_max_difference=difference,region_flag_centres=joint_flags.tolist(),
        local_dirty_centres=dirty,all16_local_scratch_centres=scratch,
        unresolved_processor_plus_rounding_residual=residual,
        additional_preparation_error_upper=str(eta),complete_energy_error_upper=str(energy),
        positive_position_halfturn_TV_upper=str(spatial),active_dephasing_TV_lower=str(control),
        pair_to_single_interval=ratio,pair_keV_at511_anchor_interval=[511*x for x in ratio],
        no_unit_prefix_preserved=True,all_failure_prefixes_preserved=True,
        literal_new_native_calls=0,new_fullflow_executions=0,full_row_award=False,other_group_G6=False)
    assert energy<F(1,1000) and spatial<F(1,1000) and control>F(35,1000)
    assert ratio[0]<2<ratio[1]
    return result,dict(first_density=rho[0],second_density=rho[1])

def main():
    p=argparse.ArgumentParser();p.add_argument('--record',action='store_true');a=p.parse_args()
    for n,h in read(HERE/'CONSUMER_PREREG.json')['sha256'].items():assert sha(HERE/n)==h,n
    result,arrays=calculate()
    if a.record:
        with(HERE/'CONSUMER_RESULT.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
        np.savez_compressed(HERE/'CONSUMER_MAPS.npz',**arrays)
    else:
        assert read(HERE/'CONSUMER_RESULT.json')==result
        old=np.load(HERE/'CONSUMER_MAPS.npz')
        for k,v in arrays.items():assert np.array_equal(old[k],v)
    print(json.dumps({'actual_handoffs':24,'density_difference':result['terminal_density_max_difference'],
        'energy_error':float(F(result['complete_energy_error_upper'])),
        'field_TV':float(F(result['positive_position_halfturn_TV_upper'])),
        'pair_keV':result['pair_keV_at511_anchor_interval']}))
if __name__=='__main__':main()
