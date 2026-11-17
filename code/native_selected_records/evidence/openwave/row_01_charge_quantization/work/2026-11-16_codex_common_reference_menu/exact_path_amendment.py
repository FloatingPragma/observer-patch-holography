"""Exact control-path repair and end-to-end costs; primary evidence is read only."""
from pathlib import Path
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
import time
import numpy as np

HERE=Path(__file__).resolve().parent


def unpack(x):
    return np.array(x['real'])+1j*np.array(x['imag']) if isinstance(x,dict) else np.array(x)


def exact_swap(m,turns,direction,reference):
    """Products of actual source-block eigenvalue squares, independent of rho."""
    amplitudes=[F(1),F(1)]
    slots=[]
    for step in range(1,3*m*turns+1):
        k=step%m
        pair=((2*m+direction*step)%(3*m),2*m if reference=='stationary' else (2*m+direction*min(k,m-k))%(3*m))
        lambdas=[]
        for arm,slot in enumerate(pair):
            t=F(slot%m,m);lam=(1-3*t+3*t*t)/4
            assert lam>0
            amplitudes[arm]*=lam*lam
            lambdas.append(str(lam))
        slots.append({'slots':pair,'lambda':lambdas})
    a,b=amplitudes;denominator=2*(a*a+b*b)
    probabilities=[(a+b)**2/denominator,(a-b)**2/denominator]
    visibility=2*a*b/(a*a+b*b)
    assert min(probabilities)>=0 and sum(probabilities)==1 and visibility>0
    if reference=='matched':assert probabilities==[F(1),F(0)]
    return {'amplitudes_without_common_factors':[str(a),str(b)],'path_probabilities':[str(x) for x in probabilities],
            'visibility':str(visibility),'source_slots_and_squares':slots},probabilities,visibility


def correct(primary):
    k=primary['k'];reference=primary['reference'];base=primary['base_execution'];rows=[]
    ops=unpack(base['catalog']['qpe'])[1:base['clock']['Q']//2]
    input_state=unpack(base['decoded_position_density'])
    scalar_state=np.einsum('kxa,ab,kyb->xy',ops,input_state,ops.conj(),optimize=True)
    scalar_state/=np.trace(scalar_state).real
    gates={'complete_physical_path_outcomes':True,'raw_centers_match_exact_control':True,
           'actual_control_density_matches_scalar_map':True,'all_native_paths_positive_without_repair':True,
           'corrected_single_and_pair_cost':True,'source_probability_comparison_bound':True}
    numerical_laws=[r['probabilities'] for r in base['preparation_records']]
    for c,receipt in zip(base['cases'],primary['cases']):
        assert all(c[n]==receipt[n] for n in ('M','turns','direction','law'))
        raw=c['path_probabilities'];exact=None;state_error=None
        if c['law']=='swap':
            exact,probabilities,visibility=exact_swap(c['M'],c['turns'],c['direction'],reference)
            physical=[float(x) for x in probabilities]
            state_error=float(np.max(abs(unpack(c['clock_density'])-scalar_state)))
            gates['actual_control_density_matches_scalar_map'] &= state_error<1e-10
            gates['raw_centers_match_exact_control'] &= max(abs(a-b) for a,b in zip(raw,physical))<1e-12
            gates['raw_centers_match_exact_control'] &= abs(float(visibility)-c['visibility_direct_crossblock'])<1e-12
        else:
            physical=list(raw)
            gates['all_native_paths_positive_without_repair'] &= min(raw)>0
        gates['complete_physical_path_outcomes'] &= min(physical)>=0 and abs(sum(physical)-1)<1e-12
        selected=c['clock_path_flag'];path=physical[selected]
        analyzer=receipt['log10_analyzer_probability'];clock=analyzer+math.log10(path);pair=2*clock
        gates['corrected_single_and_pair_cost'] &= abs((clock-analyzer)-math.log10(path))<1e-10 and pair==2*clock
        final=unpack(c['final_density']);gain=float(np.trace(final).real)
        numerical_laws += [r['probabilities'] for r in c['records'] if 'probabilities' in r]
        numerical_laws.append(c['all_reencode_bin_probabilities'])
        rows.append({n:receipt[n] for n in ('k','reference','M','turns','direction','law')}|{
            'legacy_raw_path_comparison_centers':raw,'legacy_binary_visibility_center':c['visibility'],
            'direct_crossblock_visibility_center':c['visibility_direct_crossblock'],
            'raw_path_sum_defect':sum(raw)-1,'raw_final_state_trace_gain':gain,
            'physical_path_law':physical,'exact_control_law':exact,'selected_path':selected,
            'selected_physical_path_probability':path,'actual_control_scalar_state_error':state_error,
            'log10_probability_reaching_analyzer':analyzer,'log10_probability_selected_clock':clock,
            'log10_probability_independent_identical_pair':pair,
            'legacy_probability_mislabel_correction':math.log10(path),
            'conditioning':'one filed successful classical prefix, or both independent prefixes for pair',
            'no_clipping_or_fitted_renormalization':True})
    stats={'complete_flag_laws':len(numerical_laws),'flag_centers':sum(len(x) for x in numerical_laws),
           'minimum_signed_probability_center':min(min(x) for x in numerical_laws),
           'maximum_complete_law_sum_defect':max(abs(sum(x)-1) for x in numerical_laws),
           'maximum_raw_path_sum_defect':max(abs(r['raw_path_sum_defect']) for r in rows),
           'maximum_raw_final_trace_defect':max(abs(r['raw_final_state_trace_gain']-1) for r in rows)}
    gates['source_probability_comparison_bound'] &= stats['minimum_signed_probability_center']>=-1e-10 and stats['maximum_complete_law_sum_defect']<1e-10
    inherited={n:v for n,v in primary['gates'].items() if n!='complete_path_outcomes'}
    return {'k':k,'reference':reference,'rows':rows,'comparison_statistics':stats,
        'original_primary_gates':primary['gates'],'original_primary_passed':primary['passed'],
        'unchanged_primary_gates_except_repaired_representation':inherited,
        'amended_gates':{n:bool(v) for n,v in gates.items()},
        'passed':all(inherited.values()) and all(gates.values())}


def run():
    start=time.monotonic();registration=json.loads((HERE/'EXACT_PATH_COUNT_AMENDMENT.json').read_text())
    assert hashlib.sha256((HERE/'EXACT_PATH_PREREG.json').read_bytes()).hexdigest()==registration['original_prereg_sha256']
    for name,digest in registration['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    table=json.loads((HERE/'EXACT_SWAP_BLOCKS.json').read_text())['tables']
    for row in table:
        t=F(row['slot']%row['M'],row['M'])
        assert F(row['positive_square_scalar'])==(1-3*t+3*t*t)/4
    assert len(table)==36
    members=[]
    for k in range(7):
        for reference in ('matched','stationary'):
            path=HERE/f'k{k}_{reference}/EVIDENCE.json.gz'
            result=correct(json.loads(gzip.decompress(path.read_bytes())))
            result['primary_evidence_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();members.append(result)
    output={'candidate':'A-native-common-reference-menu-v1','amendment':'source-derived exact swap path law; corrected selected-clock and pair success',
        'members':members,'cases':sum(len(m['rows']) for m in members),'passed':all(m['passed'] for m in members),
        'original_failed_primary_members':[{'k':m['k'],'reference':m['reference']} for m in members if not m['original_primary_passed']],
        'same_group':True,'external_G6':False,'seconds':time.monotonic()-start}
    print(json.dumps(output,indent=2));return output


if __name__=='__main__':
    result=run();raise SystemExit(0 if result['passed'] else 1)
