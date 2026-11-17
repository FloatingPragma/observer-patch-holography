"""Common reference-menu successor; fresh complete primary records, fixed gates."""
from pathlib import Path
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
import time
import numpy as np
import core
from native_selectors import selector
from reader import read

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'2026-11-16_codex_common_refinement_family'
CANDIDATE='A-native-common-reference-menu-v1'
SCALE=511*math.sqrt(5)/2
OLD_PAIR=1024.4508487023415
EXTRA_CASE={'visibility_direct_crossblock','visibility_binary_discrepancy','normalized_crossblock','arm_weights',
            'selected_path_probability','log10_reaching_analyzer_before_selectors','log10_selected_clock_before_selectors'}
EXTRA_GATES={'native_normalized_cross','binary_crossblock_agreement','literal_stationary_reference'}


def unpack(x):
    return np.array(x['real'])+1j*np.array(x['imag']) if isinstance(x,dict) else np.array(x)


def probabilities(z):
    return {'x_plus':float((1+z.real)/2),'x_minus':float((1-z.real)/2),
            'y_plus':float((1+z.imag)/2),'y_minus':float((1-z.imag)/2)}


def exact_swap_blocks():
    """Exact rational optical square; no state or phase is installed."""
    def multiply(a,b):return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    optical=[[F(int(i==j))-F(1,3) for j in range(3)] for i in range(3)]
    tables=[]
    for m in (4,8):
        for slot in range(3*m):
            edges=core.routing(slot,m);P=[[F(0) for _ in range(3)] for _ in range(3)]
            for a,b in edges:
                for j in range(3):
                    i=b if j==a else a if j==b else j
                    P[i][j]+=F(1,len(edges))
            Q=multiply(optical,multiply(P,optical));Q2=multiply(Q,Q)
            scalar=sum(Q2[i][i] for i in range(3))/2
            assert scalar>0 and Q2==[[scalar*x for x in row] for row in optical]
            tables.append({'M':m,'slot':slot,'positive_square_scalar':str(scalar),
                           'four_target_events_scalar':str(scalar*scalar/F(9))})
    return {'tables':tables,'passed':True,'meaning':'Every accepted four-target-event swap block is a positive scalar on the whole optical plane. The shared plane filter and reencoding preserve the positive crossblock for every admitted optical density; no floating binary resolution is inferred.'}


def compare_matched(base,old):
    """Every legacy output, including full states/flags/current histories, unchanged."""
    left={k:v for k,v in base.items() if k not in ('candidate','reference_mode','seconds')}
    left['cases']=[{k:v for k,v in c.items() if k not in EXTRA_CASE} for c in base['cases']]
    left['gates']={k:v for k,v in base['gates'].items() if k not in EXTRA_GATES}
    right={k:v for k,v in old.items() if k not in ('candidate','seconds')}
    assert left==right,'matched mode differs from sealed parent beyond explicit new fields/runtime'
    return True


def execute(k,reference,folder):
    start=time.perf_counter();base=core.run(k,reference);old=json.loads(gzip.decompress((OLD/f'k{k}/EVIDENCE.json.gz').read_bytes()))
    equivalence=compare_matched(base,old['base_execution']) if reference=='matched' else None
    before=dict(core.CALLS);catalog={}
    for m in (4,8):
        for count in (2,3,5,2*m):catalog[f'{m}:{count}']=selector(count,8,core.graph(2,m))
    W=unpack(base['clock']['W']);unit=np.linalg.matrix_power(W,2**k)
    joint_words={n:np.kron(np.linalg.matrix_power(unit,n),np.linalg.matrix_power(unit,n)) for n in (1,8)}
    cases=[];cost_corrections=[]
    gates=base['gates']|{'matched_parent_equivalence':equivalence is not False,'pair_impurity':True,
        'pair_unique_branch':True,'pair_interval_contains_prediction':True,'k6_pair_excludes_old_value':True,
        'selected_path_cost_included':True,'corrected_pair_cost':True,'complete_path_outcomes':True}
    old_index={(c['M'],c['turns'],c['direction'],c['law']):c for c in old['cases']}
    for c in base['cases']:
        m=c['M'];frames=3*m*c['turns'];s3=F(catalog[f'{m}:3']['success_probability']);s5=F(catalog[f'{m}:5']['success_probability'])
        factor=s5**11*s3**frames;adjustment=math.log10(factor.numerator)-math.log10(factor.denominator)
        analyzer=c['log10_reaching_analyzer_before_selectors']+adjustment
        clock=analyzer+math.log10(c['selected_path_probability']);pair=2*clock
        rho=unpack(c['clock_density']);joint=np.kron(rho,rho)
        moments={n:complex(np.einsum('ab,ba->',joint,U)) for n,U in joint_words.items()}
        records={n:probabilities(z) for n,z in moments.items()}
        positive=c['positive_fidelity'];negative=1-positive-c['audited_off_compact_band_weight']
        off=1-positive**2-negative**2
        interval=read(records[1],records[8],eta=.01,probability_error=2e-5)
        parts=interval['frequency_components'];expected=2*base['expected_frequency']
        gates['pair_impurity'] &= -1e-10<=off<=.01
        gates['pair_unique_branch'] &= interval['unique_component']
        gates['pair_interval_contains_prediction'] &= any(a<=expected<=b for a,b in parts)
        if k==6:gates['k6_pair_excludes_old_value'] &= interval['unique_component'] and SCALE*parts[0][1]<OLD_PAIR
        gates['selected_path_cost_included'] &= abs(clock-c['log10_selected_clock_before_selectors']-adjustment)<1e-10
        gates['corrected_pair_cost'] &= abs(pair-2*(analyzer+math.log10(c['selected_path_probability'])))<1e-10
        gates['complete_path_outcomes'] &= min(c['path_probabilities'])>=0 and abs(sum(c['path_probabilities'])-1)<1e-12
        tag={key:c[key] for key in ('M','turns','direction','law')}
        result=tag|{'k':k,'reference':reference,'visibility_binary':c['visibility'],'visibility_direct':c['visibility_direct_crossblock'],
             'normalized_crossblock':c['normalized_crossblock'],'arm_weights':c['arm_weights'],
             'binary_crossblock_discrepancy':c['visibility_binary_discrepancy'],
             'path_probabilities':c['path_probabilities'],'selected_path':c['clock_path_flag'],
             'selector_success_factor':str(factor),'log10_analyzer_probability':analyzer,
             'log10_selected_clock_probability':clock,'log10_independent_pair_preparation_probability':pair,
             'conditioning':'filed successful classical prefix for one, both independent filed prefixes for pair',
             'pair_records':records,'pair_reader':interval,'pair_off_band':off,
             'pair_frequency_source_comparison':expected,'pair_energy_interval_keV':[[SCALE*a,SCALE*b] for a,b in parts],
             'settings_cost':'Each X/Y and stride is a separate measurement setting; pair preparation probability is per prepared pair, not the success of collecting all settings.',
             'R_events_per_first_binary_setting':2*2**k,'R_events_per_stride8_binary_setting':16*2**k}
        if reference=='matched':
            prior=old_index[(c['M'],c['turns'],c['direction'],c['law'])]['log10_success_conditional_on_classical_prefix']
            assert prior==analyzer
            cost_corrections.append(tag|{'k':k,'old_mislabelled_success':prior,'probability_scope_of_old_field':'reaching analyzer before selecting returned clock',
                'selected_path_probability':c['selected_path_probability'],'log10_selected_clock':clock,'log10_pair':pair,
                'amendment':math.log10(c['selected_path_probability'])})
        cases.append(result)
    output={'candidate':CANDIDATE,'k':k,'reference':reference,'base_execution':base,'selector_catalog':catalog,
        'cases':cases,'historical_cost_corrections':cost_corrections,'gates':{x:bool(y) for x,y in gates.items()},
        'passed':all(gates.values()),'literal_all_calls':dict(core.CALLS),
        'literal_selector_calls':{law:core.CALLS[law]-before[law] for law in before},
        'seconds':time.perf_counter()-start,'same_group':True,'external_G6':False,'physical_award':False}
    folder.mkdir(parents=True,exist_ok=True)
    payload=gzip.compress((json.dumps(output,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0)
    (folder/'EVIDENCE.json.gz').write_bytes(payload)
    receipt={key:output[key] for key in ('candidate','k','reference','cases','historical_cost_corrections','gates','passed','literal_all_calls','literal_selector_calls','seconds','same_group','external_G6','physical_award')}
    receipt['evidence_sha256']=hashlib.sha256(payload).hexdigest()
    (folder/'RECEIPT.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    return receipt


def run():
    prereg=json.loads((HERE/'PREREG.json').read_text())
    for name,digest in prereg['source_files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    inputs=json.loads((HERE/'INPUTS.json').read_text())
    for name,digest in inputs['files'].items():assert hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest()==digest,name
    exact=exact_swap_blocks();(HERE/'EXACT_SWAP_BLOCKS.json').write_text(json.dumps(exact,indent=2)+'\n')
    receipts=[]
    for k in range(7):
        for reference in ('matched','stationary'):
            result=execute(k,reference,HERE/f'k{k}_{reference}');receipts.append(result)
            print(json.dumps({'k':k,'reference':reference,'passed':result['passed'],'seconds':result['seconds'],'gates':result['gates']}),flush=True)
    summary={'candidate':CANDIDATE,'cases':sum(len(r['cases']) for r in receipts),'matched_cases':112,'stationary_cases':112,
        'passed':all(r['passed'] for r in receipts),'historical_cost_corrections':sum(len(r['historical_cost_corrections']) for r in receipts),
        'member_receipt_hashes':{f'k{r["k"]}_{r["reference"]}/RECEIPT.json':hashlib.sha256((HERE/f'k{r["k"]}_{r["reference"]}/RECEIPT.json').read_bytes()).hexdigest() for r in receipts},
        'calculation_seconds':sum(r['seconds'] for r in receipts),'exact_swap_block_count':len(exact['tables']),
        'same_group':True,'external_G6':False,'physical_award':False}
    (HERE/'RECEIPT.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2));return summary


if __name__=='__main__':
    result=run();raise SystemExit(0 if result['passed'] else 1)
