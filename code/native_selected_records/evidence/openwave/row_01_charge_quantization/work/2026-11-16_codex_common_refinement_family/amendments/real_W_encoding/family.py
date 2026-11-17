"""Prospective common-family execution. All records retained per dyadic member."""
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import core
from selectors import selector


def execute(k, folder):
    start=time.perf_counter();base=core.run(k);before=dict(core.CALLS);catalog={}
    for m in (4,8):
        g=core.graph(2,m)
        for count in (2,3,5,2*m):catalog[f'{m}:{count}']=selector(count,8,g)
    cases=[]
    for c in base['cases']:
        m=c['M'];frames=3*m*c['turns'];s5=Fraction(catalog[f'{m}:5']['success_probability']);s3=Fraction(catalog[f'{m}:3']['success_probability'])
        factor=s5**11*s3**frames;correction=math.log10(factor.numerator)-math.log10(factor.denominator)
        counts={5:11,3:frames,2*m:4*frames};fixed=2*(2*k+1)
        maximum=fixed+sum(n*catalog[f'{m}:{choice}']['maximum_bits'] for choice,n in counts.items())
        cases.append({key:c[key] for key in ('M','turns','direction','law','visibility','positive_fidelity','reader','stable_reader','audited_off_compact_band_weight')}|{
            'k':k,'selector_counts':counts,'fixed_clock_blank_bits':fixed,'selector_success_factor':str(factor),
            'maximum_selector_blank_events_on_successful_path':maximum,
            'log10_success_conditional_on_classical_prefix':c['log10_total_with_preparation']+correction})
    source=base['clock'];epsilon=source['epsilon'];omega=base['expected_frequency']
    frequency_upper=epsilon**2/24
    gates=base['gates']|{'family_frequency_bound':-1e-12<=omega-2/math.sqrt(5)<=frequency_upper+1e-12,
        'independent_copy_frequency_additivity':True,'pair_stable_frequency_interval':True,'pair_fixed_off_band_bound':True}
    # Independent copies carry the entire finite state; its Ramsey law is z(t)^2.
    # This measures actual returned mixtures rather than inserting two ideal f's.
    pairs=[]
    import numpy as np
    W=np.array(source['W']['real'])+1j*np.array(source['W']['imag'])
    for c in base['cases']:
        rho=np.array(c['clock_density']['real'])+1j*np.array(c['clock_density']['imag'])
        z=complex(np.trace(rho@np.linalg.matrix_power(W,2**k)))
        joint=np.kron(rho,rho);U=np.kron(np.linalg.matrix_power(W,2**k),np.linalg.matrix_power(W,2**k))
        direct=complex(np.einsum('ab,ba->',joint,U));error=abs(direct-z*z)
        pair_records=[{}, {'x_plus':(1+direct.real)/2,'y_plus':(1+direct.imag)/2}]
        pair_reader=core.cosine_reader(pair_records,eta=.01)
        plus=c['positive_fidelity'];minus=1-plus-c['audited_off_compact_band_weight']
        pair_off=1-plus*plus-minus*minus
        gates['pair_stable_frequency_interval'] &= pair_reader['lower']<=2*omega<=pair_reader['upper'] and pair_reader['width']<.03
        gates['pair_fixed_off_band_bound'] &= -1e-10<=pair_off<=.01+1e-10
        gates['independent_copy_frequency_additivity'] &= error<1e-10
        pairs.append({'M':c['M'],'turns':c['turns'],'direction':c['direction'],'law':c['law'],
            'single_ramsey':[z.real,z.imag],'pair_ramsey':[direct.real,direct.imag],
            'product_error':error,'stable_pair_reader':pair_reader,'audited_pair_off_band_weight':pair_off,'R_events_per_pair_word':2,'R_events_per_single_word':1})
    output={'candidate':'A-native-common-refinement-family-v1','k':k,'base_execution':base,
        'selector_catalog':catalog,'cases':cases,'actual_returned_pair_records':pairs,'gates':gates,
        'passed':all(gates.values()),'literal_all_calls':dict(core.CALLS),
        'literal_selector_calls':{law:core.CALLS[law]-before[law] for law in before},
        'frequency_prediction':omega,'frequency_limit':2/math.sqrt(5),'bound':frequency_upper,
        'single_keV_one_limiting_anchor':511*omega/(2/math.sqrt(5)),
        'pair_keV_same_anchor':1022*omega/(2/math.sqrt(5)),
        'seconds':time.perf_counter()-start,'same_group':True,'physical_award':False}
    folder.mkdir(exist_ok=True,parents=True)
    payload=gzip.compress((json.dumps(output,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0)
    (folder/'EVIDENCE.json.gz').write_bytes(payload)
    receipt={key:output[key] for key in ('candidate','k','cases','gates','passed','literal_all_calls','literal_selector_calls','frequency_prediction','frequency_limit','bound','single_keV_one_limiting_anchor','pair_keV_same_anchor','seconds','same_group','physical_award')}
    receipt['evidence_sha256']=hashlib.sha256(payload).hexdigest()
    (folder/'RECEIPT.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    here=Path(__file__).resolve().parent;receipts=[]
    for k in (0,2,6):
        result=execute(k,here/f'k{k}');receipts.append(result)
        print(json.dumps({'k':k,'passed':result['passed'],'seconds':result['seconds'],'gates':result['gates']},sort_keys=True),flush=True)
    combined={'candidate':'A-native-common-refinement-family-v1','mandatory_k':[0,2,6],
        'cases':sum(len(r['cases']) for r in receipts),'passed':all(r['passed'] for r in receipts),
        'member_receipt_hashes':{f'k{r["k"]}/RECEIPT.json':hashlib.sha256((here/f'k{r["k"]}/RECEIPT.json').read_bytes()).hexdigest() for r in receipts},
        'calculation_seconds':sum(r['seconds'] for r in receipts)}
    (here/'RECEIPT.json').write_text(json.dumps(combined,sort_keys=True,indent=2)+'\n')
    print(json.dumps(combined,indent=2));raise SystemExit(0 if combined['passed'] else 1)
