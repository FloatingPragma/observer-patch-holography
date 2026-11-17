"""Native coherent finite selectors and a charged common-source successor."""
from fractions import Fraction
import gzip
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import time
import numpy as np
import common


def selector(count,cap,graph):
    width=(count-1).bit_length();base=2**width;good=tuple(range(count));bad=tuple(range(count,base))
    words=[]
    for word in range(base):
        bits=[(word>>(width-1-j))&1 for j in range(width)];records=[]
        for bit in bits:
            out,record=common.native([0]*graph['N'],tuple(graph['blank']),bit)
            assert not any(out) and record[3]==0;records.append(record)
        words.append({'word':word,'bits':bits,'native_records':records,'probability':str(Fraction(1,base)),
                      'predicate':'accept' if word<count else 'retry'})
    # Attempt index is a coarse flag; all rejected words and the accepted word
    # stay as coherent memory, never replaced by a classical uniform mixture.
    bell=np.array([1.,0.,0.,1.],complex)/math.sqrt(2)
    target=np.outer(bell,bell.conj());attempts=[];success=Fraction(0);max_error=0.
    for t in range(1,cap+1):
        prefixes=list(product(bad,repeat=t-1))
        if not prefixes:break
        p=Fraction(len(prefixes)*count,base**t);success+=p
        amps=np.ones((len(prefixes)*count,1),complex)/math.sqrt(base**t)
        joint=amps*bell[None,:]
        retained=joint.T@joint.conj()
        error=float(np.max(abs(retained-float(p)*target)));max_error=max(max_error,error)
        assert error<1e-12
        attempts.append({'attempt':t,'prior_rejected_word_prefixes':prefixes,
            'accepted_words_coherent':good,'each_joint_word_amplitude_squared':str(Fraction(1,base**t)),
            'coarse_outcome_probability':str(p),'reference_crossblock_error':error})
    timeouts=list(product(bad,repeat=cap)) if bad else []
    timeout=Fraction(len(timeouts),base**cap)
    assert success+timeout==1
    expected_bits=sum(Fraction(width)*Fraction(len(bad),base)**n for n in range(cap))
    conditional_bits=sum(Fraction(a['coarse_outcome_probability'])*a['attempt']*width for a in attempts)/success
    return {'count':count,'width':width,'base':base,'cap':cap,'attempt_word_catalog':words,
        'accepted_attempts':attempts,'timeout_word_prefixes':timeouts,'timeout_probability':str(timeout),
        'success_probability':str(success),'expected_bits_per_entered_selector':str(expected_bits),
        'mean_bits_conditioned_on_success':str(conditional_bits),'maximum_bits':width*(cap if bad else 1),
        'coherent_reference_max_error':max_error,
        'conditional_choice_amplitude_squared':str(Fraction(1,count))}


def run():
    start=time.perf_counter();base=common.run();pre_counts=dict(common.CALLS)
    catalog={};tables=[]
    for m in (4,8):
        graph=common.graph(2,m)
        for count in (2,3,5,16,2*m):
            key=f'{m}:{count}'
            if key not in catalog:catalog[key]=selector(count,8,graph)
        # Every accepted target word supplies the actual ordered star edge in
        # BOTH arms. The word survives the event; address workspace uncomputes.
        for slot in range(3*m):
            v,k=divmod(slot,m);table=[]
            for word in range(2*m):
                ticket,edgebit=divmod(word,2);star=v if ticket<m-k else (v+1)%3
                local=[e for e in common.TRI if star in e][edgebit]
                assert local==common.routing(slot,m)[word]
                physical=[graph['target'][1],graph['target'][2],0]
                table.append({'word':word,'ticket':ticket,'edgebit':edgebit,'local_edge':local,
                              'physical_edge':[physical[x] for x in local]})
            tables.append({'M':m,'slot':slot,'coherent_word_routing':table})
    cases=[]
    for case in base['cases']:
        m=case['M'];frames=3*m*case['turns']
        s5=Fraction(catalog[f'{m}:5']['success_probability']);s3=Fraction(catalog[f'{m}:3']['success_probability'])
        # Nine weighted-triangle cooling preparations and two clock-A encoders.
        counts={5:11,3:frames,2*m:4*frames}
        factor=s5**11*s3**frames
        fixed_clock_bits=2*(4+1)  # four epsilon bits + one source-coin bit per A
        max_bits=fixed_clock_bits+sum(n*catalog[f'{m}:{choice}']['maximum_bits'] for choice,n in counts.items())
        mean_bits=Fraction(fixed_clock_bits)+sum(n*Fraction(catalog[f'{m}:{choice}']['mean_bits_conditioned_on_success']) for choice,n in counts.items())
        correction=math.log10(factor.numerator)-math.log10(factor.denominator)
        # Selector registers are independent of path/reference; its scalar
        # success multiplies the COMPLETE conditional output, not just traces.
        rho=np.array(case['final_density']['real'])+1j*np.array(case['final_density']['imag'])
        charged=float(factor)*rho
        assert np.max(abs(charged/np.trace(charged)-rho))<1e-12
        cases.append({k:case[k] for k in ('M','turns','direction','law','visibility','positive_fidelity','reader')}|{
            'selector_counts':counts,'clock_fixed_blank_bits':fixed_clock_bits,'selector_success_factor':str(factor),
            'log10_selector_factor':correction,'log10_success_conditional_on_classical_prefix':case['log10_total_with_preparation']+correction,
            'maximum_selector_blank_events_on_successful_path':max_bits,
            'mean_selector_blank_events_conditioned_on_all_selectors_success':str(mean_bits),
            'unchanged_full_conditional_density_error':float(np.max(abs(charged/np.trace(charged)-rho)))})
    added={k:common.CALLS[k]-pre_counts[k] for k in common.CALLS}
    return {'candidate':'A-native-selected-source-menu-v2','base_execution':base,'native_selector_catalog':catalog,
        'target_word_tables':tables,'cases':cases,'passed':base['passed'],
        'literal_selector_enumeration_calls':added,'literal_all_calls':dict(common.CALLS),
        'complete_attempt_register_kept':True,'reflection_and_fourier_memory_services_still_declared':True,
        'seconds':time.perf_counter()-start,'same_group':True,'physical_award':False}


if __name__=='__main__':
    output=run();folder=Path(__file__).resolve().parent/'v2_run';folder.mkdir(exist_ok=True)
    blob=gzip.compress((json.dumps(output,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0)
    (folder/'EVIDENCE.json.gz').write_bytes(blob)
    receipt={k:output[k] for k in ('candidate','passed','cases','literal_selector_enumeration_calls','literal_all_calls','seconds','same_group','physical_award')}
    receipt['evidence_sha256']=hashlib.sha256(blob).hexdigest()
    (folder/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,indent=2,sort_keys=True));assert output['passed']
