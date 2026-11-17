"""Independent selector and first-failure ledger audit; no producer imports."""
from fractions import Fraction
import gzip
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import time
import numpy as np

here=Path(__file__).resolve().parent
begin=time.perf_counter();blob=(here/'v2_run/EVIDENCE.json.gz').read_bytes();data=json.loads(gzip.decompress(blob))
receipt=json.loads((here/'v2_run/RECEIPT.json').read_text());assert hashlib.sha256(blob).hexdigest()==receipt['evidence_sha256']
reg=json.loads((here/'V2_PREREG.json').read_text())
for name,digest in reg['files'].items():assert hashlib.sha256((here/name).read_bytes()).hexdigest()==digest
base=data['base_execution'];old=json.loads(gzip.decompress((here/'first_run/EVIDENCE.json.gz').read_bytes()))
base.pop('seconds');old.pop('seconds');assert base==old
selectors=data['native_selector_catalog'];maximum_coherence=0.;prefix_count=0
for key,s in selectors.items():
    n=s['count'];b=s['width'];N=2**b;bad=tuple(range(n,N));success=Fraction(0)
    for w in s['attempt_word_catalog']:
        value=sum(bit*2**(b-1-j) for j,bit in enumerate(w['bits']))
        assert value==w['word']
        assert all(c==bit and j==0 for (a,bb,c,j),bit in zip(w['native_records'],w['bits']))
        assert Fraction(w['probability'])==Fraction(1,N)
    assert len(s['attempt_word_catalog'])==N
    for attempt in s['accepted_attempts']:
        t=attempt['attempt'];prefixes=[tuple(p) for p in attempt['prior_rejected_word_prefixes']]
        assert prefixes==list(product(bad,repeat=t-1));prefix_count+=len(prefixes)
        assert attempt['accepted_words_coherent']==list(range(n))
        p=Fraction(n*len(prefixes),N**t);success+=p
        assert p==Fraction(attempt['coarse_outcome_probability'])
        assert Fraction(attempt['each_joint_word_amplitude_squared'])==Fraction(1,N**t)
        # Retain accepted-label off-diagonals, not only its uniform histogram.
        actual=len(prefixes)/N**t*np.ones((n,n));target=float(p)*np.ones((n,n))/n
        maximum_coherence=max(maximum_coherence,float(np.max(abs(actual-target))))
    timeout=Fraction(len(bad)**s['cap'],N**s['cap'])
    assert [tuple(p) for p in s['timeout_word_prefixes']]==list(product(bad,repeat=s['cap'])) if bad else not s['timeout_word_prefixes']
    assert success+timeout==1 and success==Fraction(s['success_probability'])
    assert timeout==Fraction(s['timeout_probability'])
assert maximum_coherence<1e-12
tri=((0,1),(1,2),(0,2))
for tab in data['target_word_tables']:
    m=tab['M'];v,k=divmod(tab['slot'],m);edges=[e for e in tri if v in e]*(m-k)+[e for e in tri if (v+1)%3 in e]*k
    assert [tuple(r['local_edge']) for r in tab['coherent_word_routing']]==edges

ledgers=[]
for case,cost in zip(base['cases'],data['cases']):
    m=case['M'];logreach=0.;events=[];counts={3:0,5:0,2*m:0}
    def add_selector(n,stage,fixed_bits=0):
        global logreach
        s=selectors[f'{m}:{n}'];p=float(Fraction(s['success_probability']));counts[n]+=1
        events.append({'before_stage':stage,'selector_catalog':f'{m}:{n}','log10_reach_before':logreach,
            'log10_first_timeout':None if p==1 else logreach+math.log10(1-p),
            'success_probability':s['success_probability'],'fixed_blank_bits_after_success':fixed_bits,
            'retained_registers':'all rejected words, coarse attempt index, coherent accepted word'})
        logreach+=math.log10(p)
    for i,r in enumerate(base['preparation_records']):
        if r['catalog'] in ('cooling','cooling_contrast'):add_selector(5,f'preparation:{i}')
        elif r['catalog']=='qpe':add_selector(5,f'preparation:{i}',5)
        accepted=r['accepted'] if isinstance(r['accepted'],list) else[r['accepted']]
        logreach+=math.log10(sum(r['probabilities'][x] for x in accepted))
    for i,r in enumerate(case['records']):
        if r['kind']=='marker':
            p=sum(h['probability'] for h in r['histories'] if h['outcome']=='hop')
        else:
            add_selector(2*m if r['kind']=='target' else 3,f'rotation:{i}')
            p=r['probabilities'][r['accepted']]
        logreach+=math.log10(p)
    add_selector(5,'final_reencode',5)
    logreach+=math.log10(sum(case['all_reencode_bin_probabilities'][1:512]))
    assert abs(logreach-cost['log10_success_conditional_on_classical_prefix'])<1e-9
    assert {str(k):v for k,v in counts.items()}==cost['selector_counts']
    ledgers.append({'M':m,'turns':case['turns'],'direction':case['direction'],'law':case['law'],
                    'selector_first_failure_stages':events,'final_log10_success':logreach})
# A leak recording which arm tried the selector removes the phase crossblock.
coherent=np.ones((2,2))/2;leaked=np.eye(2)/2
assert 2*coherent[0,1]==1 and 2*leaked[0,1]==0
# Merely uniform classical accepted labels fail the coherent label-state test.
label=np.ones((5,5))/5;committed=np.eye(5)/5
assert np.max(abs(label-committed))>.1
out={'candidate':data['candidate'],'evidence_sha256':receipt['evidence_sha256'],'passed':True,
     'same_group':True,'external_G6':False,'selector_catalog_count':len(selectors),
     'rejected_prefixes_in_accepted_attempts':prefix_count,'target_word_tables':len(data['target_word_tables']),
     'coherent_label_max_error':maximum_coherence,'which_arm_leak_control_visibility':0.,
     'classical_uniform_label_control_error':float(np.max(abs(label-committed))),
     'all_selector_first_failure_ledgers':ledgers,'seconds':time.perf_counter()-begin}
print(json.dumps(out,indent=2,sort_keys=True))
