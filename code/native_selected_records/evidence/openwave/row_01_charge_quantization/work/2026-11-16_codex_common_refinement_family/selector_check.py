"""Independent family native-word, prefix and conditional-cost verification."""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time


def run(folder):
    begin=time.perf_counter();words=cases=prefixes=0
    for k in range(7):
        d=json.loads(gzip.decompress((folder/f'k{k}/EVIDENCE.json.gz').read_bytes()))
        for e in d['base_execution']['embeddings']:
            assert len(e['blank_selector_tickets'])==4**k
            for word in e['blank_selector_tickets']:
                assert len(word['bits'])==len(word['records'])==2*k
                assert word['active']==(not any(word['bits'])) and word['probability']==1/(4**k)
                assert all(c==bit and j==0 for (a,b,c,j),bit in zip(word['records'],word['bits']));words+=1
        for name,s in d['selector_catalog'].items():
            n=s['count'];N=2**s['width'];total=Fraction(0)
            for w in s['attempt_word_catalog']:
                assert sum(b*2**(s['width']-1-i) for i,b in enumerate(w['bits']))==w['word']
                assert all(c==bit and j==0 for (a,b,c,j),bit in zip(w['native_records'],w['bits']));words+=1
            for a in s['accepted_attempts']:
                t=a['attempt'];history=a['prior_rejected_word_prefixes'];prefixes+=len(history)
                assert len(history)==(N-n)**(t-1) and len(set(map(tuple,history)))==len(history)
                assert all(len(h)==t-1 and all(n<=w<N for w in h) for h in history)
                assert a['accepted_words_coherent']==list(range(n))
                probability=Fraction(n*(N-n)**(t-1),N**t)
                assert Fraction(a['coarse_outcome_probability'])==probability;total+=probability
            assert total==Fraction(s['success_probability']) and total+Fraction(s['timeout_probability'])==1
        for c,original in zip(d['cases'],d['base_execution']['cases']):
            m=c['M'];frames=3*m*c['turns'];cat=d['selector_catalog']
            factor=Fraction(cat[f'{m}:5']['success_probability'])**11*Fraction(cat[f'{m}:3']['success_probability'])**frames
            assert factor==Fraction(c['selector_success_factor'])
            assert c['fixed_clock_blank_bits']==2*(2*k+1)
            assert c['selector_counts']=={'5':11,'3':frames,str(2*m):4*frames};cases+=1
    return {'passed':True,'same_group':True,'external_G6':False,'native_word_records':words,
            'coherent_rejected_prefixes':prefixes,'charged_cases':cases,'seconds':time.perf_counter()-begin}


if __name__=='__main__':print(json.dumps(run(Path(__file__).resolve().parent),indent=2,sort_keys=True))
