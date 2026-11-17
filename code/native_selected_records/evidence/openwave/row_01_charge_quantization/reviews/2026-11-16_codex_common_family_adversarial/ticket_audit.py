"""Independent exact clock-ticket and coherent-selector record census."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import gzip,json,time

B=Path(__file__).resolve().parents[2]/'work/2026-11-16_codex_common_refinement_family'
start=time.monotonic();counts={'refinement_words':0,'refinement_blank_records':0,'source_selector_words':0,'retained_rejected_prefixes':0};members=[]
for k in range(7):
    r=json.loads(gzip.decompress((B/f'k{k}/EVIDENCE.json.gz').read_bytes()))
    for emb,g in zip(r['base_execution']['embeddings'],r['base_execution']['graphs']):
        tickets=emb['blank_selector_tickets'];assert len(tickets)==4**k
        words=set();active=0
        for row in tickets:
            bits=row['bits'];assert len(bits)==2*k and all(c in (0,1) for c in bits)
            words.add(tuple(bits));assert len(row['records'])==2*k
            assert row['probability']==1/(4**k)
            assert row['active']==(not any(bits));active+=row['active']
            for bit,record in zip(bits,row['records']):
                assert record==[*g['blank'],bit,0]
                counts['refinement_blank_records']+=1
            counts['refinement_words']+=1
        assert len(words)==4**k and active==1
    for key,s in r['selector_catalog'].items():
        n=s['count'];width=(n-1).bit_length();base=2**width;cap=8;rejects=base-n
        assert s['base']==base and s['width']==width and s['cap']==cap
        assert len(s['attempt_word_catalog'])==base
        for number,row in enumerate(s['attempt_word_catalog']):
            assert row['word']==number and row['bits']==[(number>>(width-1-j))&1 for j in range(width)]
            assert F(row['probability'])==F(1,base)
            assert row['predicate']==('accept' if number<n else 'retry')
            for bit,record in zip(row['bits'],row['native_records']):assert record==[29,30,bit,0]
            assert len(row['native_records'])==width;counts['source_selector_words']+=1
        survival=F(0)
        for row in s['accepted_attempts']:
            t=row['attempt'];prefixes=list(product(range(n,base),repeat=t-1))
            assert [tuple(x) for x in row['prior_rejected_word_prefixes']]==prefixes
            assert row['accepted_words_coherent']==list(range(n))
            mass=F(rejects**(t-1)*n,base**t)
            assert F(row['coarse_outcome_probability'])==mass
            assert F(row['each_joint_word_amplitude_squared'])==F(1,base**t)
            assert row['reference_crossblock_error']<1e-12
            survival+=mass;counts['retained_rejected_prefixes']+=len(prefixes)
        assert survival==1-F(rejects,base)**cap==F(s['success_probability'])
        assert F(s['timeout_probability'])==1-survival
        assert [tuple(x) for x in s['timeout_word_prefixes']]==list(product(range(n,base),repeat=cap))
        assert F(s['conditional_choice_amplitude_squared'])==F(1,n)
        counts['retained_rejected_prefixes']+=len(s['timeout_word_prefixes'])
    members.append({'k':k,'active_probability':str(F(1,4**k)),'folded_choice_dimension':11})
print(json.dumps({'verdict':'PASS','counts':counts,'members':members,'same_group_only':True,
 'producer_imports':False,'seconds':time.monotonic()-start,
 'coherent_map_scope':'Accepted label and rejected prefixes remain coherent/retained as specified by the read source; these record counts alone do not infer coherence from uniform histograms.'},indent=2))
