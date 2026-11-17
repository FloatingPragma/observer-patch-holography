"""Integrity/sensitivity checks, separate from scientific alternative laws."""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import time
from audit import run

here=Path(__file__).resolve().parent
original=json.loads(gzip.decompress((here/'first_run/EVIDENCE.json.gz').read_bytes()))

def break_crossblock(d):
    matrix=d['cases'][0]['final_density']
    for part in ('real','imag'):
        for a in range(33):
            for b in range(33,66):matrix[part][a][b]=matrix[part][b][a]=0.

def wrong_index(d):
    matrix=d['cases'][0]['reencode_input']
    matrix['real'][0],matrix['real'][1]=matrix['real'][1],matrix['real'][0]

def false_current(d):d['charge_records'][2]['body'][4][3]+=1

def nonedge(d):d['cases'][0]['records'][0]['physical_edge'][0]=0

def false_word(d):d['transport']['traces'][0][0][0]=17

results=[];begin=time.perf_counter()
for name,mutate in [('lost_reference_crossblocks',break_crossblock),('wrong_position_index',wrong_index),
                    ('false_current',false_current),('marker_nonedge',nonedge),('wrong_work_word',false_word)]:
    with tempfile.TemporaryDirectory(prefix='common-source-mutation-') as tmp:
        folder=Path(tmp);(folder/'first_run').mkdir()
        for filename in ('common.py','DESIGN.md','TASK.md','LEDGER.md','PREREG.json'):
            shutil.copy2(here/filename,folder/filename)
        d=deepcopy(original);mutate(d)
        blob=gzip.compress(json.dumps(d,separators=(',',':')).encode(),mtime=0)
        (folder/'first_run/EVIDENCE.json.gz').write_bytes(blob)
        receipt=json.loads((here/'first_run/RECEIPT.json').read_text());receipt['evidence_sha256']=hashlib.sha256(blob).hexdigest()
        (folder/'first_run/RECEIPT.json').write_text(json.dumps(receipt))
        rejected=False
        try:run(folder)
        except AssertionError:rejected=True
        results.append({'mutation':name,'custody_recomputed_so_rejection_is_semantic':True,'rejected':rejected})
        assert rejected,name
print(json.dumps({'type':'integrity and detector tests, not alternative-law evidence','results':results,
                  'all_rejected':True,'seconds':time.perf_counter()-begin},indent=2))
