"""Amendment verification of the complete pair source and its saved receipt."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
from validate_descendants import HERE, invoke, snapshot

output=HERE/'validation_pair_receipt';output.mkdir(exist_ok=False);results=[]
with tempfile.TemporaryDirectory(prefix='contract-pair-receipt-') as temp:
    folder=Path(temp);shutil.copytree(HERE/'descendants/producers/pair',folder/'pair')
    script=folder/'pair/pair_source.py';saved=folder/'pair/saved/PAIR_RECEIPT.json'
    results.append(invoke([sys.executable,'-B','-O',str(script)],folder,'complete_saved_verify',output,True))
    original=saved.read_bytes();value=json.loads(original);value['literal_native_calls']+=1;saved.write_text(json.dumps(value))
    before=snapshot(folder/'pair/saved')
    results.append(invoke([sys.executable,'-B',str(script)],folder,'changed_saved_receipt',output,False))
    if before!=snapshot(folder/'pair/saved'):raise RuntimeError('Saved receipt overwritten')
    saved.write_bytes(original);first=folder/'first'
    results.append(invoke([sys.executable,'-B',str(script),'--record-to',str(first)],folder,'explicit_first_record',output,True))
    before=snapshot(first)
    results.append(invoke([sys.executable,'-B',str(script),'--verify',str(first)],folder,'first_record_readonly_verify',output,True))
    results.append(invoke([sys.executable,'-B',str(script),'--record-to',str(first)],folder,'first_record_overwrite_rejected',output,False))
    if before!=snapshot(first):raise RuntimeError('First record changed during verify or refused overwrite')
summary={'passed':all(r['passed'] for r in results),'results':[{k:v for k,v in r.items() if k not in ('stdout','stderr')} for r in results]}
(output/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
raise SystemExit(0 if summary['passed'] else 1)
