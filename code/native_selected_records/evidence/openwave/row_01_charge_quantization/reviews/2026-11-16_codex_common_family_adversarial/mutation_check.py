"""Mutate temporary record copies; leave all producer evidence untouched."""
from pathlib import Path
import copy,gzip,json,subprocess,sys,tempfile,time

HERE=Path(__file__).resolve().parent
B=HERE.parents[1]/'work/2026-11-16_codex_common_refinement_family'
original=json.loads(gzip.decompress((B/'k0/EVIDENCE.json.gz').read_bytes()))
mutations=[]
with tempfile.TemporaryDirectory() as td:
    target=Path(td);(target/'k0').mkdir()
    for label in ('native-current','preparation-coherence','outer-reference-crossblock','counter-complement','pair-off-band'):
        r=copy.deepcopy(original);b=r['base_execution']
        if label=='native-current':b['charge_records'][0]['body'][0][3]+=1
        elif label=='preparation-coherence':b['preparation_records'][1]['input_density'][0][1]+=.01
        elif label=='outer-reference-crossblock':
            row=next(x for x in b['cases'][0]['records'] if x['kind']=='target')
            row['input_density']['real'][0][4]+=.01
        elif label=='counter-complement':b['cases'][0]['all_reencode_bin_probabilities'][0]+=.001
        else:r['actual_returned_pair_records'][0]['audited_pair_off_band_weight']=0.
        (target/'k0/EVIDENCE.json.gz').write_bytes(gzip.compress(json.dumps(r).encode()))
        start=time.monotonic();run=subprocess.run([sys.executable,'-B',str(HERE/'audit_optimized.py'),str(target),'0'],capture_output=True,text=True,timeout=60)
        assert run.returncode!=0,label
        mutations.append({'mutation':label,'rejected':True,'exit_code':run.returncode,
                          'reason':run.stderr.strip().splitlines()[-1],'seconds':time.monotonic()-start})
print(json.dumps({'passed':True,'temporary_only':True,'mutations':mutations},indent=2))
