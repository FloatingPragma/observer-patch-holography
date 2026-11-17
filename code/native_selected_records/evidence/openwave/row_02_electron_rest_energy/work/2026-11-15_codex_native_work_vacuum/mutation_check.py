"""Read-only to the scientific packet; mutations are isolated in temporary copies."""
from pathlib import Path
import copy,json,subprocess,sys,tempfile,time
BASE=Path(__file__).resolve().parent
old=json.loads((BASE/'OUTPUT.json').read_text())
new=json.loads((BASE/'refinement/OUTPUT.json').read_text())
mutations=[('invalid_native_coin','old',lambda d:d['source_records'][0].__setitem__(3,2)),
 ('first_work_current','new',lambda d:d['work_traces'][0][0].__setitem__(3,7)),
 ('selector_coin','new',lambda d:d['selector_full_records'][1]['tickets'][0]['bits'].__setitem__(0,2)),
 ('unique_fourier_bin','new',lambda d:d['refinements'][6]['all_counter_probabilities'].__setitem__(117,0.2)),
 ('negative_frequency_erased','new',lambda d:d['refinements'][4]['record_only_reader']['poles'][0].__setitem__('phase',1.0)),
 ('wrong_clock_count','new',lambda d:d['refinements'][6]['ramsey_records'][4].__setitem__('native_W_words',1)),
 ('inflated_native_count','new',lambda d:d.__setitem__('literal_native_function_calls',999999))]
results=[]
for name,which,mutate in [('clean','old',lambda d:None)]+mutations:
    a,b=copy.deepcopy(old),copy.deepcopy(new);mutate(a if which=='old' else b)
    with tempfile.TemporaryDirectory(prefix='row02-review-') as td:
        p=Path(td);(p/'refinement').mkdir();(p/'OUTPUT.json').write_text(json.dumps(a));(p/'refinement/OUTPUT.json').write_text(json.dumps(b))
        start=time.perf_counter();r=subprocess.run([sys.executable,'-B',str(BASE/'check.py'),td],capture_output=True,text=True,timeout=60)
        results.append({'name':name,'exit_code':r.returncode,'seconds':time.perf_counter()-start,'last_error_line':r.stderr.splitlines()[-1] if r.stderr else None})
assert results[0]['exit_code']==0 and all(r['exit_code']!=0 for r in results[1:])
print(json.dumps({'scope':'Integrity/detector mutation checks, not emergence evidence','results':results},indent=2))
