"""Independent registration, exact-parent and targeted mutation checks."""
from pathlib import Path
import copy,gzip,hashlib,json,subprocess,sys,tempfile,time
HERE=Path(__file__).resolve().parent
P=HERE.parents[1]/'work/2026-11-16_codex_common_reference_menu'
OLD=HERE.parents[1]/'work/2026-11-16_codex_common_refinement_family'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
start=time.monotonic();bindings={}
# The preserved files, rather than the current amended source, bind old freezes.
registrations=[('PREREG.json','source_files'),('AMENDMENT_PREREG.json','source_files'),('EXACT_PATH_PREREG.json','files'),('EXACT_PATH_COUNT_AMENDMENT.json','files')]
archived={p.name:[] for p in P.rglob('*') if p.is_file()}
for p in P.rglob('*'):
 if p.is_file():archived.setdefault(p.name,[]).append(p)
for name,key in registrations:
 registration=json.loads((P/name).read_text());files={}
 for rel,digest in registration[key].items():
  candidates=[P/rel]+archived.get(Path(rel).name,[])
  matches=[str(p.relative_to(P)) for p in candidates if p.exists() and sha(p)==digest]
  assert matches,(name,rel,digest)
  files[rel]={'sha256':digest,'matching_preserved_paths':sorted(set(matches))}
 bindings[name]={'sha256':sha(P/name),'files':files,'utc':registration['utc']}
extra={'visibility_direct_crossblock','visibility_binary_discrepancy','normalized_crossblock','arm_weights','selected_path_probability','log10_reaching_analyzer_before_selectors','log10_selected_clock_before_selectors'}
extra_gates={'native_normalized_cross','binary_crossblock_agreement','literal_stationary_reference'}
for k in range(7):
 fresh=json.loads(gzip.decompress((P/f'k{k}_matched/EVIDENCE.json.gz').read_bytes()))['base_execution']
 prior=json.loads(gzip.decompress((OLD/f'k{k}/EVIDENCE.json.gz').read_bytes()))['base_execution']
 for key in ('candidate','seconds','reference_mode'):fresh.pop(key,None);prior.pop(key,None)
 for c in fresh['cases']:
  for key in extra:c.pop(key,None)
 for key in extra_gates:fresh['gates'].pop(key,None)
 assert fresh==prior,('exact legacy payload preservation',k)
primary=json.loads(gzip.decompress((P/'k0_stationary/EVIDENCE.json.gz').read_bytes()))
correction=json.loads((P/'EXACT_PATH_CORRECTION.json').read_text())
member_index=next(i for i,m in enumerate(correction['members']) if m['k']==0 and m['reference']=='stationary')
def mutate_hold(x,y):next(r for r in x['base_execution']['cases'][0]['records'] if r['kind']=='marker_hold')['native_event']=True
def mutate_cross(x,y):x['base_execution']['cases'][0]['final_density']['real'][0][33]+=.001
def mutate_cost(x,y):y['members'][member_index]['rows'][0]['log10_probability_independent_identical_pair']+=.01
def mutate_physical(x,y):y['members'][member_index]['rows'][1]['physical_path_law'][0]+=.001
def mutate_pair(x,y):x['cases'][0]['pair_records']['8']['x_plus']+=.001
mutations=[]
for name,mutate in [('stationary_native_event',mutate_hold),('reference_crossblock',mutate_cross),('pair_cost_factor',mutate_cost),('physical_control_law',mutate_physical),('stride8_record',mutate_pair)]:
 x=copy.deepcopy(primary);y=copy.deepcopy(correction);mutate(x,y)
 with tempfile.TemporaryDirectory(prefix='oph-reference-review-') as tmp:
  path=Path(tmp);(path/'k0_stationary').mkdir()
  (path/'k0_stationary/EVIDENCE.json.gz').write_bytes(gzip.compress(json.dumps(x).encode()))
  (path/'EXACT_PATH_CORRECTION.json').write_text(json.dumps(y))
  result=subprocess.run([sys.executable,str(HERE/'audit.py'),str(path),'0','stationary'],capture_output=True,text=True)
  assert result.returncode!=0,name
  mutations.append({'name':name,'rejected':True,'returncode':result.returncode,'error':result.stderr.splitlines()[-1]})
result={'passed':True,'registrations':bindings,'exact_legacy_cases_preserved':112,'mutations':mutations,'seconds':time.monotonic()-start,'same_group_only':True,'producer_imports':False}
print(json.dumps(result,indent=2))
