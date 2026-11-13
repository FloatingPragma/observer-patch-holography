from pathlib import Path
import hashlib,json,os,re,subprocess,sys,difflib
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
OUT=Path(os.environ['AUDIT_OUTPUT']);OUT.mkdir(parents=True,exist_ok=True)
manifest=json.loads((HERE/'manifest.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for path,h in manifest['input_sha256'].items():assert sha(ROOT/path)==h,('source bytes changed',path)
report={'phase':sys.argv[1],'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'current_source_files':len(manifest['input_sha256']),'all_source_bytes_exact':True}
if sys.argv[1]!='before':
 fresh=json.loads((ROOT/manifest['mutable_outputs'][0]).read_text());before=json.loads((HERE/'before.json').read_text())
 assert fresh.keys()==before.keys(),('top-level keys changed',fresh.keys()-before.keys())
 sys.path.insert(0,str(ROOT))
 from tools.native_physical_results import collect
 assert fresh['native_spatial_and_field_records']==collect(ROOT),'new source collector differs'
 changes=[];allowed=set(manifest['input_sha256'].values())
 def walk(a,b,path=''):
  assert type(a)is type(b),('type',path)
  if isinstance(a,dict):
   assert a.keys()==b.keys(),('keys',path)
   for k in a:walk(a[k],b[k],path+'/'+k)
  elif isinstance(a,list):
   assert len(a)==len(b),('length',path)
   for i,(x,y)in enumerate(zip(a,b)):walk(x,y,path+'/'+str(i))
  elif a!=b:
   assert isinstance(a,str)and re.fullmatch('[0-9a-f]{64}',a)and re.fullmatch('[0-9a-f]{64}',b)and b in allowed,('old scientific value changed',path,a,b)
   matches=[p for p,h in manifest['input_sha256'].items()if h==b]
   changes.append({'path':path,'before':a,'after':b,'authenticated_files':matches})
 walk({k:v for k,v in before.items()if k!='native_spatial_and_field_records'},{k:v for k,v in fresh.items()if k!='native_spatial_and_field_records'})
 old_native=before['native_spatial_and_field_records']; new_native=fresh['native_spatial_and_field_records']
 assert old_native.keys() <= new_native.keys()
 for k,v in old_native.items():
  if k=='receipts_sha256':
   assert v.keys()<=new_native[k].keys(), 'inherited native receipts removed'
   for receipt_path,h in v.items():walk(h,new_native[k][receipt_path],'/native_spatial_and_field_records/receipts_sha256/'+receipt_path)
  else:walk(v,new_native[k],'/native_spatial_and_field_records/'+k)
 report.update({'old_fields_keys_types_values_exact_except_authenticated_provenance':True,'old_changed_hash_leaves':changes,'new_section_matches_actual_collector':True,'output_sha256':{p:sha(ROOT/p)for p in manifest['mutable_outputs']}})
 (OUT/'postdiction_ledger.json').write_bytes((ROOT/manifest['mutable_outputs'][0]).read_bytes())
 (OUT/'POSTDICTION_LEDGER.md').write_bytes((ROOT/manifest['mutable_outputs'][1]).read_bytes())
 (OUT/'markdown.diff').write_text(''.join(difflib.unified_diff((HERE/'before.md').read_text().splitlines(True),(ROOT/manifest['mutable_outputs'][1]).read_text().splitlines(True),fromfile='committed',tofile='actual-producer')))
(OUT/(sys.argv[1]+'-report.json')).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,sort_keys=True),flush=True)
