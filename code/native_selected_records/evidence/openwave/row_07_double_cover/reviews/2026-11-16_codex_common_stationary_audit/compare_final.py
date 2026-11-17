from pathlib import Path
import json,hashlib,math
from fractions import Fraction as F
root=Path(__file__).resolve().parents[4]; src=root/'openwave/row_01_charge_quantization/work/2026-11-16_codex_common_reference_menu'; dst=root/'openwave/row_07_double_cover/reviews/2026-11-16_codex_common_stationary_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(n,d):(dst/n).write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
manifest=src/'MANIFEST.sha256';assert sha(manifest)=='a9ced490f0bc6fee5ed08aacb6ef9c8613758ee1af243ce3e98aa153928ebf2c'
leaves=0
for line in manifest.read_text().splitlines():
 h,p=line.split('  ',1);assert sha(src/p)==h,p;leaves+=1
primary=json.loads((src/'EXACT_PATH_CORRECTION.json').read_text()); own=json.loads((dst/'SPECTRAL_RESULT.json').read_text()); key=lambda r:tuple(r[n]for n in ['k','reference','M','turns','direction','law']); lookup={key(r):r for r in own['cases']};maxima={k:0. for k in ['visibility','path_probability','log_interference','log_clock','pair_log']};controls=0
for member in primary['members']:
 for r in member['rows']:
  a=lookup[key(r)]; pairs=[('visibility','direct_crossblock_visibility_center'),('path_probability','selected_physical_path_probability'),('log_interference','log10_probability_reaching_analyzer'),('log_clock','log10_probability_selected_clock')]
  for x,y in pairs:maxima[x]=max(maxima[x],abs(a[x]-r[y]));assert abs(a[x]-r[y])<1e-9
  maxima['pair_log']=max(maxima['pair_log'],abs(2*a['log_clock']-r['log10_probability_independent_identical_pair']))
  if r['law']=='swap':
   m,t,d=r['M'],r['turns'],r['direction']; amps=[F(1),F(1)]
   for step in range(1,3*m*t+1):
    q=step%m; rot=(2*m+d*step)%(3*m);ref=(2*m+d*min(q,m-q))%(3*m) if r['reference']=='matched' else 2*m
    for i,slot in enumerate((rot,ref)):
     u=F(slot%m,m);square=(1-3*u+3*u*u)/4;amps[i]*=square**2
   pplus=(amps[0]+amps[1])**2/(2*(amps[0]**2+amps[1]**2));pminus=1-pplus
   law=r['exact_control_law'];assert [F(x)for x in law['amplitudes_without_common_factors']]==amps
   assert [F(x)for x in law['path_probabilities']]==[pplus,pminus]
   assert F(law['visibility'])==pplus-pminus and pplus>pminus>=0
   assert r['selected_physical_path_probability']==float(pplus);controls+=1
assert len(lookup)==224 and controls==112
out={'primary_candidate':primary['candidate'],'primary_manifest_sha256':sha(manifest),'manifest_leaves_checked':leaves,'final_receipt_sha256':sha(src/'FINAL_RECEIPT.json'),'exact_path_correction_sha256':sha(src/'EXACT_PATH_CORRECTION.json'),'own_spectral_result_sha256':sha(dst/'SPECTRAL_RESULT.json'),'cases':224,'rational_control_laws_independently_equal':controls,'maximum_spectral_differences':maxima,'passed':True,'external_G6':False};dump('FINAL_COMPARISON.json',out);print(json.dumps(out,indent=2))
