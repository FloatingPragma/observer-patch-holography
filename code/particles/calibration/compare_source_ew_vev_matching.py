"""Descriptive comparison of every prescribed electroweak matching control.

Targets are consumed only here, after forward matching. The primary scale,
orders and top method are fixed before comparison; none is chosen by agreement.
No joint significance or model-uncertainty estimate is constructed.
"""
from __future__ import annotations
import argparse,copy,hashlib,json
from decimal import Decimal,localcontext
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PACKET=HERE.parent/'runs/calibration/source_ew_vev_matching.json'
TARGETS=HERE/'source_ew_comparison_targets.json'
MATCHING_SPEC=HERE/'source_ew_vev_matching_spec.json'
PRIOR_TOP=HERE.parent/'data/pdg_2025_conditional_quark_comparison.json'
OUTPUT=HERE.parent/'runs/calibration/source_ew_vev_matching_comparison.json'
REVIEWED_TARGETS_SHA256='4c740d14565123a8fc8534c004309ccc986ae86efb5523917339c4dc7152cc83'
sys.path.insert(0,str(HERE))

def require(ok,message):
 if not ok:raise ValueError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def serialized(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def strict_load(path):
 def pairs(items):
  result={}
  for k,v in items:
   require(k not in result,'duplicate JSON key');result[k]=v
  return result
 def invalid(value):raise ValueError('nonfinite JSON value')
 return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=invalid)
def read_targets():
 require(sha(TARGETS)==REVIEWED_TARGETS_SHA256,'unreviewed target fixture')
 fixture=strict_load(TARGETS)
 require(sha(PRIOR_TOP)==fixture['inherited_top_fixture']['sha256'],'inherited top fixture changed')
 top=next(r for r in strict_load(PRIOR_TOP)['masses'] if r['quark']=='t')
 row=next(r for r in fixture['targets'] if r['quantity']=='top')
 require(row['central']==top['central'] and row['quoted_uncertainty']==top['plus']==top['minus'],'top source row mismatch')
 return fixture

def comparison_value(value,target):
 require(type(value) in (int,float),'comparison value must be numeric')
 value=Decimal(str(value));central=Decimal(target['central'])
 require(value.is_finite() and value>0 and central.is_finite() and central>0,'positive finite comparison values')
 difference=value-central
 return {'computed':str(value),'reference':copy.deepcopy(target),'unit':target['unit'],
         'signed_difference':str(difference),'absolute_difference':str(abs(difference)),
         'relative_difference':str(difference/central),'relative_percent':str(100*difference/central)}

def compare(packet,fixture):
 require(serialized(fixture)==serialized(read_targets()),'comparison targets or interpretation changed')
 spec=strict_load(MATCHING_SPEC)
 require(serialized(packet.get('specification'))==serialized(spec),'matching specification changed')
 require(serialized(packet.get('scope'))==serialized(spec['scope']),'matching physical scope changed')
 forward=packet['forward']
 require(hashlib.sha256(serialized(forward)).hexdigest()==spec['retained_forward_payload']['sha256'],'matching forward payload changed')
 rows=forward['rows'];inventory=[(r['scale_multiplier'],r['matching_order']) for r in rows]
 require(inventory==[(s,o) for s in (1,2,4) for o in (0,1,2)],'all nine prescribed scale/order rows required')
 require(spec['calculation']['primary_scale_multiplier']==1 and spec['calculation']['primary_GF_order']==2
         and spec['calculation']['primary_boson_order']==2 and spec['calculation']['primary_top']=='method1, pure-QCD order4, remaining order2','predefined primary matching changed')
 target={r['quantity']:r for r in fixture['targets']}
 comparisons=[];primary=[]
 with localcontext() as ctx:
  ctx.prec=32
  for source_row in rows:
   scale,order=source_row['scale_multiplier'],source_row['matching_order']
   entries=[('W','W','BW_mass_GeV'),('Z','Z','BW_mass_GeV'),('H','higgs','mass_GeV'),('GF',None,'GF_local_zero_muon_GeVm2'),('top','top_method0','mass_GeV'),('top','top_method1','mass_GeV')]
   require(('top_QCD4_EW2' in source_row)==(order==2),'extra top control inventory changed')
   if order==2:entries.append(('top','top_QCD4_EW2','mass_GeV'))
   for quantity,method,field in entries:
    field_path=f'{method}/{field}' if method else field
    value=source_row[method][field] if method else source_row[field]
    is_primary=scale==1 and order==2 and (quantity!='top' or method=='top_QCD4_EW2')
    entry={'scale_multiplier':scale,'matching_order':order,'quantity':quantity,'method':method,
           'source_field':field_path,'primary':is_primary,**comparison_value(value,target[quantity])}
    comparisons.append(entry)
    if is_primary:primary.append(copy.deepcopy(entry))
 require(len(comparisons)==57 and [r['quantity'] for r in primary]==['W','Z','H','GF','top'],'complete comparison inventory')
 return {'schema':'oph.source_ew_vev_matching_comparison.v1','status':'DESCRIPTIVE_CONDITIONAL_COMPARISON',
         'primary_selection':fixture['primary_selection'],'primary':primary,'comparisons':comparisons,
         'retained_forward_output':copy.deepcopy(forward),'source_inputs':copy.deepcopy(packet['inputs']),
         'sources':fixture['sources'],'interpretation':fixture['interpretation'],'boundary':fixture['boundary'],
         'joint_significance':None,'theory_uncertainty':None,'fit_or_success_verdict':None}

def build():
 packet=strict_load(PACKET);fixture=read_targets()
 # Check the independent source-boundary, transport and one-loop matching
 # certificate before any empirical target enters the comparison arithmetic.
 import verify_source_ew_vev_matching as verifier
 verifier.verify(packet)
 result=compare(packet,fixture)
 paths=(PACKET,TARGETS,MATCHING_SPEC,PRIOR_TOP,HERE/'verify_source_ew_vev_matching.py',HERE/'verify_source_ew_top_pole.py',Path(__file__))
 result['source_pins']={p.relative_to(ROOT).as_posix():sha(p) for p in paths}
 return result

def verify_payload(payload):
 require(serialized(payload)==serialized(build()),'EW comparison receipt mismatch')
 return {'verified':True,'matching_control_rows':9,'descriptive_comparisons':57,'primary_observables':5,'joint_significance_claimed':False}

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
 rendered=json.dumps(build(),indent=2,sort_keys=True)+'\n'
 if args.check:
  require(OUTPUT.read_text()==rendered,'comparison receipt differs; regenerate')
  print('Source EW descriptive comparison: PASS')
 else:
  OUTPUT.write_text(rendered);print(OUTPUT)
if __name__=='__main__':main()
