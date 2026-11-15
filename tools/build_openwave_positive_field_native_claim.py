#!/usr/bin/env python3
"""Bind the independently accepted finite positive-Maxwell composite source.

This source-bound consumer requires the sibling OpenWave workspace; it is not
advertised as a standalone replay of the enormous native program.
"""
from pathlib import Path
from fractions import Fraction
import argparse,hashlib,json
R=Path(__file__).resolve().parents[1];W=R.parent
RECEIPT='openwave-tuning/row_completion/ROW_13/POSITIVE_FIELD_VALIDATION.json'
PIN='d382fba4e62d7e82d5cf3445b8f902a9a82f21209cb22b9aec68aeec15660749'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
 p=W/RECEIPT;assert sha(p)==PIN;v=json.loads(p.read_text())
 assert v['row']==13 and v['full_criterion_met'] and not v['full_native_word_executed']
 assert len(v['independent_reviews'])==2 and all(r['producer_independent'] and r['full_criterion_met'] for r in v['independent_reviews'])
 for n,h in v['evidence_files'].items():assert sha(W/n)==h,n
 tv=Fraction(v['complete_native_TV_exact']);weak=Fraction(v['weak_complete_bad_probability_upper_exact']);unforced=Fraction(v['unforced_complete_bad_probability_upper_exact'])
 assert weak>=Fraction(2,25)+tv and weak<Fraction(81,1000)
 assert unforced>=Fraction(7,100)+tv and unforced<Fraction(71,1000)
 return {'schema':'oph.openwave.positive_field_native_claim.v1','claim':'Stable neutral three-quark composite coupled to the actual positive Maxwell field, with complete finite native preparation and readout bounds.','status':'Internally accepted source-specific full ROW13 alternative; official OpenWave status unchanged.','literal_criterion':v['literal_criterion'],'criterion_sha256':v['criterion_sha256'],'full_criterion_met':True,'internal_row':13,'internal_row_count_increment':0,'official_rows_changed':[],
 'acceptance':{'workspace_path':RECEIPT,'sha256':PIN,'primary_receipt_preserved':v['primary_receipt_preserved']},'source_domain':v['source_domain'],'preparation':v['preparation'],'clock':v['source_clock'],'source_energy_unit':v['source_energy_unit'],'source_stability_all_time':True,'complete_native_time_window':v['complete_native_time_window'],'native_probe_domain':'Linear profile xi(x)=x along the declared path, |f|<=1/4.',
 'complete_native_TV_exact':v['complete_native_TV_exact'],'complete_native_TV_display':float(tv),'unforced_complete_bad_probability_upper_exact':str(unforced),'unforced_complete_bad_probability_upper_display':float(unforced),'weak_complete_bad_probability_upper_exact':str(weak),'weak_complete_bad_probability_upper_display':float(weak),'outward_rounding':v['probability_bound_rounding'],
 'independent_reviews':v['independent_reviews'],'raw_native_pulses_per_independent_review':v['raw_native_pulses_per_independent_review'],'full_native_pulses_upper':v['full_native_pulses_upper'],'full_native_word_executed':False,'scope_exclusions':v['scope_exclusions'],'all31_same_source_proved':False,'evidence_file_count':len(v['evidence_files']),'evidence_files':v['evidence_files'],
 'replay_scope':'This builder verifies exact source-bound evidence and the accepted composition. Independent scientific replays are in the referenced frozen review packages. Requires the sibling openwave-tuning workspace; no standalone full-native replay claim.','portable_scientific_closure_in_RER':False}
def main():
 a=argparse.ArgumentParser();g=a.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');opts=a.parse_args();d=build();p=R/'claims/openwave_positive_field_native.json';b=json.dumps(d,indent=2,sort_keys=True)+'\n'
 if opts.write:p.write_text(b)
 else:assert p.read_text()==b,'Claim differs from exact accepted receipt'
 print('POSITIVE_FIELD_NATIVE_CLAIM_PASS',d['evidence_file_count'],d['complete_native_TV_display'])
if __name__=='__main__':main()
