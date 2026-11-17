import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_larmor_dynamics.certificate import equivalent
from .source import field_certificate,block,upper,exact_geometry
from .propagation import refined_records,coefficients
from .spectrum import spectra
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def pins():
 seeds={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md','.tex')};seeds.update(SOURCE_FILES);seeds.update(('model/openwave_dynamic_wave/PROOF.md','model/openwave_dynamic_stress/PROOF.md'))
 return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}
def certificate():
 source_errors={f'{i},{j}':str(upper(block(i,j)['error']))for i,j in((0,0),(1,0),(0,1),(1,1))};c=coefficients();records=refined_records();spectrum=spectra()
 assert records['master_same_source_arrival_delay_interval'][0]>0 and records['master_time4_bending_interval'][1]<0
 assert records['certified_compression_comparisons'][-1]['full_arrival_instrument_reference_stable_TV_upper']<.000358
 return {'schema':'oph.fixed_master_geometry_refinement.v1','verdict':'PASS_FIXED_GEOMETRY_MASTER_AND_COMPLETE_NESTED_REFINEMENT',
  'construction':'cap4 of each of two native harmonic geometry modes; the declared full source is fixed before its nested occupation compressions',
  'same_matter_field_reference_and_observer_source':True,'full_source_dimension':80000,'native_sector_count':64,'geometry_dimension':25,'wave_and_vacuum_dimension':50,
  'complete_native_source_energy_bound':32,'complete_readout_work_absolute_upper':32,
  'source_H_assembly_error_upper_by_native_block':source_errors,'all25_exact_positive_wave_root_certificates':field_certificate()['root_certificates'],
  'Chebyshev_propagator':{'step':'1/2','normalization':32,'degree':64,'coefficient_intervals':c['coefficient_intervals'],'tail_upper':str(c['tail']),
   'arithmetic_contract':'binary64 round-to-nearest, no overflow; explicit gamma_n bounds for the real-matrix/complex-vector recurrences, assembly and complete records'},
  'spectra':spectrum,**records,
  'infinite_geometry_cap_or_spatial_continuum_limit_proved':False,'physical_G_seconds_or_Einstein_limit_identified':False,
  'scope':'Finite refinement within one fixed cap-four construction. Re-exponentiating a smaller coordinate changes the action; the nested comparison compresses the fixed action and retains its omitted-output labels. Source-parameter labels index a finite event word and are not a fundamental global clock.'}
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();before=pins();r=certificate();assert pins()==before;r['source_sha256']=before;path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:assert equivalent(json.loads(path.read_text()),r)
 print(r['verdict'])
if __name__=='__main__':main()
