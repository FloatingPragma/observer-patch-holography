"""Exact positive-witness replay and full-source production records."""
import argparse,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .model import exact_ground,preparation_records
from .quotient import generate_partition,spectrum
ROOT=Path(__file__).resolve().parents[4]

def certificate():
    here=Path(__file__).parent
    seeds={str(p.relative_to(ROOT)) for p in here.iterdir() if p.suffix in ('.py','.md','.npz')}
    seeds.update(('model/requirements.txt','model/coherent_defect_wave/chiral_higgs/PROOF.md',
                  'model/coherent_defect_wave/chiral_joint/periodic_family/PROOF.md'))
    inputs=sorted(seeds|set(import_closure(seeds)))
    return {'schema':'oph.full_native_family_quantum_cell.v1',
        'ground':exact_ground(),'complete_vacuum_energy_records':spectrum(),
        'complete_production_records':preparation_records(),
        'input_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs},
        'scope':{'all_native_quark_and_lepton_Yukawas_active':True,
            'one_native_generation':True,'quantum_scalar_replaced_by_background':False,
            'global_full_carrier_ground_proved':False,'spatial_propagation_or_GeV_mass':False}}

def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    for arg in ('write','check','generate-partition'):mode.add_argument('--'+arg,action='store_true')
    args=parser.parse_args()
    if args.generate_partition:generate_partition()
    result=json.dumps(certificate(),indent=2,sort_keys=True)+'\n';path=Path(__file__).with_name('results.json')
    if args.check:
        if path.read_text()!=result:raise SystemExit('FAIL stale full-family cell receipt')
    else:path.write_text(result)
    print('PASS full-family quantum cell ground and complete production records')

if __name__=='__main__':main()
