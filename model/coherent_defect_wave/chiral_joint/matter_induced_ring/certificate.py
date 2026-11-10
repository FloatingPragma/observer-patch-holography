"""Exact ring, full-state records and validated transitive source custody."""
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from ..higgs_composites.certificate import validate_pins
from ..periodic_boson_records.model import operator_certificate
from .model import norm_certificate,split_certificate,ancestry_certificate,resolvent_certificate
from .records import complete_records

ROOT=Path(__file__).resolve().parents[4]
CONSUMED=('model/coherent_defect_wave/chiral_joint/periodic_spectrum/results.json',
          'model/coherent_defect_wave/chiral_joint/periodic_boson_records/results.json')


def input_hashes():
    here=Path(__file__).parent
    seeds={str(path.relative_to(ROOT)) for path in here.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(action.SOURCE_FILES);seeds.update(CONSUMED)
    for name in CONSUMED:
        receipt=json.loads((ROOT/name).read_text())
        validate_pins(receipt['input_sha256']);seeds.update(receipt['input_sha256'])
    seeds.update(('model/requirements.txt','model/coherent_defect_wave/chiral_joint/field_limit/PROOF.md'))
    paths=set(import_closure(seeds))|seeds
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate_data():
    loops=operator_certificate()
    return {'schema':'oph.same_source_matter_induced_ring.v1',
            'source_norm':norm_certificate(),'complete_physical_split':split_certificate(),
            'actual_native_loop_operator_checks':loops,
            'irreducible_ring_ancestry':ancestry_certificate(),
            'exact_energy_dependent_elimination':resolvent_certificate(),
            'complete_source_clock_records':complete_records(),
            'parameter_ancestry':{'idle_Hamiltonian':'unchanged unit-coupling periodic lepton/Higgs source',
                'source_hopping':'actual oriented +i sigma3/4 and -i sigma3/4',
                'projection':'one native ec spin0 at site1; every compatible scalar and link state',
                'resolvent_settings':['-32','-64'],'source_time':'1/16384',
                'empirical_inputs':[],'mass_or_gauge_coupling_calibration':False},
            'scope':{'matter_generated_gauge_invariant_two_link_kernel':True,
                'same_source_complete_state_accounting':True,'nonzero_full_H_loop_records':True,
                'autonomous_pure_gauge_dynamics':False,'two_dimensional_plaquette':False,
                'Coulomb_or_relativistic_phase':False,'physical_photon_or_W_Z_gluon_identified':False,
                'GeV_units':False,'colored_spectrum':False,'literal_Peter_Weyl_119_spectrum_transfer':False}}


def certificate():
    pins=input_hashes()
    return {**certificate_data(),'input_sha256':pins,'consumed_receipt_pins_validated':True}


def write():
    result=certificate()
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
