"""Exact square records and validated complete Peter-Weyl custody."""
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from ..higgs_composites.certificate import validate_pins
from .model import BLOCKS,geometry_certificate,physical_checks,leading_histories
from .bounds import feshbach_bound,complete_records

ROOT=Path(__file__).resolve().parents[4]
CONSUMED='model/coherent_defect_wave/chiral_joint/peter_weyl_links/results.json'


def input_hashes():
    here=Path(__file__).parent
    seeds={str(path.relative_to(ROOT)) for path in here.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(action.SOURCE_FILES);seeds.add(CONSUMED)
    source=json.loads((ROOT/CONSUMED).read_text());validate_pins(source['input_sha256'])
    seeds.update(source['input_sha256']);seeds.add('model/requirements.txt')
    paths=set(import_closure(seeds))|seeds
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate_data():
    return {'schema':'oph.same_PW_source_matter_induced_plaquette.v1',
            'geometry_and_source':geometry_certificate(),
            'physical_preparations_and_loop_operators':[physical_checks(name) for name in BLOCKS],
            'complete_minimal_histories':[leading_histories(name) for name in BLOCKS],
            'full_carrier_matter_elimination':[feshbach_bound(name) for name in BLOCKS],
            'complete_original_clock_records':complete_records(),
            'parameter_ancestry':{'idle_source':'frozen complete Peter-Weyl Hamiltonian, all unit coefficients',
                'settings':'two native blocks, equal recorded probability; same contractible square',
                'preparation_and_reference':'declared finite invariant controller, full static dual-R15 CAR reference',
                'empirical_inputs':[],'physical_gauge_coupling_or_mass_calibration':False},
            'scope':{'matter_generated_contractible_loop_on_full_quantum_source':True,
                'all_native_gauge_currents_used':True,'all_source_interactions_in_full_remainder':True,
                'added_primitive_magnetic_term':False,'autonomous_pure_gauge_Hamiltonian':False,
                'Coulomb_confining_or_relativistic_phase':False,
                'physical_photon_gluon_W_or_Z_identified':False,'GeV_units':False,
                'source_ground_or_mass_spectrum_imported':False}}


def certificate():
    pins=input_hashes()
    return {**certificate_data(),'input_sha256':pins,'consumed_receipt_pins_validated':True}


def write():
    result=certificate()
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
