"""Exact neutral instruments, complete spectral records and transitive custody."""
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from ..higgs_composites.certificate import validate_pins,replay_differences
from .model import instrument_certificate
from .readout import exact_transfer_checks,star_translation_certificate,complete_readout

ROOT=Path(__file__).resolve().parents[4]
CONSUMED=('model/coherent_defect_wave/chiral_joint/higgs_composites/results.json',
          'model/coherent_defect_wave/chiral_joint/periodic_spectrum/results.json')


def input_hashes():
    here=Path(__file__).parent
    seeds={str(path.relative_to(ROOT)) for path in here.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(action.SOURCE_FILES);seeds.update(CONSUMED)
    for name in CONSUMED:
        receipt=json.loads((ROOT/name).read_text());validate_pins(receipt['input_sha256']);seeds.update(receipt['input_sha256'])
    seeds.add('model/requirements.txt')
    paths=set(import_closure(seeds))|seeds
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate_data():
    return {'schema':'oph.neutral_quantum_higgs_composite.v1',
        'operator':'J_nu,alpha=phi0 c_alpha,11-phi1 c_alpha,10',
        'quantum_scalar_states':5,'native_full_fermion_modes_per_cell':30,
        'finite_instrument':instrument_certificate(),'exact_Gauss_transfers':exact_transfer_checks(),
        'star_selection':star_translation_certificate(),'complete_readout':complete_readout(),
        'parameter_ancestry':{'idle_source':'unchanged frozen two-cell periodic lepton/Higgs Hamiltonian',
            'readout_normalization':'declared success J/2 and positive failure square root',
            'spin_and_site':'declared spin0 and site1 apparatus setting',
            'empirical_inputs':[],'mass_calibration':False},
        'scope':{'gauge_neutral_quantum_composite_readout':True,'same_complete_periodic_ground':True,
            'complete_failures_and_zero_weight_energy_alternatives':True,
            'source_only_translation_phases_are_particle_momenta':False,
            'moving_neutrino_or_relativistic_mass_identified':False,'GeV_units':False,
            'three_observed_neutrino_generations':False,'interacting_continuum':False}}


def certificate():
    pins=input_hashes();return {**certificate_data(),'input_sha256':pins,'consumed_receipt_pins_validated':True}


def write():
    result=certificate();Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n');return result
