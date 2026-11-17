"""Reproducible same-ground scalar/link-loop records and source custody."""
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from ..higgs_composites.certificate import validate_pins,replay_differences
from .model import operator_certificate
from .readout import transfer_checks,complete_readout
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[4]))
    from tools.pin_scope import bound_pins

ROOT=Path(__file__).resolve().parents[4]
CONSUMED='model/coherent_defect_wave/chiral_joint/periodic_spectrum/results.json'


def input_hashes():
    here=Path(__file__).parent
    seeds={str(path.relative_to(ROOT)) for path in here.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(action.SOURCE_FILES);seeds.add(CONSUMED)
    receipt=json.loads((ROOT/CONSUMED).read_text());validate_pins(bound_pins(receipt['input_sha256']));seeds.update(bound_pins(receipt['input_sha256']))
    seeds.update(('model/requirements.txt','model/particle_candidate_bosons/PROOF.md',
                  'model/particle_candidate_scalar_gauge_response/PROOF.md'))
    paths=set(import_closure(seeds))|seeds
    return bound_pins({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)})


def certificate_data():
    return {'schema':'oph.periodic_scalar_link_loop_records.v1',
            'exact_operators_and_instruments':operator_certificate(),'exact_transfers':transfer_checks(),
            'complete_readout':complete_readout(),
            'parameter_ancestry':{'source_H':'unchanged periodic lepton/Higgs source',
                'readout_settings':'local scalar number and two actual closed quantum-link traces',
                'setting_probabilities':'declared equal thirds','empirical_inputs':[]},
            'scope':{'same_interacting_ground_for_scalar_and_loop_channels':True,
                'complete_degenerate_energy_projector_records':True,'elastic_and_failure_records_retained':True,
                'named_physical_bosons_identified':False,'color_spectrum_proved':False,
                'continuum_propagation_or_GeV_mass':False,
                'other_primary_Wilson_source_boson_results_transferred':False}}


def certificate():
    pins=input_hashes();return {**certificate_data(),'input_sha256':bound_pins(pins),'consumed_periodic_receipt_pins_validated':True}


def write():
    result=certificate();Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n');return result
