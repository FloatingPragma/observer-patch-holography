"""Fail-closed source custody and exact finite clock certificates."""
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from ..higgs_composites.certificate import validate_pins
from .model import polynomial_certificate, peter_weyl_transfer, family_transfer, rational_ring_records

ROOT = Path(__file__).resolve().parents[4]
CONSUMED = (
    'model/coherent_defect_wave/chiral_joint/peter_weyl_links/results.json',
    'model/coherent_defect_wave/chiral_joint/matter_induced_ring/results.json',
    'model/coherent_defect_wave/chiral_joint/family_transport/results.json',
)


def input_hashes():
    here = Path(__file__).parent
    seeds = {str(path.relative_to(ROOT)) for path in here.iterdir()
             if path.suffix in ('.py','.md')}
    seeds.update(CONSUMED)
    for name in CONSUMED:
        receipt = json.loads((ROOT/name).read_text())
        validate_pins(receipt['input_sha256'])
        seeds.update(receipt['input_sha256'])
    paths = set(import_closure(seeds)) | seeds
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def certificate():
    pins = input_hashes()
    return {'schema':'oph.same_H_finite_rational_clock.v1',
            'unitary_clock_theorem':polynomial_certificate(),
            'Peter_Weyl_neutral_arrival':peter_weyl_transfer(),
            'five_family_local_contrast':family_transfer(),
            'complete_periodic_ring_records':rational_ring_records(),
            'architecture':{
                'service':'declared exact algebraic finite matrix functional calculus',
                'one_common_original_H_per_experiment':True,
                'physical_H_shift_or_log_replacement':False,
                'efficient_local_circuit_for_huge_source_proved':False,
                'empirical_inputs':[],
            },
            'physical_scope':{'finite_unitary_records':True,
                'relativistic_particle_identification':False,'physical_GeV_scale':False},
            'input_sha256':pins, 'consumed_receipt_pins_validated':True}
