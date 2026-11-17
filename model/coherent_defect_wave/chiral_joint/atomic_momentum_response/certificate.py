"""Exact momentum records and fail-closed source/reference custody."""
import hashlib
import json
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from ..higgs_composites.certificate import validate_pins
from ..periodic_family.model import Configuration
from ..peter_weyl_links.model import descriptor
from .model import REFERENCE_RECEIPT
from .pump import GROUND_RESPONSE_RECEIPT, full_local_operator_checks, actual_composite_translation_checks, complete_pump_records
from .checks import native_translation_join, character_checks, complete_reference_records, rotation_translation_join

ROOT = Path(__file__).resolve().parents[4]
PACKAGE = 'model/coherent_defect_wave/chiral_joint/atomic_momentum_response/'


def input_hashes():
    seeds = {PACKAGE + name for name in ('__init__.py', '__main__.py', 'model.py',
        'checks.py', 'pump.py', 'certificate.py', 'test_momentum.py', 'README.md', 'PROOF.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/local_reference/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/ground_composite_response/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/peter_weyl_links/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/periodic_family/PROOF.md'))
    for name in (REFERENCE_RECEIPT, GROUND_RESPONSE_RECEIPT):
        path = ROOT / name
        if not path.is_file():
            raise FileNotFoundError('missing consumed atomic momentum receipt: ' + name)
        receipt = json.loads(path.read_text())
        validate_pins(receipt['input_sha256'])
        seeds.add(name)
        seeds.update(receipt['input_sha256'])
    for name in sorted(seeds):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing declared atomic momentum input: ' + name)
    paths = import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT / name).is_file():
            raise FileNotFoundError('missing atomic momentum import input: ' + name)
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    return {'schema': 'oph.atomic_reference_source_and_joint_momentum.v1',
        'actual_full_PW_sources': [descriptor(Configuration(side=side)) for side in (2, 4, 6)],
        'source_and_atomic_reference_join': native_translation_join(),
        'actual_cubic_translation_semidirect_action': rotation_translation_join(),
        'size_dependent_complete_source_momentum': character_checks(),
        'complete_purifying_reference_records': complete_reference_records(),
        'complete_neutral_pump_operator': full_local_operator_checks(),
        'actual_neutral_operator_translation_join': actual_composite_translation_checks(),
        'complete_unconditional_pump_records': complete_pump_records(),
        'clock_and_apparatus': {
            'source_clock': 'exp(-it H_PW), full unit hopping and all onsite source terms',
            'static_reference_clock': 'identity on every retained reference',
            'initial_source_density': 'tensor_x(P_local_ground/2) tensor all blank quantum links',
            'source_reference_preparation': 'product of normalized local ground/reference spin singlets',
            'preparation_service': 'declared exact finite isometry on the full physical carrier',
            'phase_readout': 'complete commuting translation spectral projectors, copied to orthogonal records',
            'pump_completion': 'finite Halmos unitary for the stated rectangular contraction',
            'native_local_gate_schedule_for_preparation_and_readout_proved': False,
            'energy_or_time_window_used_to_select_momentum': False,
            'all_momentum_records_hold_for_every_real_source_time': True},
        'scope': {
            'actual_three_axis_interacting_source': True,
            'all_native_quark_lepton_scalar_link_terms_retained': True,
            'static_references_are_part_of_complete_apparatus': True,
            'source_and_joint_momentum_identified': False,
            'atomic_product_used_as_full_interacting_stationary_vacuum': False,
            'energy_momentum_dispersion_or_named_particle_identification': False,
            'continuum_thermodynamic_or_GeV_claim': False},
        'input_sha256': input_hashes()}
