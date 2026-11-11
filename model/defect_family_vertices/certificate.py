"""Full native contact, scalar pulse and neutral-channel receipt."""
import hashlib
from pathlib import Path

from .source import ROOT, QUARK_RECEIPT, LEPTON_RECEIPT, encoded_value, ground_pair, scalar_readout
from .vertices import contacts, pulse_responses, scalar_double_commutator
from .neutral import calculate as neutral, RECEIPT as NEUTRAL_RECEIPT


def source_bindings():
    from model.particle_common_contract.provenance import import_closure
    here = Path(__file__).parent
    seeds = {path.relative_to(ROOT).as_posix() for path in here.iterdir()
             if path.suffix in ('.py', '.md')}
    seeds.update({QUARK_RECEIPT, LEPTON_RECEIPT, NEUTRAL_RECEIPT,
        'model/defect_scalar_response/results.json',
        'model/coherent_defect_wave/chiral_joint/local_reference/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/ground_composite_response/PROOF.md',
        'model/observer_particles/neutrino_action_extension.py'})
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(seeds))}


def certificate():
    from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import field_data
    bindings = source_bindings()
    ground = ground_pair()
    value = {'schema': 'oph.defect_family_vertices.v1',
        'construction': {'local_CAR_modes': 30, 'scalar_states': 5,
            'primitive_coefficients': {'up': '1', 'down': '1', 'charged_lepton': '1'},
            'reference': 'actual local physical Gauss-ground projector divided by two',
            'clock': 'exp(-it H_local)', 'source_dimension': str(5*2**30),
            'ground_energy_minimal_polynomial': str(field_data()[1].as_expr()),
            'ground_energy_selected_root': 'smallest real root',
            'ground_energy': encoded_value(ground['field'].unit),
            'squared_transfer_eigenvalue': encoded_value(ground['lam']),
            'three_family_source_constructed': False, 'measured_particle_inputs': []},
        'complete_scalar_readout': scalar_readout(),
        'independent_scalar_commutator': scalar_double_commutator(),
        'native_contact_vertices': contacts(),
        'actual_scalar_pulse_responses': pulse_responses(),
        'neutral_channels': neutral(),
        'complete_apparatus': {
            'state_selection': 'declared actual physical ground preparation law',
            'preparation_success_probability_from_arbitrary_input_supplied': False,
            'charged_probe_references_retained': True,
            'scalar_pointer_effects_complete': True,
            'scalar_no_operation_retained': True,
            'all_spectral_zero_residues_retained': True,
            'source_ground_is_unconstrained_or_boundary_compensated_ground': False},
        'scope': {'actual_scalar_insertions_and_full_CAR_contacts': True,
            'source_specific_spectral_normalization': True,
            'scalar_cost_susceptibility_relabelled_as_Higgs_vertex': False,
            'frequency_dependent_1PI_or_on_shell_vertex': False,
            'physical_particle_masses_or_neutrino_oscillations': False,
            'family_hierarchy_fitted': False}}
    if source_bindings() != bindings:
        raise ValueError('source changed during native family calculation')
    value['input_sha256'] = bindings
    return value
