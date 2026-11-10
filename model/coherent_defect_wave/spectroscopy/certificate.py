"""Source-bound exact finite phase records and their physical interpretation."""
from fractions import Fraction as F
from pathlib import Path
import hashlib

from model.particle_common_contract.provenance import import_closure
from .readout import source_experiment, fourier_branches, filter_clock_branches
from .readout import energy_branch_control, numerical_free_phase_bound, resolution_plan
from .background_readout import background_experiment

ROOT = Path(__file__).resolve().parents[3]


def input_hashes():
    local = Path(__file__).parent
    paths = {p for p in local.glob('*.py') if not p.name.startswith('test_')}
    paths.update(local/name for name in ('PROOF.md', 'README.md'))
    paths.update(ROOT/name for name in (
        'model/requirements.txt',
        'model/particle_common_contract/provenance.py',
        'model/particle_primary_joint_spectroscopy/PROOF.md',
        'model/spin_shell_readout/PROOF.md',
        'Lean/QFT/FiniteUnitaryScatteringNoGo.lean',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    ))
    paths = {ROOT/name for name in import_closure({str(p.relative_to(ROOT)) for p in paths})}
    paths = {p for p in paths if p.name in {q.name for q in p.parent.iterdir()}}
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def certificate():
    observer, initial, preparation, powers, branches, moments = source_experiment()
    rows = [{'gauss_preparation': -1, 'terminal': 'preparation_failure', 'probability': '1/2'}]
    for label, state in branches:
        rows.append({'gauss_preparation': 1, **label,
                     'probability': str(state.norm_squared()),
                     'positive_constraint_preserved': state.gauss_positive()})
    plain = [state.norm_squared()*2 for _, state in fourier_branches(powers[:4])]
    active = [F(row['probability'])*F(288, 17) for row in rows
              if row.get('filter') == 'nonflat_filter']
    # Independent moment formula for the four-bin polynomial kernel.
    exact = [F(1, 4)+(3*moments[1]+2*moments[2]+moments[3])/8,
             F(1, 4)-moments[2]/4,
             F(1, 4)+(-3*moments[1]+2*moments[2]-moments[3])/8,
             F(1, 4)-moments[2]/4]
    if plain != exact or sum(F(row['probability']) for row in rows) != 1:
        raise ArithmeticError('source vectors disagree with complete phase-kernel moments')
    identity = filter_clock_branches((preparation[0][1],)*6)
    if any(state.norm_squared() for label, state in identity if label['filter'] == 'nonflat_filter'):
        raise ArithmeticError('identity-clock countercontrol populated the nonflat filter')
    return {
        'schema': 'oph.coherent_defect_wave.spectroscopy.v1',
        'source_program': 'coherent-defect-wave-v1',
        'input_sha256': input_hashes(),
        'source_dimension': observer.dimension,
        'positive_gauss_dimension': observer.dimension//2,
        'initial_preparation': '|link=1> tensor A|port0,internal3>',
        'preparation_success_and_failure': ['1/2', '1/2'],
        'clock': {'unit': 'one unchanged joint U tick', 'bins_executed': 4,
                  'phase_centers_in_quarter_turns': [0, 1, 2, 3],
                  'maximum_actual_U_power': 5,
                  'minimum_retained_pointer_labels': 9,
                  'continuous_time_generator_selected': False},
        'normalized_Gauss_positive_moments': list(map(str, moments)),
        'source_power_support_sizes': [len(state.amplitudes) for state in powers],
        'plain_clock_conditional_probabilities': list(map(str, plain)),
        'nonflat_filter_probability_given_positive_Gauss': '17/144',
        'nonflat_filter_original_trial_probability': '17/288',
        'nonflat_clock_conditional_probabilities': list(map(str, active)),
        'complete_original_trial_records': rows,
        'nonflat_frequency_support': numerical_free_phase_bound(moments),
        'identity_clock_has_zero_nonflat_weight': True,
        'resolution_plan_not_executed': resolution_plan(),
        'positive_energy_branch_countercontrol': energy_branch_control(),
        'stationary_background_relative_readout': background_experiment(),
        'scope': {'same_U_complete_phase_instrument_executed': True,
                  'both_exact_flat_phases_removed_in_active_filter': True,
                  'positive_original_nonflat_phase_weight': True,
                  'positive_constraint_preserved_by_every_success_branch': True,
                  'first_experiment_stationary_reference_preparation': False,
                  'explicit_stationary_reference_in_second_experiment': True,
                  'physical_ground_state_preparation': False,
                  'conserved_momentum_identified': False,
                  'physical_energy_branch_selected': False,
                  'isolated_particle_pole_or_threshold_identified': False,
                  'relativistic_mass_or_GeV_value': False,
                  'Standard_Model_species': False},
    }
