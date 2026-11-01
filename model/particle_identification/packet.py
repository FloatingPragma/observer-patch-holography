"""Replay exact scoped particle identification tests without mass fitting."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import audit, energy, vacuum

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/particle_identification/audit.py',
    'model/particle_identification/energy.py',
    'model/particle_identification/vacuum.py',
    'model/particle_identification/packet.py',
    'model/particle_identification/test_audit.py',
    'model/particle_identification/test_energy.py',
    'model/particle_identification/test_vacuum.py',
    'model/particle_identification/README.md',
    'paper/tex_fragments/PARTICLE_IDENTIFICATION.tex',
    'model/cube_eta_protection/cube.py',
    'model/cube_eta_protection/repair.py',
    'model/cube_eta_protection/spectrum.py',
    'model/field_symmetry/symmetry.py',
    'model/moving_bound_pair/limit.py',
    'model/moving_bound_pair/README.md',
    'model/particle_comparison/results.json',
    'model/particle_comparison/observations.json',
    'model/particle_comparison/README.md',
)
SOURCE_PATHS = tuple(dict.fromkeys((*SOURCE_PATHS, *energy.CALIBRATION_SOURCES)))


def comparison_ancestry():
    """Audit retained known-data status; do not rerun or refit observations."""
    previous = strict_json.loads((ROOT/'model/particle_comparison/results.json').read_text())
    observations = strict_json.loads((ROOT/'model/particle_comparison/observations.json').read_text())
    assert observations['exposure'] == 'retrospective_known_data'
    assert previous['CP11']['eligible_physical_comparisons_added'] == 0
    assert previous['CP11']['eligible_physical_dictionary'] is False
    assert previous['two_candidate_assignment']['independently_predicted_absolute_masses'] == 0
    assert previous['particle_identity_verdict']['electron_positron_identification_established'] is False
    return {'exposure': observations['exposure'],
            'previous_mass_scale': 'one supplied electron mass calibration',
            'independently_predicted_absolute_masses': 0,
            'previous_dimer_equal_gap_test': 'conditional retrospective equality; not physical identity',
            'previous_dimer_excluded_assignments': previous['three_lepton_exhaustive_comparison']['total_cases'],
            'previous_dimer_exclusion_applied_to_new_cube': False,
            'new_observations_or_mass_fits': 0,
            'new_eligible_physical_comparisons': 0,
            'physical_species_assignments': 0}


def run():
    gap = audit.low_momentum_gap(1, F(1, 4), F(1, 10), F(1, 10))
    return {'schema': 'oph.particle_identification.v1',
            'all_pass': True, 'gates': ['CP-07', 'CP-12'],
            'identified_model_energies_and_inertial_mass': energy.certificate(),
            'grading': audit.grading_certificate(),
            'vacuum': vacuum.certificate(),
            'same_number_energy_origin': audit.neutral_number_shift_certificate(),
            'mass_shell': audit.shell_certificate(),
            'nonrelativistic_control_example': {k: str(v) if isinstance(v, F) else v for k, v in gap.items()},
            'comparison_ancestry': comparison_ancestry(),
            'physical_target_sources': {
                'electron_spin': 'https://pdg.lbl.gov/2025/listings/rpp2025-list-electron.pdf',
                'quark_fermion_quantum_numbers': 'https://pdg.lbl.gov/2025/reviews/rpp2025-rev-quark-model.pdf',
                'Lorentz_mass_shell': 'https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf'},
            'remaining_identity_requirements': [
                'physical vacuum and common Hamiltonian/clock identification',
                'physical electric and other gauge charges with current couplings',
                'Lorentz spin/helicity, chirality and appropriate exchange law',
                'propagating mass shells and stable or resonant pole interpretation',
                'common units independent of fitting the target species mass',
                'complete interactions and independent dimensionless physical comparisons'],
            'scope': {'accepted_physical_gates_added': 0,
                      'whole_OPH_exclusion': False,
                      'other_odd_excitations_or_scalings_excluded': False,
                      'nonrelativistic_composite_interpretation_excluded': False,
                      'source_readable_laws_admitted': True},
            'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        OUTPUT.write_text(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text()), result)
    print(json.dumps({k: result[k] for k in ('all_pass', 'comparison_ancestry', 'scope')}, indent=2))


if __name__ == '__main__':
    main()
