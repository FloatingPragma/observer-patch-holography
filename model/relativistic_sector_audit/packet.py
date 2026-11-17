"""Replay the source-bound doubled Weyl sector and its exact error budgets."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import sector as s

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/relativistic_sector_audit/__init__.py',
    'model/relativistic_sector_audit/__main__.py',
    'model/relativistic_sector_audit/sector.py',
    'model/relativistic_sector_audit/packet.py',
    'model/relativistic_sector_audit/test_sector.py',
    'model/relativistic_sector_audit/README.md',
    'paper/tex_fragments/RELATIVISTIC_SECTOR_AUDIT.tex',
    'model/interacting_field/lattice.py',
    'model/interacting_field/DERIVATION.md',
    'model/periodic_momentum/momentum.py',
    'model/dilute_momentum_ward/ward.py',
    'model/dilute_momentum_ward/README.md',
    'model/contact_scattering/onshell.py',
    'model/coarse_energy_ward/ward.py',
    'model/lattice_momentum_balance/momentum.py',
)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(x) for key, x in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(x) for x in value]
    return value


def run():
    budgets = [encode(s.continuum_budget(F(1, denominator), 1, 2, 1))
               for denominator in (8, 16, 32, 64, 128)]
    example = s.continuum_budget(F(1, 64), 1, 2, 1)
    if example['complete_vector_error_upper'] >= F(1237, 10**6):
        raise ArithmeticError('advertised continuum example bound failed')
    return {'schema': 'oph.relativistic_sector_audit.v1', 'all_pass': True,
            'issues': [1142, 1143, 1145], 'gates': ['CP-06', 'CP-07', 'CP-09'],
            'domain': {'action': 'H_a=dGamma[-(v/a) sigma.sin(k)]+U(a)D on the unchanged infinite lattice',
                      'declared_scaling': 'c=v/a, v>0 fixed, finite U(a)>=0',
                      'momentum_domain': 'eight disjoint nu+a[-Lambda,Lambda]^3; a Lambda<pi/2',
                      'calculator_domain': 'exact rational positive a,Lambda,v; a Lambda<=1; Nmax<=4096',
                      'state_domain': 'all N<=Nmax in initial spin-complete eight-valley band, arbitrary untouched references',
                      'no_energy_shift_or_time_reversal': True,
                      'no_projection_during_evolution': True},
            'source_valleys': s.valley_certificate(),
            'relativistic_algebra': s.relativistic_algebra_certificate(),
            'theorem': {'one_particle_embedding': '(W_a f)(nu+a p)=a^(-3/2)f_nu(p), zero outside eight patches',
                        'Fock_embedding': 'direct sum of wedge^n W_a for 0<=n<=Nmax',
                        'continuum_generator': 'dGamma(direct sum_nu -v sigma.diag((-1)^nu) p)',
                        'contact_fraction': 'f=8(a Lambda/pi)^3 in infinite-lattice Haar measure',
                        'contact_bound': '||D Omega||<=floor(Nmax/2)*sqrt(f), countable positive-form proof',
                        'full_dynamical_error': 'min(2,abs(t)*(Nmax*v*a^2*Lambda^3/2+abs(U)*floor(Nmax/2)*sqrt(f)))',
                        'projector_error': 'a^2 Lambda^2/6 away from p=0',
                        'coupling_limit_domain': 'U(a)=o(a^(-3/2)); Nmax,Lambda,v and compact time interval fixed',
                        'N0_N1': 'contact identically zero for every finite U(a)',
                        'limiting_mass': '0 in all eight valleys',
                        'limiting_interaction': 'free on the proved domain'},
            'source_unit_energy_table': s.energy_table(),
            'finite_two_particle_contact_control': s.finite_contact_check(),
            'rational_convergence_examples': budgets,
            'reported_example_vector_error_upper': '1237/1000000',
            'physical_boundary': {'physical_species_assignments': 0,
                                  'independently_predicted_SI_masses': 0,
                                  'massive_relativistic_sector_constructed': False,
                                  'positive_energy_interacting_vacuum_constructed': False,
                                  'finite_density_limit_constructed': False,
                                  'full_Fock_uniform_limit_constructed': False,
                                  'source_action_or_readout_laws_inferred_from_uniqueness': False,
                                  'new_native_preparation_or_detector_execution': False,
                                  'existing_negative_branch_retained': True,
                                  'protected_cube_repair_invariance_transferred': False,
                                  'protected_cube_response_lines_identified_with_valleys': False,
                                  'critical_pair_scaling_covered': False,
                                  'physical_Einstein_dynamics_constructed': False},
            'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        OUTPUT.write_text(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text()), result)
    print(json.dumps({key: result[key] for key in ('all_pass', 'source_unit_energy_table', 'physical_boundary')}, indent=2))


if __name__ == '__main__':
    main()
