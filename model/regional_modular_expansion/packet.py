"""Replay the exact interacting regional-modular certificate."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import expansion as e


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/regional_modular_expansion/expansion.py',
    'model/regional_modular_expansion/packet.py',
    'model/regional_modular_expansion/test_expansion.py',
    'model/regional_modular_expansion/README.md',
    'model/field_symmetry/symmetry.py',
    'model/field_stress/energy.py',
    'model/full_fock_field/field.py',
    'model/interacting_field/lattice.py',
)


def entries(matrix):
    return [[i, j, str(matrix[i, j])] for i in range(matrix.rows)
            for j in range(matrix.cols) if matrix[i, j] != 0]


def numerical_comparison(data, beta=e.BETA):
    """An independent spectral/log replay, explicitly a floating diagnostic."""
    from ..field_thermal.thermal import car_restriction, restrict
    h = data['H_numerator'].matrix(data['scale']).toarray()
    eigenvalues, vectors = np.linalg.eigh(h)
    full = (vectors*np.exp(-float(beta)*eigenvalues))@vectors.conj().T
    reduced = restrict(full, car_restriction(8, e.KEPT))/16
    weights, basis = np.linalg.eigh(reduced)
    if weights.min() <= 0:
        raise ArithmeticError('faithful complete regional reference required')
    k = -(basis*np.log(weights))@basis.conj().T
    polynomial = np.asarray(e.polynomial(data, beta), dtype=complex)
    error = float(np.linalg.norm(k-polynomial, 2))
    bare = k-float(beta)*np.asarray(data['coefficients'][0], dtype=complex)
    tail = e.logarithm_tail(beta, data['bound'])
    if error > float(tail)+1e-13:
        raise ArithmeticError('independent modular replay exceeds the exact tail')
    return {'minimum_unnormalized_reference_weight': float(weights.min()),
            'full_operator_polynomial_error': error,
            'full_minus_empty_bare_difference': float((bare[15, 15]-bare[0, 0]).real),
            'floating_residual_admission': 1e-13,
            'numerics_are_not_the_outward_certificate': True}


def run():
    data = e.plaquette()
    free = e.plaquette(0)
    isolated = e.plaquette(boundary=False)
    c1, c2, c3 = data['coefficients'][:3]
    if not e.centered(c2).is_zero_matrix:
        raise ArithmeticError('the exact quadratic boundary term stopped being scalar')
    commutator = c1*c3-c3*c1
    commutator_squared = e.s.trace(commutator.conjugate().T*commutator)
    if commutator_squared <= 0:
        raise ArithmeticError('the interacting cubic correction should not commute with local energy')
    if any(not e.centered(c).is_zero_matrix for c in isolated['coefficients'][1:]):
        raise ArithmeticError('disconnected Gibbs factors acquired a nonscalar correction')
    proof = e.witness(data)
    if proof['exact_difference_upper'] >= 0 or proof['bare_generator_distance_lower'] <= 0:
        raise ArithmeticError('the exact all-order bare-generator separation was lost')
    encoded = [entries(c) for c in data['coefficients']]
    result = {
        'schema': 'oph.interacting_regional_modular.v1',
        'issues': [1145],
        'gates': ['CP-09'],
        'all_pass': True,
        'domain': {
            'action': 'unchanged full-Fock Weyl/contact H, c=1/8, U=1; only H-I scalar centering',
            'spatial_sites': 4, 'CAR_modes': 8, 'full_dimension': 256,
            'number_sector_dimensions': [1, 8, 28, 56, 70, 56, 28, 8, 1],
            'regional_modes': list(e.KEPT), 'regional_dimension': 16,
            'discarded_dimension': 16, 'norm_bound_exact': str(data['bound']),
            'beta_exact': str(e.BETA), 'general_analytic_domain': 'abs(beta)*B < log(2)',
            'rational_checker_domain': '0 <= beta*B <= 1/2',
            'order': e.ORDER, 'Hamiltonian_model_time_only': True,
            'plaquette_Gibbs_preparation_executed_in_source': False,
        },
        'exact_coefficients': encoded,
        'coefficient_sha256': hashlib.sha256(canonical(encoded).encode()).hexdigest(),
        'leading_correction': {
            'quadratic_is_scalar': True, 'quadratic_scalar': str(e.s.trace(c2)/16),
            'cubic_centered_Pauli': e.pauli_terms(e.centered(c3)),
            'cubic_local_energy_commutator_HS_squared': str(commutator_squared),
            'free_cubic_centered_Pauli': e.pauli_terms(e.centered(free['coefficients'][2])),
            'removed_crossing_bonds_leave_only_scalar_corrections': True,
        },
        'certificate': {key: str(value) for key, value in proof.items()},
        'entropy_tangent': {
            'domain': 'all Hermitian trace-zero regional variations, including noncommuting ones',
            'bound': '|DS_rho(delta)-Tr(delta P8)| <= remainder * ||delta||_1',
            'arbitrary_untouched_reference_extension': 'the modular operator bound tensors with identity without growth',
            'normalization': 'K_A = -log(E exp(-beta(H-I))) + log Tr_A(E exp(-beta(H-I))) I; scalar drops only on trace-zero tangents',
        },
        'independent_numerical_replay': numerical_comparison(data),
        'scope': {
            'closed': 'exact same-action regional modular polynomial and uniform all-state first-law error; certified bare local-Gibbs mismatch',
            'remaining': [
                'Geometric cap weights, Lorentz stress/null translations and a common physical clock.',
                'Edge/area normalization, fixed-proper-volume generalized-entropy stationarity and positive Newton scale.',
                'An interacting continuum/scaling theorem and a joint Einstein-register inhabitant.',
                'Source preparation of this full plaquette Gibbs reference; finite action compilation alone does not prepare it.',
            ],
            'no_Einstein_or_physical_particle_gate_accepted': True,
        },
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                          for path in SOURCE_PATHS},
    }
    return result


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
    print(json.dumps({key: result[key] for key in ('all_pass', 'domain', 'certificate', 'scope')}, indent=2))


if __name__ == '__main__':
    main()
