"""Replay exact conserved periodic momenta and their finite-range boundary."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as s

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..full_fock_field.field import FieldGeometry
from ..coarse_energy_ward import ward as w
from ..particle_response.dispersion import PeriodicSpectrum
from . import momentum as m

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/periodic_momentum/momentum.py', 'model/periodic_momentum/packet.py',
    'model/periodic_momentum/test_momentum.py', 'model/periodic_momentum/README.md',
    'model/coarse_energy_ward/ward.py', 'model/lattice_momentum_balance/momentum.py',
    'model/full_fock_field/field.py', 'model/interacting_field/lattice.py',
    'model/particle_response/dispersion.py', 'model/bound_pair/protection.py',
    'Lean/QFT/PeriodicMomentum.lean',
)


def encoded(matrix):
    return [[i, j, str(v)] for (i, j), v in sorted(m.cleaned(matrix).todok().items())]


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def torus_certificate():
    geometry = m.PeriodicGeometry()
    h, h0, v = m.action(geometry)
    shifts = tuple(geometry.shift(axis) for axis in range(3))
    for axis, permutation in enumerate(shifts):
        for operator in (h, h0, v):
            need(m.permute_clifford(operator, permutation) == operator,
                 'translation changed the actual periodic source action')
        iterate = tuple(range(geometry.modes))
        for _ in range(geometry.sides[axis]):
            iterate = tuple(permutation[j] for j in iterate)
        need(iterate == tuple(range(geometry.modes)), 'wrong translation order')
    for p in shifts:
        for q in shifts:
            need(tuple(p[j] for j in q) == tuple(q[j] for j in p),
                 'Cartesian translations failed to commute')
    payload = [[str(mask), str(a), str(b)] for mask, (a, b) in sorted(h.terms.items())]
    return {'periodic_sides': list(geometry.sides), 'CAR_modes': geometry.modes,
            'full_Fock_dimension': 2**geometry.modes, 'action_majorana_terms': len(h.terms),
            'exact_translation_covariances': 9, 'commuting_translation_pairs': 9,
            'translation_orders': list(geometry.sides),
            'three_conserved_principal_spectral_logs': True,
            'complete_Fock_coverage_is_symbolic_not_state_execution': True,
            'action_sha256': hashlib.sha256(canonical(payload).encode()).hexdigest()}


def ring_data():
    geometry = m.PeriodicGeometry((1, 1, 4))
    h, h0, v = m.action(geometry)
    shift = m.graded_matrix(geometry.shift(2))
    p, projectors = m.spectral_log(shift, 4)
    hm, vm = m.clifford_matrix(h, 8), m.clifford_matrix(v, 8)
    one = [1 << j for j in range(8)]
    q = p.extract(one, one)
    additive = m.lift_matrix(q)
    winding = m.cleaned((additive-p)/2)
    return geometry, shift, p, projectors, hm, vm, q, additive, winding


def ring_certificate(data=None):
    geometry, shift, p, projectors, h, v, q, additive, winding = ring_data() if data is None else data
    identity = s.eye(256)
    need(shift**4 == identity, 'graded translation order failed')
    need(m.is_zero(p-p.conjugate().T) and m.is_zero(h*p-p*h),
         'principal momentum is not Hermitian and conserved')
    need(m.is_zero(sum(projectors, s.zeros(256))-identity), 'lost translation sectors')
    for j, projector in enumerate(projectors):
        need(m.is_zero(projector-projector.conjugate().T), 'non-Hermitian translation projector')
        need(m.is_zero(shift*projector-s.I**j*projector), 'wrong exact translation eigenvalue')
        need(m.is_zero(p*projector-F(m.branch_labels(4)[j], 2)*projector),
             'wrong exact principal logarithm label')
        for k, other in enumerate(projectors):
            need(m.is_zero(projector*other-(projector if j == k else s.zeros(256))),
                 'translation spectral projectors are not complete orthogonal projections')
    negative_cut = m.cleaned(p-2*projectors[2])
    need(m.is_zero(negative_cut*projectors[2]+projectors[2]) and
         m.is_zero(h*negative_cut-negative_cut*h), 'alternate -pi branch failed')
    need(m.is_zero(winding**3-winding), 'winding operator lost its integer spectrum')
    need(m.is_zero(h*additive-additive*h-2*(h*winding-winding*h)),
         'contact winding accounts for the additive-log force incorrectly')
    labels = [j for j in m.branch_labels(4) for _ in range(2)]
    histogram = Counter((sum(labels[j] for j in range(8) if mask >> j & 1)-
                         m.centered_integer(sum(labels[j] for j in range(8) if mask >> j & 1), 4))//4
                        for mask in range(256))
    for k in (-1, 0, 1):
        projector = (winding*winding+k*winding)/2 if k else identity-winding*winding
        need(s.expand(s.trace(projector)) == histogram[k], 'winding multiplicity mismatch')
    initial, final = m.fourier_pair(4, 1, 1), m.fourier_pair(4, -1, -1)
    need((initial.conjugate().T*initial)[0] == 1 and (final.conjugate().T*final)[0] == 1,
         'Umklapp pair states are not normalized')
    vertex = s.expand((final.conjugate().T*v*initial)[0])
    force = s.expand((final.conjugate().T*(h*additive-additive*h)*initial)[0])
    need(vertex == F(1, 4) and force == F(1, 2), 'exact contact Umklapp witness failed')
    for vector in (initial, final):
        need(m.is_zero(p*vector-vector), 'Umklapp changed total crystal momentum')
    need(m.is_zero(additive*initial-initial) and m.is_zero(additive*final+final),
         'additive logarithm did not cross the retained branch cut')
    low_pair, outside_pair = m.fourier_pair(4, 0, 0), m.fourier_pair(4, 1, -1)
    smooth_leakage = s.expand((outside_pair.conjugate().T*v*low_pair)[0])
    need(smooth_leakage == F(1, 4), 'contact smooth-sector leakage witness failed')
    wrong = m.graded_matrix(geometry.shift(2), grading=False)
    wrong_commutator = m.cleaned(h*wrong-wrong*h)
    open_h = m.clifford_matrix(w.total(w.energies(FieldGeometry((1, 1, 4)))), 8)
    open_commutator = m.cleaned(open_h*shift-shift*open_h)
    need(not m.is_zero(wrong_commutator) and not m.is_zero(open_commutator),
         'a boundary or fermionic-grading negative control failed')
    word = m.adjacent_word(geometry.shift(2))
    need(m.word_columns(8, word) == m.graded_columns(geometry.shift(2)),
         'complete adjacent FSWAP word differs from graded translation')
    # Equality of binary-rational matrices, not admission of an eigensolver residual.
    existing = PeriodicSpectrum()
    need(np.array_equal(np.asarray(h, dtype=complex), existing.hamiltonian) and
         np.array_equal(np.asarray(shift, dtype=complex), existing.shift),
         'ring action or translation differs from the retained periodic source')
    return {'sites': 4, 'CAR_modes': 8, 'full_Fock_dimension': 256,
            'original_periodic_source_matrix_equal': True,
            'exact_spectral_ranks': [int(s.trace(projector)) for projector in projectors],
            'principal_momentum_labels_in_pi_units': ['0', '1/2', '1', '-1/2'],
            'minus_one_eigenvalue_uses_plus_pi': True,
            'alternate_minus_pi_log_is_also_conserved': True,
            'winding_identity': 'dGamma(log(t)/i) = P + 2*pi*W',
            'exact_winding_spectrum_multiplicities': {str(k): histogram[k] for k in sorted(histogram)},
            'winding_minimal_polynomial': 'W^3-W=0',
            'contact_pair_matrix_element': str(vertex),
            'zero_momentum_pair_to_outside_smooth_sector_matrix_element': str(smooth_leakage),
            'additive_log_commutator_matrix_element_in_pi_units': str(force),
            'Umklapp_pairs': {'initial': [1, 1], 'final': [-1, -1],
                             'conserved_P_in_pi_units': '1', 'additive_Q_in_pi_units': ['1', '-1']},
            'ungraded_translation_commutator_HS_squared': str(s.trace(wrong_commutator.conjugate().T*wrong_commutator)),
            'open_boundary_translation_commutator_HS_squared': str(s.trace(open_commutator.conjugate().T*open_commutator)),
            'ideal_adjacent_FSWAP_count': len(word), 'all_Fock_columns_checked': 256,
            'ideal_word_is_not_literal_pulse_error_certificate': True,
            'source_operator_sha256': hashlib.sha256(canonical({name: encoded(value) for name, value in
                [('H', h), ('T', shift), ('P_over_pi', p), ('Q_over_pi', additive), ('W', winding)]}).encode()).hexdigest()}


def run():
    prime_checks = [m.prime_branch_certificate(n) for n in (3, 5, 7, 11, 17, 31, 97)]
    return {'schema': 'oph.periodic_momentum.v1', 'issues': [1145], 'gates': ['CP-09'], 'all_pass': True,
            'domain': {'action': 'same local Weyl hopping and onsite contact, declared periodic boundary',
                      'speed': '1/8', 'interaction': '1', 'all_number_sectors_and_untouched_references': True,
                      'principal_branch': 'Arg in (-pi,pi], with +pi at -1',
                      'physical_momentum_or_duration_identified': False},
            'three_dimensional_exact': torus_certificate(), 'complete_ring_exact': ring_certificate(),
            'controlled_smooth_sector': {'maximum_number': 2, 'individual_momentum_cutoff': '1/4',
                'model_duration': '1/100', 'interaction': '1',
                'bounds': {k: str(v) for k, v in m.smooth_sector_budget(2, F(1, 4), F(1, 100)).items()},
                'general_initial_bound': 'Nmax*kappa^3/6, provided Nmax*kappa<pi',
                'general_leakage_bound': 'min(1, abs(t)*U*floor(Nmax/2))',
                'general_evolved_bound': 'Nmax*kappa^3/6 + (pi+Nmax)*leakage_bound',
                'uniform_in_volume_at_fixed_Nmax': True, 'interacting_sector_is_invariant': False,
                'not_a_finite_density_or_continuum_limit': True},
            'all_branch_locality': {'odd_prime_period_checks': [int(row['period']) for row in prime_checks],
                'spin_block_dimension': 2, 'integer_coefficient': '2*m+L*n_m, n_m arbitrary integers',
                'constant_remainder_mod_L': [int(row['constant_remainder_mod_period']) for row in prime_checks],
                'prime_root_minimal_polynomial': '1+z+...+z^(L-1)',
                'no_number_preserving_Hermitian_log_has_finite_range_on_a_longer_prime_ring': True,
                'higher_body_terms_in_a_local_density_are_included': True,
                'principal_kernel_absolute_value': 'pi/(L*sin(pi*r/L)), 1<=r<=floor(L/2)',
                'principal_operator_error_lower': '1/(R+1), for integer R>=0 and L>2*(R+1)',
                'approximant_domain': 'any sum of operators supported on site sets of torus diameter<=R, full-Fock operator norm',
                'arbitrary_branch_quantitative_lower_claimed': False},
            'scope': {'closed': 'three conserved global periodic lattice momenta, exact contact winding, and all-branch finite-range obstruction on odd-prime rings',
                'remaining': ['Local conserved momentum/stress densities or a controlled selected-state continuum replacement.',
                              'A common physical time/chart, symmetric Lorentz stress tensor and its scaling law.',
                              'Geometric modular weights, entropy/area stationarity, Newton normalization and an Einstein-register inhabitant.'],
                'algebraic_tails_and_selected_low_energy_continuum_limits_excluded': False,
                'no_local_Einstein_momentum_column_accepted': True},
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}}


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
    print(json.dumps({k: result[k] for k in ('all_pass', 'three_dimensional_exact', 'complete_ring_exact', 'all_branch_locality', 'scope')}, indent=2))


if __name__ == '__main__':
    main()
