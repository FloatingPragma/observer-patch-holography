"""Rational response intervals from the complete ordered-spin cube witness."""
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from model.native_ordered_ground import certificate, source
from model.centered_particle_spectrum.certificate import decimal_bounds

ROOT = Path(__file__).resolve().parents[2]
WITNESS = 'model/native_ordered_ground/witness.npz'


def finite_response():
    with np.load(ROOT/WITNESS, allow_pickle=False) as archive:
        witness = {k: archive[k] for k in archive.files}
    verified = certificate.verify(witness)
    vectors = witness['vectors'].astype(object)
    ground = vectors[:, 0]
    norm = sum(int(x)**2 for x in ground)
    scale = 1 << certificate.BITS
    energies = [F(int(x), scale) for x in witness['energies']]
    origin = energies[0]
    e0lo, e0hi = map(F, verified['ground_energy_interval'])
    delta = F(verified['complete_polar_Weyl_radius'])+max(origin-e0lo, e0hi-origin)
    gram = F(verified['Gram_defect'])
    state_error = F(verified['normalized_ground_vector_error'])
    gap_lower = F(verified['finite_field_gap_interval'][0])
    if gap_lower <= 0:
        raise ArithmeticError('finite source ground is not isolated')
    odd = set(source.cube()[1])
    probes = {
        'local_Z0': ([1-2*(mask & 1) for mask in range(256)], 1),
        'uniform_Z_over_8': ([sum(1-2*((mask >> x) & 1) for x in range(8))
                                for mask in range(256)], 8),
        'staggered_Z_over_8': ([sum((-1 if x in odd else 1)*(1-2*((mask >> x) & 1))
                                     for x in range(8)) for mask in range(256)], 8),
    }
    h, _ = source.matrix()
    flip = np.eye(256, dtype=np.int64)[np.arange(256) ^ 255]
    if not np.array_equal(h@flip, flip@h):
        raise ArithmeticError('ground transverse-parity symmetry changed')
    rows = []
    for name, (diagonal, divisor) in probes.items():
        if any(diagonal[x^255] != -diagonal[x] for x in range(256)):
            raise ArithmeticError('probe has no exact vanishing ground mean')
        coefficients = vectors.T @ (ground*np.array(diagonal, dtype=object))
        weights = [F(int(x)**2, scale**2*norm*divisor**2) for x in coefficients]
        for eta in (F(1, 4), F(1, 16), F(1, 64)):
            shifted = [e-origin+eta for e in energies]
            chi = sum((w/e for w, e in zip(weights, shifted)), F())
            inverse2 = sum((w/(e*e) for w, e in zip(weights, shifted)), F())
            numerator = chi-eta*inverse2
            # Resolvent identities compare the exact polar spectral operator
            # with H-E0; Gram and input-vector errors are bounded separately.
            shape = gram*(2+gram)+2*state_error
            chi_error = delta/(eta*eta)+shape/eta
            second_error = 2*delta/(eta**3)+shape/(eta*eta)
            energy_error = 3*delta/(eta*eta)+shape/eta
            clo, chi_hi = max(F(), chi-chi_error), chi+chi_error
            nlo, nhi = max(F(), numerator-energy_error), numerator+energy_error
            plo, phi = max(F(), inverse2-second_error), inverse2+second_error
            if clo <= 0 or nlo <= 0 or plo <= 0:
                raise ArithmeticError('no certified positive relaxed response')
            # The exact odd probe has no ground component. Every positive
            # excitation is at least the independently certified source gap.
            # One more exact resolvent term removes the linear eta bias:
            # chi = chi_eta + eta*psi_eta + eta^2 integral dnu/[E(E+eta)^2].
            exact_chi_lo = clo+eta*plo
            exact_chi_hi = chi_hi+eta*phi+eta**2*phi/gap_lower
            rows.append({'probe': name, 'eta': str(eta),
                'regularized_inverse_interval': decimal_bounds(clo, chi_hi, 10),
                'unregularized_inverse_interval': decimal_bounds(exact_chi_lo, exact_chi_hi, 10),
                'minimum_unit_tangent_energy_interval': decimal_bounds(1/(4*exact_chi_hi), 1/(4*exact_chi_lo), 10),
                'feasible_regularized_tangent_energy_interval': decimal_bounds(nlo/(4*chi_hi**2), nhi/(4*clo**2), 10),
                'resolvent_filter_success_interval': decimal_bounds(eta**2*plo, eta**2*phi, 10),
                'all_complementary_reject_probability_interval': decimal_bounds(1-eta**2*phi, 1-eta**2*plo, 10),
                'probe_LCU_normalization': '1',
                'ground_mean_zero_from_exact_parity': True})
    return {'source': 'unchanged open eight-site cube, J=1, B=1/4, unique H_B ground',
            'complete_source_dimension': 256,
            'complete_integer_witness_verified': True,
            'finite_field_ground_gap_lower': str(gap_lower),
            'rows': rows,
            'finite_ground_is_the_zero_field_infinite_ordered_state': False,
            'large_native_resolvent_tape_executed': False,
            'all_success_and_complementary_rejection_weights_retained': True,
            'physical_stiffness_speed_or_particle_mass_inferred': False}
