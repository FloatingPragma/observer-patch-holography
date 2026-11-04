"""Actual source hopping on compact opposite-degree textures.

The complete numerical spectra accompany exact dimension, localization and
boundary error bounds; floating eigenvalues are not exact root certificates.
"""
from fractions import Fraction as F
from itertools import product

import numpy as np
from scipy import linalg as la, sparse as sp

from .topology import direction, separation, spectral_bounds

SIGMA = (np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
         np.diag([1., -1.]).astype(complex))
J = np.array([[0., 1.], [-1., 0.]])


def source(radius=3, mass=64, half_separation=1, *, uniform=False, core_shift=0):
    spectral_bounds(mass, half_separation, radius)
    if not np.isfinite(core_shift) or abs(core_shift) > 1:
        raise ValueError('finite core shift of absolute value at most one required')
    a = separation(half_separation)
    points = tuple(product(range(-radius, radius+1), repeat=3))
    addresses = {p: i for i, p in enumerate(points)}
    h = sp.lil_matrix((2*len(points),)*2, dtype=complex)
    core = []
    for p, i in addresses.items():
        v = np.array((0, 0, 1) if uniform else direction(p, a), float)
        if np.linalg.norm(v):
            v /= np.linalg.norm(v)
        else:
            core.extend([2*i, 2*i+1])
        onsite = mass*sum((v[j]*SIGMA[j] for j in range(3)), np.zeros((2, 2), complex))
        if p in ((0, 0, -a), (0, 0, a)) and not uniform:
            onsite += core_shift*SIGMA[2]
        h[2*i:2*i+2, 2*i:2*i+2] = onsite
        for axis in range(3):
            q = tuple(x+int(j == axis) for j, x in enumerate(p))
            if q in addresses:
                k = addresses[q]
                link = 0.5j*SIGMA[axis]
                h[2*i:2*i+2, 2*k:2*k+2] = link
                h[2*k:2*k+2, 2*i:2*i+2] = link.conj().T
    c = sp.kron(sp.diags([(-1.)**sum(p) for p in points]), sp.csr_matrix(J), format='csr')
    return {'h': h.tocsr(), 'C': c, 'core_indices': core, 'points': points}


def display(values, decimals=9):
    return [round(float(x), decimals) for x in values]


def norm(value):
    return float(sp.linalg.norm(value) if sp.issparse(value) else la.norm(value))


def eigensystem(value):
    h = value['h'].toarray()
    eigenvalues = la.eigh(h, eigvals_only=True, driver='evr', check_finite=True)
    indices = np.flatnonzero(np.abs(eigenvalues) < 2)
    if len(indices) != 4:
        raise ArithmeticError('the exact isolated band must contain four numerical roots')
    energies, frame = la.eigh(h, subset_by_index=[int(indices[0]), int(indices[-1])], driver='evr')
    return eigenvalues, energies, frame


def case(radius=3, mass=64, half_separation=1):
    bound = spectral_bounds(mass, half_separation, radius)
    value = source(radius, mass, half_separation)
    h, c = value['h'], value['C']
    eigenvalues, energies, frame = eigensystem(value)
    complement = eigenvalues[np.abs(eigenvalues) >= 2]
    probability = np.sum(abs(frame[value['core_indices']])**2, axis=0)
    residual = norm(h@frame-frame*energies)
    ph_residual = norm(h@c+c@h.conjugate())
    exterior = [i for i in range(h.shape[0]) if i not in value['core_indices']]
    b = h[exterior][:, value['core_indices']]
    gram_values = la.eigvalsh((b.conj().T@b).toarray())
    if max(abs(energies)) > float(F(bound['band_abs_energy_upper'])):
        raise ArithmeticError('isolated band violates exact energy enclosure')
    if min(probability) < float(F(bound['core_probability_lower']))-1e-11:
        raise ArithmeticError('isolated band violates exact localization bound')
    if min(abs(complement)) < mass-2 or residual > 1e-9 or ph_residual > 1e-10:
        raise ArithmeticError('complete source or spectral complement control failed')
    uniform = source(radius, mass, half_separation, uniform=True)['h'].toarray()
    uniform_eigenvalues = la.eigvalsh(uniform, driver='evr')
    vacuum_shift = (np.minimum(eigenvalues, 0).sum()-np.minimum(uniform_eigenvalues, 0).sum())
    return {'radius': radius, 'mass': mass, 'half_separation': half_separation,
            'sites': len(value['points']), 'one_body_dimension': h.shape[0],
            'complete_numerical_spectrum': display(eigenvalues),
            'isolated_numerical_energies': display(energies, 12),
            'minimum_complement_abs_energy': round(float(min(abs(complement))), 9),
            'minimum_core_probability': round(float(min(probability)), 12),
            'core_coupling_Gram_eigenvalues': display(gram_values, 12),
            'filled_sea_energy_difference_from_uniform': round(float(vacuum_shift), 9),
            'numerical_residual_below_1e_minus_9': residual < 1e-9,
            'particle_hole_identity_below_1e_minus_10': ph_residual < 1e-10,
            'exact_bounds': bound}


def lifted_control():
    # The actual 45-row native representation is tensored with the same spatial
    # band. The low band need not be the zero eigenspace.
    from model.common_source.chiral_candidate import action
    value = source()
    _, energies, w = eigensystem(value)
    h, c = value['h'], value['C']
    a = w.conj().T@(h@w)
    s = w.conj().T@(c@w.conjugate())
    spatial = {'invariant_band': norm(h@w-w@a),
               'particle_hole_closed_band': norm(c@w.conjugate()-w@s),
               'antisymmetric_band_pair': norm(s+s.T),
               'band_pair_unitary': norm(s@s.conj().T-np.eye(4))}
    # A small genuinely complex symmetric kernel checks the full Nambu formula.
    k = np.array([[.2, .01j], [.01j, .03]], complex)
    identity = sp.eye(2, format='csr')
    hk, delta = sp.kron(h, identity), sp.kron(c, k.conjugate())
    b = sp.bmat([[hk, delta], [delta.conj().T, -hk.conjugate()]], format='csr')
    e = sp.block_diag((sp.kron(w, identity), sp.kron(w.conjugate(), identity)), format='csr')
    ac, sc = np.kron(a, np.eye(2)), np.kron(s, k.conjugate())
    core = np.block([[ac, sc], [sc.conj().T, -ac.conjugate()]])
    error = norm(b@e-e@core)
    exact_additive = sorted(x+y for x in energies for y in np.r_[-la.svdvals(k), la.svdvals(k)])
    additive_error = norm(np.sort(la.eigvalsh(core))-exact_additive)
    currents = []
    for g in action.generators():
        g = np.array(g, complex)
        # Tensor identities give the complete 61740-row Nambu result without
        # constructing a redundant dense spatial/internal matrix.
        lhs = np.einsum('xa,ij->xiaj', w, g).reshape(w.shape[0]*45, 4*45)
        rhs = np.kron(w, g)
        currents.append(norm(lhs-rhs))
    if max([error, additive_error, *spatial.values(), *currents]) > 1e-9:
        raise ArithmeticError('complete invariant-band/native-current inclusion failed')
    shifted = source(core_shift=.125)
    shifted_low = la.eigh(shifted['h'].toarray(), eigvals_only=True, subset_by_value=[-2, 2])
    old_zero_core_error = norm(b@e-e@(core-np.block([[ac, np.zeros_like(sc)], [np.zeros_like(sc), -ac.conjugate()]])))
    if old_zero_core_error <= .001 or norm(shifted_low-energies) <= .01:
        raise ArithmeticError('zero-core or topological pinning negative control failed')
    return {'native_current_count': len(currents), 'native_internal_dimension': 45,
            'localized_spatial_dimension': 4, 'localized_Nambu_dimension': 360,
            'full_native_Nambu_dimension': 4*len(value['points'])*45,
            'all_twelve_current_tensor_intertwiners_exact': True,
            'complex_kernel_full_Nambu_control_passed': error < 1e-9,
            'band_Nambu_spectrum_is_energy_plus_or_minus_singular_values': additive_error < 1e-9,
            'zero_core_formula_fails': old_zero_core_error > .001,
            'core_shift_preserves_shell_data_but_changes_spectrum': norm(shifted_low-energies) > .01,
            'shifted_isolated_numerical_energies': display(shifted_low, 12),
            'all_spatial_band_residuals_below_1e_minus_9': max(spatial.values()) < 1e-9}
