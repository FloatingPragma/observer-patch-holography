"""Replay collision-suppressed local and coarse momentum Ward bounds."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as s

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..coarse_energy_ward import ward as w
from ..periodic_momentum import momentum as p
from . import ward as d

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_PATHS = (
    'model/dilute_momentum_ward/ward.py', 'model/dilute_momentum_ward/packet.py',
    'model/dilute_momentum_ward/test_ward.py', 'model/dilute_momentum_ward/README.md',
    'model/periodic_momentum/momentum.py', 'model/lattice_momentum_balance/momentum.py',
    'model/coarse_energy_ward/ward.py', 'model/full_fock_field/field.py',
    'model/interacting_field/lattice.py', 'Lean/QFT/DiluteMomentumWard.lean',
)


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def finite_two_particle_certificate():
    geometry = p.PeriodicGeometry((1, 1, 4))
    data = d.columns(geometry)
    potential, masks = d.occupation_matrix(data['V'], 8)
    force, _ = d.occupation_matrix(w.total(data['columns'][2]['force']), 8)
    index = {mask: j for j, mask in enumerate(masks)}
    shifted = p.graded_columns(geometry.shift(2))
    T = s.ImmutableSparseMatrix(len(masks), len(masks),
        {(index[shifted[mask][0]], j): shifted[mask][1] for j, mask in enumerate(masks)})
    _, projectors = p.spectral_log(T, 4)
    force_squared = []
    for k, projector in enumerate(projectors):
        sigma2 = (0, 1, 2, 1)[k]
        f = p.cleaned(projector*force*projector)
        need(p.is_zero(f**3-sigma2*f), 'sharp contact force polynomial failed')
        need(s.trace(p.cleaned(f*f)) == 2*sigma2, 'sharp force rank or norm failed')
        force_squared.append(sigma2)
    # Independent exterior-square projection of spin-complete momenta {0,1}.
    orbitals = s.ImmutableSparseMatrix(8, 4,
        {(2*x+spin, 2*k+spin): s.I**((-k*x) % 4)/2
         for x in range(4) for k in (0, 1) for spin in (0, 1)})
    one = p.cleaned(orbitals*orbitals.conjugate().T)
    pairs = [tuple(j for j in range(8) if mask >> j & 1) for mask in masks]
    omega = p.cleaned(s.ImmutableSparseMatrix(len(masks), len(masks),
        {(i, j): one[a, c]*one[b, e]-one[a, e]*one[b, c]
         for i, (a, b) in enumerate(pairs) for j, (c, e) in enumerate(pairs)}))
    need(p.is_zero(omega*omega-omega) and s.trace(omega) == 6, 'band exterior projection failed')
    values = (F(1, 4), F(1, 2), F(1, 4), F())
    for projector, value in zip(projectors, values, strict=True):
        compressed = p.cleaned(projector*omega*potential*omega*projector)
        need(p.is_zero(compressed*compressed-value*compressed) and
             s.trace(compressed) == value, 'contact principal angle differs from convolution count')
        escape = p.cleaned((s.eye(len(masks))-omega)*potential*omega*projector)
        need(s.trace(p.cleaned(escape.conjugate().T*escape)) == value*(1-value),
             'initial contact escape norm lost its outside-band factor')
    return {'sites': 4, 'full_N2_dimension': len(masks),
            'complete_four_total_momentum_sectors': True,
            'force_norm_squared_in_U_squared_units': force_squared,
            'force_norm_formula': 'sqrt(2)*U*abs(sin(K_j/2))',
            'nonzero_force_rank_per_nonzero_sector': 2,
            'band_momenta': [0, 1], 'band_N2_dimension': 6,
            'contact_principal_angles_squared': [str(x) for x in values],
            'escape_squared_per_K': [str(x*(1-x)) for x in values]}


def cell_diagnostic():
    geometry = p.PeriodicGeometry((1, 1, 8))
    data = d.columns(geometry)
    h, masks = d.occupation_matrix(data['H'], 16)
    potential, _ = d.occupation_matrix(data['V'], 16)
    initial = d.zero_pair(masks, 8)
    need((initial.conjugate().T*initial)[0] == 1, 'zero pair not normalized')
    projected = potential*initial
    need((projected.conjugate().T*projected)[0] == F(1, 8), 'exact band collision norm is incorrect')
    outside = projected-initial*(initial.conjugate().T*projected)[0]
    need((outside.conjugate().T*outside)[0] == F(7, 64), 'exact initial escape norm is incorrect')
    h = np.asarray(h, dtype=complex)
    values, vectors = np.linalg.eigh(h)
    differences = values[:, None]-values[None, :]
    time = .5
    psi = np.asarray(initial, dtype=complex)
    def matrix(operator):
        return np.asarray(d.occupation_matrix(operator, 16)[0], dtype=complex)
    def evolve(operator):
        return vectors@(np.exp(1j*time*differences)*(vectors.conj().T@operator@vectors))@vectors.conj().T
    def integral(operator):
        factor = time*np.exp(1j*time*differences/2)*np.sinc(time*differences/(2*np.pi))
        return vectors@(factor*(vectors.conj().T@operator@vectors))@vectors.conj().T
    column = data['columns'][2]
    divergence = w.divergence(column['stress'], 8)
    identity_residual, unforced = 0., 0.
    for local, flux, force in zip(column['local'], divergence, column['force'], strict=True):
        local = matrix(local)
        residual = evolve(local)-local+integral(matrix(flux))
        identity_residual = max(identity_residual, float(np.linalg.norm((residual-integral(matrix(force)))@psi)))
        unforced = max(unforced, float(np.linalg.norm(residual@psi)))
    bound = d.budget(F(1, 8), 2, F(1, 2))['local_Ward_cell_defect_upper']
    need(identity_residual < 1e-12 and 1e-5 < unforced < float(bound), 'finite cell identity or nonzero-force control failed')
    return {'sites': 8, 'N2_dimension': len(masks), 'initial_band': 'both particles at momentum zero',
            'exact_initial_collision_norm_squared': '1/8', 'exact_initial_escape_norm_squared': '7/64',
            'duration': time, 'complete_forced_cell_identity_residual': identity_residual,
            'unforced_cell_defect_norm': unforced, 'rigorous_local_cell_bound': str(bound),
            'diagnostic_tolerance': 1e-12, 'floating_diagnostic_is_not_the_proof': True}


def run():
    geometry = p.PeriodicGeometry()
    data = d.columns(geometry)
    band = d.band_counts()
    bounds = d.budget(band['fraction'], 2, F(1, 8), cutoff=F(11, 112), block_sites=8)
    payload = [[[str(mask), str(a), str(b)] for mask, (a, b) in sorted(force.terms.items())]
               for column in data['columns'] for force in column['force']]
    return {'schema': 'oph.dilute_momentum_ward.v1', 'issues': [1145], 'gates': ['CP-09'], 'all_pass': True,
            'domain': {'action': 'unchanged periodic complete-Fock Weyl/contact action',
                      'speed': '1/8', 'interaction': '1', 'initial_state_domain': 'spin-complete band B and N<=Nmax',
                      'band_projector_is_not_applied_during_evolution': True,
                      'arbitrary_untouched_references': True, 'physical_time_or_momentum_identified': False},
            'exact_local_attachment': {'sides': [3, 3, 3], 'CAR_modes': 54, 'full_Fock_dimension': 2**54,
                'local_Ward_identities': 81, 'geometric_oriented_faces': len(data['oriented_edges']),
                'nonzero_stress_faces_per_column': [len(c['stress']) for c in data['columns']],
                'rectangular_coarse_partitions_per_column': 27,
                'source_force_sha256': hashlib.sha256(canonical(payload).encode()).hexdigest()},
            'analytic_theorem': {'projected_contact_bound': 'Omega D Omega <= f*floor(Nmax/2)',
                'initial_interaction_norm_upper': 'U*floor(Nmax/2)*sqrt(f)',
                'local_force_initial_norm_upper': 'min(U,4*U*f)',
                'local_cell_defect_upper': 'min(U*T,4*U*f*T+U^2*floor(Nmax/2)*sqrt(f)*T^2/2)',
                'block_bound': 'sum of constituent site bounds, internal kinetic faces cancel',
                'two_particle_momentum_error_upper': 'kappa^3/3+4*kappa*min(1,T*U*sqrt(f)), 2*kappa<pi',
                'zero_or_one_particle_interaction_force': 'exactly zero'},
            'finite_band_example': {'sides': [64, 64, 64], 'band_side_labels': [-1, 0, 1],
                'spatial_volume': band['volume'], 'band_modes_without_spin': len(band['band']),
                'band_fraction': str(band['fraction']), 'nonzero_contact_K_sectors': len(band['contact_eigenvalues']),
                'maximum_contact_angle_squared': str(band['contact_norm_squared']),
                'exact_max_initial_escape_squared': str(band['initial_escape_norm_squared']),
                'cutoff_upper': '11/112 > pi/32', 'Nmax': 2, 'duration': '1/8', 'block_sites': 8,
                'bounds': {k: str(v) for k, v in bounds.items()},
                'analytic_resource_example_not_large_state_execution': True},
            'exact_two_particle_force_and_contact': finite_two_particle_certificate(),
            'finite_cell_diagnostic': cell_diagnostic(),
            'scope': {'closed': 'state-dependent local momentum Ward defects with exact coarse aggregation and collision-suppressed leakage',
                'fixed_U_fixed_time_dilute_long_wavelength_error_tends_to_zero': True,
                'critical_band_edge_scaling_covered': False,
                'critical_scaling_boundary': 'f=O(a^3), U=O(a^-2) makes U*sqrt(f)=O(a^-1/2); this bound does not control that limit',
                'remaining': ['finite-density or critically scaled interacting local stress limit',
                              'symmetric Lorentz stress, geometric modular weights and Einstein-register inhabitant',
                              'preparation/readout and physical time/unit identification']},
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
    print(json.dumps({k: result[k] for k in ('all_pass', 'exact_local_attachment', 'finite_band_example', 'scope')}, indent=2))


if __name__ == '__main__':
    main()
