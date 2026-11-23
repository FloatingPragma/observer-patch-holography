"""Cyclic label transport and the charged thresholds of a verified three-port guide."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.bands.measure import (
    DEGREE_FIVE_COST, HOP, bloch, digest, verify_route,
)


def quotient(x):
    x = np.asarray(x)
    return float(x.sum()/np.sqrt(x).sum()**2)


def twisted_spectrum(n, holonomy, character):
    shift = (holonomy*character) % 3
    momenta = 2*np.pi*(np.arange(n) + shift/3)/n
    energies = np.array([bloch(k)[1] for k in momenta])
    i = np.argmin(energies[:, 0])
    momentum = float(np.angle(np.exp(1j*momenta[i])))
    return energies, momentum


def assemble(shifts):
    """Explicit position-by-label contacts in the three-state basis, for a signed end."""
    n_ports = len(shifts)
    index = np.arange(3*n_ports).reshape(n_ports, 3)
    h = np.diag(np.repeat(2*HOP + DEGREE_FIVE_COST*(np.arange(n_ports)%3 == 1), 3))
    for p, shift in enumerate(shifts):
        dst = index[(p+1)%n_ports, (np.arange(3)+shift)%3]
        h[dst, index[p]] = -HOP
        h[index[p], dst] = -HOP
    return h


def direct_controls(n, holonomy, spectra, rng):
    length = 3*n
    arbitrary = rng.integers(0, 3, size=length)
    arbitrary[-1] = (holonomy - arbitrary[:-1].sum()) % 3
    boundary = np.zeros(length, dtype=int)
    boundary[-1] = holonomy
    # Canonical label at port p is the original label less the accumulated shift to p.
    accumulated = np.concatenate(([0], np.cumsum(arbitrary[:-1])%3))
    mapping = (3*np.arange(length)[:, None] + (np.arange(3)-accumulated[:, None])%3).ravel()
    h, h_boundary = assemble(arbitrary), assemble(boundary)
    np.testing.assert_array_equal(h, h_boundary[np.ix_(mapping, mapping)])
    roll = (3*np.arange(length)[:, None] + (np.arange(3)+1)%3).ravel()
    np.testing.assert_array_equal(h, h[np.ix_(roll, roll)])
    # Character convention chosen so a forward R^a contributes exp(-2πija/3).
    f = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    blocks = np.einsum('aj,parb,bk->pjrk', f.conj(), h.reshape(length, 3, length, 3), f, optimize=True)
    off_character = max(float(np.max(abs(blocks[:, j, :, k]))) for j in range(3) for k in range(3) if j != k)
    spectral_error = max(float(np.max(abs(np.linalg.eigvalsh(blocks[:, j, :, j]) - np.sort(spectra[j].ravel())))) for j in range(3))
    assert off_character < 1e-12 and spectral_error < 1e-11
    return {'all_seam_placements_gauge_proof_checked_on_nonuniform_representative': True,
            'integer_global_cyclic_shift_symmetry': True, 'real_space_dimension': len(h),
            'max_off_character_matrix_element': off_character, 'max_whole_character_spectrum_error': spectral_error}


def additive_control(e0):
    # Unit a and arbitrary interior phase choose a nonsingular, inequivalent example with no empirical input.
    phase = np.pi/24
    roots = 1 + np.sqrt(2)*np.cos(phase+2*np.pi*np.arange(3)/3)
    powers = roots**2
    np.testing.assert_allclose(quotient(powers), 2/3, atol=1e-14)
    scales = [0.01, 0.1, 1., 10., 100., 1e4]
    qs = [quotient(e0+s*powers) for s in scales]
    assert all(a < b for a, b in zip([1/3]+qs, qs+[2/3]))
    # The historical 432-state single-load compression has a=1/(3 sqrt(3)).
    # Its exclusive-read phase is 1/12 + 1/6 - 2/72. Coupling has no keV calibration.
    counted_phase = 1/12 + 1/6 - 2/72
    counted_roots = (1 + np.sqrt(2)*np.cos(counted_phase + 2*np.pi*np.arange(3)/3))/(3*np.sqrt(3))
    assert np.all(counted_roots > 0)
    counted = np.sort(counted_roots**2)
    np.testing.assert_allclose(counted.sum(), 2/9, atol=1e-14)
    np.testing.assert_allclose(quotient(counted), 2/3, atol=1e-14)
    combined = e0 + counted
    return {'register_phase_example_radians': phase, 'register_a_example': 1.,
            'register_effect_eigenvalues': powers.tolist(), 'field_charged_floor': e0,
            'scales': scales, 'Q_of_total_thresholds': qs, 'Q_of_isolated_balanced_register': quotient(powers),
            'historical_single_load_register': {
                'read_phase': counted_phase, 'isolated_effects_in_increasing_order': counted.tolist(),
                'isolated_ratios': (counted/counted[0]).tolist(), 'isolated_Q': quotient(counted),
                'illustrative_unit_coupling_total_ratios': (combined/combined[0]).tolist(),
                'illustrative_unit_coupling_total_Q': quotient(combined)},
            'small_scale_limit': 1/3, 'large_scale_limit': 2/3, 'inter_character_transitions': 0}


def measure():
    physical_n, geometry = verify_route()
    rng = np.random.default_rng(49173)
    rows = []
    for n in sorted({12, 24, physical_n, 96, 192, 384}):
        cases = []
        for holonomy in range(3):
            data = [twisted_spectrum(n, holonomy, j) for j in range(3)]
            floors = [float(e.min()) for e, _ in data]
            np.testing.assert_allclose(floors[1], floors[2], atol=1e-13, rtol=0)
            row = {'total_seam_shift_modulo_three': holonomy, 'three_character_thresholds': floors,
                   'threshold_ratios': [x/floors[0] for x in floors], 'threshold_Q': quotient(floors),
                   'minimum_momenta': [k for _, k in data], 'splitting_from_trivial_character': floors[1]-floors[0]}
            if n == physical_n:
                row['real_space_controls'] = direct_controls(n, holonomy, [e for e, _ in data], rng)
            cases.append(row)
        rows.append({'cells': n, 'verified_physical_route': n == physical_n, 'total_shifts': cases})
    e0 = float(bloch(0)[1][0]); curvature = float(bloch(0)[-1][0])
    return {'geometry': geometry, 'model_status': 'comparison with an explicitly supplied C3 label',
            'sha256': {str(p.relative_to(HERE.parent)): digest(p) for p in (HERE/'measure.py', HERE/'DESIGN.md', HERE.parent/'bands'/'measure.py', HERE.parent/'bands'/'route.json')},
            'hop': HOP, 'degree_five_energy': DEGREE_FIVE_COST, 'neutral_reference': 0, 'unshifted_guide_threshold': e0,
            'low_band_curvature_at_zero': curvature, 'asymptotic_n_squared_nontrivial_splitting': curvature/2*(2*np.pi/3)**2,
            'periodic_guides': rows, 'separate_additive_register': additive_control(e0)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = measure()
    with args.output.open('x') as out:
        json.dump(report, out, allow_nan=False, indent=2)
        out.write('\n')
