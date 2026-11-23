"""Clock the written signed guide with the historical register's positive event effect."""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.bands.measure import (
    bloch, digest, evolve, packet, paired, pair_space, verify_route,
)


def register():
    """Rebuild the historical compression from its literal integer coefficient matrices."""
    phase = Fraction(1, 12) + Fraction(1, 6) - 2*Fraction(1, 72)
    shift = np.roll(np.eye(3), 1, axis=0)
    z = np.exp(1j*float(phase))
    amplitude = (32*np.eye(3) + 16*np.sqrt(2)*(z*shift + z.conjugate()*shift.T))/(96*np.sqrt(3))
    effect = amplitude.conj().T @ amplitude
    f = np.exp(-2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    roots = np.sqrt(2)/3*(1/np.sqrt(6) + np.cos(float(phase)+2*np.pi*np.arange(3)/3)/np.sqrt(3))
    np.testing.assert_allclose(amplitude@f, f*roots, atol=1e-15, rtol=0)
    np.testing.assert_allclose(effect@f, f*roots**2, atol=1e-15, rtol=0)
    assert np.all(roots > 0)
    return f, effect, roots**2, {'underlying_register_dimension': 432, 'compressed_effect_dimension': 3,
                                  'literal_amplitude_coefficients': {'identity': 32, 'each_oriented_shift_times_sqrt2': 16,
                                                                    'denominator_times_sqrt3': 96},
                                  'repair_plus_probability': '1/12', 'hold_probability': '1/6',
                                  'joint_probability': '1/72', 'exclusive_clock_angle_radians': str(phase),
                                  'amplitude_on_character_0_1_2': roots.tolist(), 'effect_on_character_0_1_2': (roots**2).tolist()}


def signed_release(n, ks, spatial, effect_weights):
    a, b, h, exchange = pair_space(n)
    vectors = np.array([s[2] for s in spatial])
    velocities = np.array([s[3] for s in spatial])
    home, away = n//4, 3*n//4
    projectile, speed = packet(ks, vectors, velocities, 0, np.pi/2, home, True)
    spectator = packet(ks, vectors, velocities, 0, 0., away, True)[0]
    initial = paired(a, b, projectile, spectator)
    offsets = np.stack([((port//3-center+n//2)%n)-n//2 for port, center in [(a,home),(b,away)]], axis=1)
    swapped_offsets = np.stack([((port//3-center+n//2)%n)-n//2 for port, center in [(b,home),(a,away)]], axis=1)
    results = []
    time = 4.
    # Compare a four-cell flight at each character's elapsed laboratory time.
    # The exact product law makes one spatial propagation sufficient for all three.
    reference_flight = evolve(h, initial, time, 2)
    reference_probability = abs(reference_flight)**2
    reference_displacement = (reference_probability[-1]-reference_probability[0])@offsets
    np.testing.assert_allclose(reference_displacement/time, [speed, 0.], atol=1e-6, rtol=0)
    assert reference_displacement[0] > 4.
    for j, w in enumerate(effect_weights):
        waves = evolve(h*w, initial, time, 2)
        probability = abs(waves)**2
        np.testing.assert_allclose(probability.sum(axis=1), 1., atol=1e-12, rtol=0)
        displacement = (probability[-1]-probability[0])@offsets
        opposite = (probability[-1,exchange]-probability[0,exchange])@swapped_offsets
        np.testing.assert_allclose(displacement, opposite, atol=1e-12, rtol=0)
        np.testing.assert_allclose(displacement/time, [w*speed, 0.], atol=1e-8, rtol=0)
        results.append({'character': j, 'shared_elapsed_time': time, 'effective_spatial_time': w*time,
                        'positive_projectile_displacement_cells': displacement[0], 'negative_projectile_displacement_cells_swapped': opposite[0],
                        'spectator_displacement_cells': displacement[1], 'measured_mean_projectile_speed': displacement[0]/time,
                        'predicted_wavepacket_mean_speed': speed*w, 'total_probability': probability[-1].sum(),
                        'time_for_four_native_spatial_units': time/w,
                        'projectile_displacement_cells_at_that_time': reference_displacement[0]})
    return {'dimensions_signed_pair': len(a), 'compensator_initial_separation_cells': n//2,
            'source_band': 1, 'central_momentum_over_pi': .5, 'logical_clock_shared_by_this_pair_fixture': True, 'reads': results}


def measure():
    n, route = verify_route()
    chars, effect, weights, raw = register()
    ks = 2*np.pi*np.arange(n)/n
    spatial = [bloch(k) for k in ks]
    dk = 1e-3
    nine_error, dark_error, slope_error, curvature_error = 0., 0., 0., 0.
    for h, energies, _, _, _ in spatial:
        actual = np.kron(h, effect)
        expected = np.sort(np.outer(energies, weights).ravel())
        nine_error = max(nine_error, float(np.max(abs(np.linalg.eigvalsh(actual)-expected))))
        for j in range(3):
            test = np.kron(np.eye(3), chars[:, [j]])
            dark_error = max(dark_error, float(np.linalg.norm(actual@test-test@(weights[j]*h), ord=2)))
    assert nine_error < 1e-12 and dark_error < 1e-13
    rows = []
    for j in np.argsort(weights):
        w = weights[j]
        internal = np.kron(np.eye(3), chars[:, [j]])
        def direct(k):
            full = np.kron(bloch(k)[0], effect)
            return np.linalg.eigvalsh(internal.conj().T@full@internal)
        for b in range(3):
            energies = [w*x[1][b] for x in spatial]
            lo, hi = int(np.argmin(energies)), int(np.argmax(energies))
            floor, curve = energies[lo], w*spatial[lo][4][b]
            forward = w*spatial[n//4][3][b]
            slope_numeric = (direct(np.pi/2+dk)[b] - direct(np.pi/2-dk)[b])/(2*dk)
            curve_numeric = (direct(ks[lo]+dk)[b]-2*direct(ks[lo])[b]+direct(ks[lo]-dk)[b])/dk**2
            slope_error = max(slope_error, abs(slope_numeric-forward)/w)
            curvature_error = max(curvature_error, abs(curve_numeric-curve)/w)
            rows.append({'character': int(j), 'orbital_band': b+1, 'threshold': floor, 'threshold_momentum_over_pi': ks[lo]/np.pi,
                         'maximum': energies[hi], 'energy_at_zero_momentum': energies[0], 'threshold_curvature': curve,
                         'threshold_times_curvature': floor*curve, 'speed_at_pi_over_two': forward,
                         'energy_at_each_of_48_momenta': energies})
    assert slope_error < 1e-5 and curvature_error < 1e-4
    lowest = [r for r in rows if r['orbital_band'] == 1]
    summary = {}
    for key in ('threshold', 'threshold_curvature', 'threshold_times_curvature', 'speed_at_pi_over_two'):
        v = np.array([r[key] for r in lowest])
        summary[key+'_ratios'] = (v/v[0]).tolist()
    rest = np.array([r['threshold'] for r in lowest])
    summary['Q'] = float(rest.sum()/np.sqrt(rest).sum()**2)
    return {'model_status': 'comparison Hamiltonian obtained by promoting an auxiliary event effect to a spatial clock',
            'route': route, 'neutral_reference_energy': 0., 'physical_scale': None,
            'sha256': {str(p.relative_to(HERE.parent)): digest(p) for p in [HERE/'measure.py', HERE/'DESIGN.md', HERE.parent/'bands'/'measure.py', HERE.parent/'bands'/'route.json']},
            'historical_register_compression': raw, 'momenta_radians_per_cell': ks.tolist(), 'nine_bands': rows, 'homologous_lowest_orbital': summary,
            'max_full_nine_eigenvalue_error': nine_error, 'max_invariant_character_residual': dark_error,
            'max_scaled_finite_difference_velocity_error': slope_error, 'max_scaled_finite_difference_curvature_error': curvature_error,
            'finite_signed_release': signed_release(n, ks, spatial, weights)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    report = measure()
    with args.output.open('x') as out:
        json.dump(report, out, allow_nan=False, indent=2)
        out.write('\n')
