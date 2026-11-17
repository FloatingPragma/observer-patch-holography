"""Independent interval generating-function trace; no producer imports."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'species_leptons/chiral_ring4_fidelity'
mp.iv.dps = 45
iv = mp.iv


def num(x):
    return iv.mpf(str(x))


def upper(x):
    return math.nextafter(float(x.b), math.inf)


def main():
    inputs = [SOURCE/'RESULT.json', SOURCE/'COMPLETE_TRACE_INPUTS.npz',
              SOURCE/'inputs/prepared_results.json', SOURCE/'inputs/pointwise_controls.json',
              SOURCE/'inputs/coulomb_results.json', SOURCE/'inputs/source.py']
    result = json.loads(inputs[0].read_text())
    saved = np.load(inputs[1])
    prepared = json.loads(inputs[2].read_text())
    points = {r['grid_index']: r for r in json.loads(inputs[3].read_text()) if r['N'] == 8}
    coulomb = json.loads(inputs[4].read_text())
    anchor = num(result['trace_anchor_energy'])
    beta = iv.mpf(5)/4
    grid = 1024
    assert result['holonomy_grid'] == grid and len(points) == grid
    # Independently build the boundary derivative only; all entries are
    # Gaussian integers of small exact magnitude in binary64 storage.
    X = np.array([[0, 1], [1, 0]], complex)
    Z = np.diag([1, -1]).astype(complex)
    C = np.zeros((16, 16), complex)
    for wall in range(2):
        C[12+2*wall:14+2*wall, 2*wall:2+2*wall] = -4*(Z+1j*X)
    J = 1j*(C-C.conj().T)
    assert np.array_equal(J, J.conj().T)
    assert np.array_equal(J@J@J, 64*J)
    assert np.trace(J) == 0 and np.trace(J@J) == 256
    # Hermiticity + J^3=64J + the two traces gives exactly two +/-8
    # eigenvalues each. N8 occupies both positive levels and six zeros.
    derivative_bound = 16
    all_rows = []
    max_band_error = 0.
    total = iv.mpf(0)
    producer_beta_index = list(saved['betas']).index(1.25)
    for j in range(grid):
        positive = []
        theta = 2*iv.pi*j/grid
        for k in range(4):
            momentum = (theta+2*iv.pi*k)/4
            w = 1-iv.cos(momentum)
            a = iv.sqrt(1+4*w*w)
            small = 2*w*w/(a+1)
            large = (a+1)/2
            for singular in (small, large):
                positive.append(8*iv.sqrt(iv.sin(momentum)**2+singular**2))
        levels = positive+[-x for x in positive]
        diagnostic = np.sort([float(x.mid) for x in levels])
        band_error = float(np.max(np.abs(diagnostic-saved['one_particle_bands'][j])))
        assert band_error < 1e-11
        max_band_error = max(max_band_error, band_error)
        coefficients = [iv.mpf(1)]+[iv.mpf(0) for _ in range(8)]
        for e in levels:
            weight = iv.exp(beta*(anchor/8-e))
            for degree in range(8, 0, -1):
                coefficients[degree] += weight*coefficients[degree-1]
        free_ground = iv.exp(beta*(anchor+sum(positive, iv.mpf(0))))
        ground_lower = num(points[j]['improved_pointwise_lower'])-iv.mpf('1e-8')
        # Min-max bounds every remaining interacting ordered eigenvalue
        # below in energy by its free counterpart, because B>=0.
        corrected = coefficients[8]-free_ground+iv.exp(beta*(anchor-ground_lower))
        assert float(corrected.a) > 0
        assert upper(corrected) <= saved['pointwise_scaled_trace_upper'][j, producer_beta_index]*(1+1e-9)
        total += corrected
        all_rows.append({'grid_index': j, 'scaled_trace_upper': upper(corrected),
                         'analytic_band_comparison_error': band_error})
        if j % 256 == 0:
            print('independent interval trace grid', j, flush=True)
    mean = total/grid
    rotor = 1+2*sum((iv.exp(-beta*n*n) for n in range(1, 17)), iv.mpf(0))
    rotor += 2*iv.exp(-beta*17**2)/(1-iv.exp(-beta*35))
    old_quadrature = iv.exp(beta*16*iv.sqrt(2)*iv.pi/grid)
    exact_quadrature = iv.exp(beta*derivative_bound*iv.pi/grid)
    counts = []
    for case in prepared:
        mass_shift = 4*abs(num(case['mass']))
        rayleigh_upper = num(case['variational_energy_upper'])
        cut = rayleigh_upper+iv.mpf(1)/8
        common = mean*rotor*iv.exp(beta*(cut-anchor+mass_shift))
        old_count, new_count = upper(common*old_quadrature), upper(common*exact_quadrature)
        source_lowers = {r['N']: num(r['complete_source_lower'])-(0 if r['N'] in (0, 16) else mass_shift)
                         for r in coulomb['every_number_sector']}
        excluded = all(float(v.a) > float(cut.b) for n, v in source_lowers.items() if n != 8)
        assert new_count <= old_count < 2 and excluded
        produced = next(r for r in result['vacuum_consumers'] if r['case'] == case['case'])
        assert old_count <= produced['count_witness']['eigenvalue_count_upper']*(1+1e-9)
        residual = num(case['full_preparation_spectral_residual'])+iv.mpf('1e-8')
        weight = min(1., upper((8*residual)**2))
        assert abs((1-weight)-produced['true_ground_fidelity_lower']) < 1e-12
        # Exact spectral inequality: weight above cut <= residual^2/(cut-Rayleigh)^2.
        # The recorded Rayleigh upper includes padding, so offset1/8 is conservative.
        energy_lower = num(case['variational_energy'])-iv.mpf('1e-6')-8*residual**2
        assert abs(float(energy_lower.a)-produced['actual_ground_energy_interval'][0]) < 1e-11
        counts.append({'case': case['case'], 'same_Lipschitz_count_upper': old_count,
                       'exact_boundary_derivative_count_upper': new_count,
                       'other_number_sectors_excluded': excluded,
                       'ground_fidelity_lower_rechecked': 1-weight})
    out = {'schema': 'oph.root.vacuum_trace_review.v1', 'all_checks_pass': True,
        'all_holonomies': grid, 'all_Fock_states_retained_by_generating_polynomial': math.comb(16, 8),
        'independent_analytic_band_max_error': max_band_error,
        'exact_boundary_current_N8_norm': derivative_bound,
        'canonical_trace_mean_upper': upper(mean), 'rotor_trace_upper': upper(rotor),
        'producer_quadrature_upper': upper(old_quadrature),
        'stronger_quadrature_upper': upper(exact_quadrature), 'vacuum_counts': counts,
        'consumed_not_rederived': ['Coulomb pointwise ground lower bounds',
            'number-sector lower bounds', 'full preparation Rayleigh and infinite-flux residual',
            'prepared charged Green residual bounds'],
        'physical_lepton_poles_solved': False,
        'source_pins': {str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (HERE/'TRACE_ROWS.json').write_text(json.dumps(all_rows, indent=2)+'\n')
    (HERE/'RESULT.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2), flush=True)


if __name__ == '__main__':
    main()
