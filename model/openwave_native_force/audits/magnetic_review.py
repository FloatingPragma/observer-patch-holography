"""Independent link, spectrum and matrix-function audit of magnetic response."""
from itertools import product
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
import sympy as s

from model.openwave_magnetic_response import response as producer

ROOT = Path(__file__).resolve().parents[3]


def internal_matrices():
    """Rebuild the folded source from its bit labels, without its constructor."""
    paulis = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]),
              s.diag(1, -1))
    labels = list(product((0, 1), repeat=3))
    gamma = [s.diag(*[-(-1)**b[j]*paulis[j] for b in labels]) for j in range(3)]
    beta = s.zeros(16)
    for i, b in enumerate(labels):
        k = labels.index(tuple(1-v for v in b))
        beta[2*i:2*i+2, 2*k:2*k+2] = s.eye(2)
    spin = -s.I*gamma[0]*gamma[1]
    sectors = {}
    for energy, orientation in product((-1, 1), repeat=2):
        p = (s.eye(16)+energy*beta)*(s.eye(16)+orientation*spin)/4
        assert p*p == p and p.H == p and p.rank() == 4
        sectors[str((energy, orientation))] = p.rank()
    assert all(beta*g+g*beta == s.zeros(16) for g in gamma)
    return [np.asarray(g, complex) for g in gamma], np.asarray(beta, complex), np.asarray(spin, complex), sectors


def link_momenta(n, spacing, field, charge):
    """Stamp every directed bond independently in the declared open box."""
    matrices = [sparse.lil_matrix((n*n, n*n), dtype=complex) for _ in range(2)]
    for x, y in product(range(n), repeat=2):
        start = x*n+y
        if x+1 < n:
            target = (x+1)*n+y
            matrices[0][start, target] = 1j/(2*spacing)
            matrices[0][target, start] = -1j/(2*spacing)
        if y+1 < n:
            target = start+1
            phase = np.exp(1j*charge*spacing*field*(x-n//2)*spacing)
            matrices[1][start, target] = 1j*phase/(2*spacing)
            matrices[1][target, start] = -1j*phase.conjugate()/(2*spacing)
    return [p.tocsr() for p in matrices]


def sqrt_remainder():
    worst = mp.mpf(0)
    with mp.workdps(90):
        coefficients = [mp.binomial(mp.mpf('0.5'), k) for k in range(7)]
        assert coefficients[1:] == [mp.mpf(1)/2, -mp.mpf(1)/8, mp.mpf(1)/16,
                                     -mp.mpf(5)/128, mp.mpf(7)/256, -mp.mpf(21)/1024]
        for i in range(1, 1001):
            x = mp.mpf(i)/1000
            delta = mp.sqrt(1+x)-sum(c*x**k for k, c in enumerate(coefficients))
            bound = mp.mpf(33)/2048*x**7
            assert 0 < delta <= bound
            worst = max(worst, delta/bound)
        # Omitting the sixth-order term does not retain a seventh-order bound.
        x = mp.mpf(1)/10
        mutated = abs(mp.sqrt(1+x)-sum(c*x**k for k, c in enumerate(coefficients[:6])))
        assert mutated > mp.mpf(33)/2048*x**7
    return {'scalar_interval_samples': 1000, 'worst_error_over_bound': float(worst),
            'sixth_order_term_deletion_detected': True}


def run():
    paths = list((ROOT/'model/openwave_magnetic_response').glob('*.py'))
    paths += list((ROOT/'model/openwave_magnetic_response').glob('*.md'))
    paths += [ROOT/'model/openwave_magnetic_response/results.json',
              ROOT/'model/native_mass_channels/algebra.py', Path(__file__)]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    gamma, beta, spin, sectors = internal_matrices()
    max_links = max_spectrum = max_square = max_current = max_matrix_tail = 0.
    checked_cases = 0
    for n, a, field, charge, mass in ((3, .6, .09, -1.7, 4.), (5, .75, -.04, .5, 5.)):
        px, py = link_momenta(n, a, field, charge)
        actual = producer.covariant_momenta(n, a, field, charge)
        max_links = max(max_links, *(float(np.max(np.abs((p-q).toarray()))) for p, q in zip((px, py), actual[:2])))
        p, q = px.toarray(), py.toarray()
        d = n*n
        orbital = p@p+q@q
        curvature = 1j*(p@q-q@p)
        h = np.kron(gamma[0], p)+np.kron(gamma[1], q)+mass*np.kron(beta, np.eye(d))
        expected_square = np.kron(np.eye(16), orbital+mass*mass*np.eye(d))+np.kron(spin, curvature)
        max_square = max(max_square, float(np.max(np.abs(h@h-expected_square))))
        blocks = [orbital+sign*curvature for sign in (-1, 1)]
        union = []
        for block in blocks:
            values, vectors = eigh(block)
            assert values.min() > -1e-13 and values.max() <= 4/a**2
            roots = np.sqrt(mass*mass+values)
            union.extend(np.repeat(roots, 4)); union.extend(np.repeat(-roots, 4))
            # Stable direct spectral function avoids subtracting two large terms.
            exact = (vectors*(values/(roots+mass)))@vectors.conj().T
            polynomial, bound = producer.positive_energy_minus_mass(sparse.csr_matrix(block), mass, a)
            err = float(np.linalg.norm(exact-polynomial.toarray(), 2))
            assert err <= bound+3e-13
            max_matrix_tail = max(max_matrix_tail, err/bound)
        full = np.linalg.eigvalsh(h)
        max_spectrum = max(max_spectrum, float(np.max(np.abs(full-np.sort(union)))))
        assert np.linalg.norm(h@h-np.kron(np.eye(16), orbital+mass*mass*np.eye(d))) > .02
        # Challenge full-field current derivative for a different electric charge.
        step = 1e-6
        upper = link_momenta(n, a, field+step, charge)[1]
        lower = link_momenta(n, a, field-step, charge)[1]
        derivative = sparse.lil_matrix((d, d), dtype=complex)
        for x, y in product(range(n), range(n-1)):
            i = x*n+y; j = i+1
            coordinate = (x-n//2)*a
            value = -.5*charge*coordinate*np.exp(1j*charge*a*field*coordinate)
            derivative[i, j] = value; derivative[j, i] = value.conjugate()
        max_current = max(max_current, float(np.max(np.abs(((upper-lower)/(2*step)-derivative).toarray()))))
        assert np.linalg.norm(derivative.toarray()) > .1
        checked_cases += 1
    assert max_links < 1e-15 and max_square < 2e-14 and max_spectrum < 8e-14 and max_current < 5e-10

    # One emitted row is checked by exact dense spectral functional calculus,
    # independently of the producer's sixth-order matrix polynomial.
    n = 33; a = 1.; field = 1/512; mass = 64.
    p, q = [m.toarray() for m in link_momenta(n, a, field, 1.)]
    axis = np.arange(n)-n//2
    xx, yy = np.meshgrid(axis, axis, indexing='ij')
    state = (np.exp(-(xx*xx+yy*yy)/64)*np.cos(np.pi*xx/32)*np.cos(np.pi*yy/32)).ravel()
    state /= np.linalg.norm(state)
    orbital = p@p+q@q; curvature = 1j*(p@q-q@p)
    means = []
    for sign in (1, -1):
        values, vectors = eigh(orbital+sign*curvature)
        roots = values/(np.sqrt(mass*mass+values)+mass)
        means.append(float(np.dot(np.abs(vectors.conj().T@state)**2, roots)))
    exact_g = 2*mass*(means[0]-means[1])/field
    recorded = json.loads((ROOT/'model/openwave_magnetic_response/results.json').read_text())
    stored_g = recorded['refinement_rows'][0]['g_from_positive_branch']
    assert abs(exact_g-stored_g) < 1e-9
    assert 1.95 < exact_g < 1.97
    # Spin normalization is S=Sigma/2. Using Sigma itself changes the g convention.
    assert np.array_equal(spin@spin, np.eye(16))
    assert np.allclose((spin/2)@(spin/2), np.eye(16)/4)
    assert np.linalg.norm(np.kron(np.eye(16), curvature)-np.kron(spin, curvature)) > .01
    if pins != {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}:
        raise ValueError('magnetic inputs changed during independent audit')
    return {'schema': 'oph.independent_magnetic_response_audit.v1',
            'verdict': 'PASS_CONDITIONAL_SPIN_CURVATURE_RESPONSE', 'source_sha256': dict(sorted(pins.items())),
            'internal_dimension': 16, 'joint_sector_ranks': sectors,
            'independent_complete_spectrum_cases': checked_cases,
            'largest_unreduced_direct_matrix': 400,
            'independent_link_matrix_max_error': max_links,
            'complete_square_max_error': max_square,
            'full_signed_spectrum_max_error': max_spectrum,
            'nonunit_charge_current_derivative_error': max_current,
            'matrix_function_error_to_analytic_bound_max_ratio': max_matrix_tail,
            'scalar_sqrt_tail': sqrt_remainder(),
            'direct_1089_row_spectral_g': exact_g,
            'reported_g_error': abs(exact_g-stored_g),
            'controls': {'deleted_curvature_rejected': True, 'scalar_spin_current_rejected': True,
                         'all_signed_branches_and_valley_multiplicities_retained': True,
                         'nonunit_charge_and_field_sign_changes_checked': True},
            'scope': {'declared_gauged_envelope_response_confirmed': True,
                      'spin_pointer_is_in_positive_branch_frame': True,
                      'self_sourced_stationary_dipole_computed': False,
                      'physical_electron_or_Bohr_magneton_identified': False,
                      'anomalous_magnetic_moment_computed': False,
                      'analytic_tail_includes_floating_roundoff': False,
                      'finest_16641_row_response_independently_repeated': False}}


if __name__ == '__main__':
    output = Path(__file__).resolve().parents[1]/'magnetic_review.json'
    output.write_text(json.dumps(run(), indent=2, sort_keys=True, allow_nan=False)+'\n')
    print('PASS: independently reconstructed magnetic current, branches and response')
