"""Independent matrix spectra, exact limits and degenerate source controls."""
from fractions import Fraction as F
from itertools import product

import mpmath as mp
import numpy as np
import pytest
import sympy as s

from model.native_mass_channels import algebra as source
from . import algebra, bounds


def test_all_source_columns_and_exact_projector_resolution():
    assert algebra.source_join()['source_rows'] == 128
    for kinetic, masses in [((1, 2, 3), (2, 1, 2, 2)), ((0, 0, 0), (1, 0, 0, 1)),
                            ((0, 0, 0), (0, 0, 0, 0)), ((1, 0, 2), (2, 0, 0, 0))]:
        rows = algebra.spectral_projectors(kinetic, masses)
        h = source.folded_hamiltonian(kinetic, masses)
        total, resolution = s.zeros(16), s.zeros(16)
        for _, energy, projector in rows:
            assert (projector*projector-projector).applyfunc(s.simplify) == s.zeros(16)
            assert (h*projector-energy*projector).applyfunc(s.simplify) == s.zeros(16)
            total += projector; resolution += energy*projector
        assert total.applyfunc(s.simplify) == s.eye(16)
        assert (resolution-h).applyfunc(s.simplify) == s.zeros(16)
    assert algebra.clifford_certificate()['Clifford_checks'] == 16


@pytest.mark.parametrize('masses', [(0, 0, 0, 0), (2, 0, 0, 0), (1, 0, 0, 1),
                                    (0, 1, 2, 2), (2, 1, 1, 0), (-1, 0, 0, 1)])
def test_energy_intervals_against_high_precision_and_literal_matrix(masses):
    mp.mp.dps = 70
    for momentum in ((0, 0, 0), (1, 0, 0), (F(-1, 3), F(2, 5), F(3, 4))):
        a, v = F(1, 8), F(3, 2)
        report = bounds.energy_rows(momentum, a, v, masses)
        exact = lambda x: mp.mpf(F(x).numerator)/F(x).denominator
        radial = mp.sqrt(sum(exact(x)**2 for x in masses[1:]))
        kinetic = sum((exact(v)/exact(a)*mp.sin(exact(a)*exact(x)))**2 for x in momentum)
        expected = []
        for row in report['rows']:
            value = row['energy_sign']*mp.sqrt(kinetic+(exact(masses[0])+row['mass_sector']*radial)**2)
            lo, hi = map(exact, row['lattice_energy_interval'])
            assert lo <= value <= hi
            expected.extend([float(value)]*row['rank'])
        # Build every16-row entry independently from valley flips and Pauli matrices.
        pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
        bitrows = tuple(product((0, 1), repeat=3))
        h = np.zeros((16, 16), complex)
        for j, bits in enumerate(bitrows):
            sl = slice(2*j, 2*j+2)
            h[sl, sl] = -float(v/a)*sum((-1)**bits[t]*np.sin(float(a*momentum[t]))*pauli[t] for t in range(3))
            opposite = tuple(1-b for b in bits); row = bitrows.index(opposite)
            h[2*row:2*row+2, sl] += float(masses[0])*np.eye(2)
            for t in range(3):
                flipped = tuple(b ^ int(t == q) for q, b in enumerate(bits)); row = bitrows.index(flipped)
                h[2*row:2*row+2, sl] += float(masses[t+1])*pauli[t]
        # Floating eigensolver admission is separate from the70-digit checks
        # of unchanged exact rational intervals above.
        assert np.max(np.abs(np.linalg.eigvalsh(h)-np.sort(expected))) < 1e-12*max(1, np.linalg.norm(h, 2))


def test_uniform_shell_errors_and_quadratic_refinement():
    p = (F(1, 3), F(2, 3), F(1))
    defects = []
    for a in (F(1, 2), F(1, 4), F(1, 8), F(1, 16)):
        report = bounds.energy_rows(p, a, 1, (0, 0, 0, 3))
        row = next(r for r in report['rows'] if r['energy_sign'] == 1)
        lattice = np.sqrt(9+sum(np.sin(float(a*x))**2/float(a*a) for x in p))
        continuum = np.sqrt(9+sum(float(x*x) for x in p))
        defect = lattice*lattice-continuum*continuum
        assert float(F(row['mass_shell_defect_interval'][0]))-1e-14 <= defect <= 1e-14
        assert continuum-lattice <= float(F(row['absolute_energy_error_upper']))
        defects.append(-defect)
    assert all(3.8 < x/y < 4.01 for x, y in zip(defects, defects[1:]))


@pytest.mark.parametrize('bad', [True, .5, float('nan'), '1/2'])
def test_inputs_do_not_silently_round_source_parameters(bad):
    with pytest.raises(ValueError):
        bounds.energy_rows((0, 0, 0), bad)
    with pytest.raises(ValueError):
        bounds.energy_rows((0, 0, bad))


def test_zero_mass_tuning_retains_entire_zero_subspace():
    rows = bounds.energy_rows((0, 0, 0), coefficients=(F(1, 3), 0, 0, F(1, 3)))['rows']
    zero = [r for r in rows if r['energy_sign'] == 0]
    assert len(zero) == 1 and zero[0]['rank'] == 8
    assert zero[0]['lattice_energy_interval'] == ['0', '0']
    assert sum(r['rank'] for r in rows) == 16
