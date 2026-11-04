from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import null_space

from ..topological_defect_spectrum import sampling
from . import algebra, finite


def test_complete_original_source_matrix_against_frozen_literal_sampling():
    points = np.array(finite.sites(1))
    addresses = {tuple(point): j for j, point in enumerate(points)}
    rng = np.random.default_rng(314)
    columns = rng.normal(size=(len(points), 2, 3))+1j*rng.normal(size=(len(points), 2, 3))
    a = F(1, 2)
    def field(x):
        result = np.zeros((len(x), 16, 3), complex)
        for j, point in enumerate(np.rint(x/float(a)).astype(int)):
            address = addresses.get(tuple(point))
            if address is not None:
                result[j, :2] = columns[address]/float(a)**1.5
        return result
    expected = sampling.lattice_action(points, a, field, velocity=2, mass=3).reshape(-1, 3)
    actual = finite.source_matrix(1, a, 2, 3)@columns.reshape(-1, 3)
    np.testing.assert_allclose(actual, expected, atol=3e-14, rtol=0)


@pytest.mark.parametrize('radius', [1, 2])
def test_full_inversion_spaces_and_every_homotopy_endpoint(radius):
    plus, minus = finite.inversion_bases(radius)
    unitary = np.column_stack((plus.toarray(), minus.toarray()))
    np.testing.assert_allclose(unitary.conj().T@unitary, np.eye(unitary.shape[0]), atol=3e-16, rtol=0)
    assert plus.shape[1]-minus.shape[1] == 2
    for sign in (-1, 1):
        for fraction in (0, F(3, 10), 1):
            h = finite.source_matrix(radius, sign=sign, hopping_fraction=fraction)
            assert np.linalg.norm((h-h.conj().T).toarray()) == 0
            assert np.linalg.norm((plus.conj().T@h@plus).toarray()) < 1e-14
            assert np.linalg.norm((minus.conj().T@h@minus).toarray()) < 1e-14
            off = (minus.conj().T@h@plus).toarray()
            zero = plus@null_space(off)
            assert zero.shape[1] >= 2
            assert np.linalg.norm(h@zero) < 1e-12
            values = np.linalg.eigvalsh(h.toarray())
            assert np.count_nonzero(abs(values) < 1e-12) >= 2


def test_strong_compact_odd_perturbations_and_scalar_countercontrols():
    h = finite.source_matrix(2)
    plus, minus = finite.inversion_bases(2)
    for coupling in (3, 100):
        changed = h+finite.compact_odd_density(2, coupling)
        assert np.linalg.norm((plus.conj().T@changed@plus).toarray()) < 1e-12
        assert np.linalg.norm((minus.conj().T@changed@minus).toarray()) < 1e-12
        values = np.linalg.eigvalsh(changed.toarray())
        assert np.count_nonzero(abs(values) < 1e-10) >= 2
    zero = plus@null_space((minus.conj().T@h@plus).toarray())
    shifted = h.toarray()+np.eye(h.shape[0])/7
    assert np.linalg.norm(shifted@zero-zero/7) < 1e-12
    assert np.linalg.norm(shifted@zero) > .1
    # The compact fixed-site countercontrol has an exactly known atomic
    # spectrum: two levels at 1/3, every other level at +/-1, and no zero.
    atomic = finite.source_matrix(1, hopping_fraction=0).toarray()
    origin = finite.sites(1).index((0, 0, 0))
    atomic[2*origin:2*origin+2, 2*origin:2*origin+2] = np.eye(2)/3
    values = np.linalg.eigvalsh(atomic)
    assert np.min(abs(values)) > .333333333333
    assert np.count_nonzero(abs(values-1/3) < 1e-14) == 2


def test_exact_source_certificates_keep_the_two_different_indices_distinct():
    value = algebra.certificate()
    assert value['exact_Fredholm_index'] == 2
    assert value['atomic']['reflection_index_plus_hedgehog'] == 2
    assert value['atomic']['reflection_index_antihedgehog'] == 2
    assert value['atomic']['identified_with_continuum_internal_chiral_index_or_degree'] is False
    assert value['exact_kernel_exhaustion_proved'] is False
    assert value['native_repair_preserves_required_operator_symmetry_proved'] is False
    for invalid in (True, 0, 2, 1.):
        with pytest.raises(ValueError):
            algebra.orientation(invalid)
    with pytest.raises(ValueError):
        algebra.onsite((0, True, 1))
    for invalid in (True, -1, F(3, 2)):
        with pytest.raises(ValueError):
            finite.source_matrix(1, hopping_fraction=invalid)
