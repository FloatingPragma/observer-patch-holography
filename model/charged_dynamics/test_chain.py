"""Independent occupation and complete-carrier controls for the joined action."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import chain as c
from ..particle_response import dimer
from ..full_fock_field.field import FieldGeometry


def annihilator(modes, mode):
    out = np.zeros((1 << modes, 1 << modes), complex)
    for mask in range(1 << modes):
        if mask >> mode & 1:
            out[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return out


def independent_generator(axis, coupling):
    sigma = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))[axis]
    ladders = [annihilator(8, j) for j in range(8)]
    out = np.zeros((256, 256), complex)
    for site in range(3):
        for a in range(2):
            for b in range(2):
                term = (1j/16)*sigma[a, b]*(float(coupling) if site == 1 else 1)*ladders[2*site+a].conj().T@ladders[2*(site+1)+b]
                out += term+term.conj().T
    out += np.diag([sum((mask >> (2*j) & 3) == 3 for j in range(4))-mask.bit_count()/2 for mask in range(256)])
    return out


@pytest.mark.parametrize('axis', range(3))
@pytest.mark.parametrize('coupling', (F(0), F(1, 3), F(1)))
def test_complete_action_matches_independent_jordan_wigner_entries(axis, coupling):
    model = c.setup(axis, coupling)
    expected = independent_generator(axis, coupling)
    np.testing.assert_array_equal(model.K, expected)
    assert len(model.energies) == 256
    assert sorted(model.sectors.tolist()) == sorted(mask.bit_count() for mask in range(256))
    np.testing.assert_allclose(model.K@model.vectors, model.vectors*model.energies, atol=2e-13, rtol=0)
    assert model.sectors[0] == 4 and model.energies[1]-model.energies[0] > 1/128
    np.testing.assert_allclose(model.ground_projector@model.ground_projector, model.ground_projector, atol=2e-14, rtol=0)


@pytest.mark.parametrize('axis', range(3))
def test_local_columns_are_the_existing_complete_dimer_projectors(axis):
    sides = [1, 1, 1]; sides[axis] = 2
    spectrum = dimer.setup(FieldGeometry(tuple(sides)))
    ground = c.ground_column(axis)
    np.testing.assert_allclose(np.outer(ground, ground.conj()), spectrum.reference, atol=2e-14, rtol=0)
    assert ground[3].real < 0 and ground[3].imag == 0
    projectors = dimer.candidate_projectors(spectrum)
    for kind, name in (('hole', 'charge_minus_one_hole'), ('particle', 'charge_plus_one_particle')):
        columns = c.local_columns(kind, axis)
        np.testing.assert_allclose(columns@columns.conj().T, projectors[name], atol=2e-14, rtol=0)
        bare = c.bare_columns(kind, axis)
        np.testing.assert_allclose(bare.conj().T@bare, np.eye(4), atol=2e-14, rtol=0)
        model = c.setup(axis, F(0))
        energy = np.trace(bare.conj().T@model.K@bare)/4
        np.testing.assert_allclose(model.K@bare, energy*bare, atol=2e-14, rtol=0)
        np.testing.assert_allclose(c.band(kind, axis, F(0)).projector, bare@bare.conj().T, atol=2e-13, rtol=0)


@pytest.mark.parametrize('axis', range(3))
def test_exact_monomial_gauge_on_every_number_sector(axis):
    transform = c.gauge(axis)
    base = np.zeros((256, 256), int); link = np.zeros_like(base)
    for up in range(5):
        for down in range(5):
            addresses, a, b = c.hubbard_integer_parts(up, down)
            base[np.ix_(addresses, addresses)] = a
            link[np.ix_(addresses, addresses)] = b
    np.testing.assert_array_equal(transform.conj().T@transform, np.eye(256))
    np.testing.assert_array_equal(transform@(base+link)@transform.conj().T/16,
                                  independent_generator(axis, F(1)))


@pytest.mark.parametrize('kind', ('hole', 'particle'))
def test_uniform_isolated_band_and_dressing_keep_the_full_complement(kind):
    model = c.setup()
    selected = np.flatnonzero(model.sectors == (3 if kind == 'hole' else 5))
    band = c.band(kind)
    assert band.gap_to_complement > 1/256
    assert abs(band.gap_to_complement-.006645502416206) < 2e-13
    assert len(selected) == 56 and len(band.energies) == 4
    dressed, bare = c.dressed_columns(kind), c.bare_columns(kind)
    np.testing.assert_allclose(dressed.conj().T@dressed, np.eye(4), atol=2e-13, rtol=0)
    np.testing.assert_allclose(dressed@dressed.conj().T, band.projector, atol=2e-13, rtol=0)
    overlap = np.linalg.eigvalsh(bare.conj().T@band.projector@bare)
    np.testing.assert_allclose(overlap, [.2704793225]*2+[.7514973527]*2, atol=5e-11, rtol=0)
    assert np.linalg.norm((np.eye(256)-band.projector)@bare) > 1


def test_mutated_returned_metadata_does_not_poison_later_models():
    first = c.setup()
    with pytest.raises(ValueError):
        first.K[0, 0] = 1
    # No cache exposes the arrays of a later model, even through reshape/base.
    assert not np.shares_memory(first.K, c.setup().K)
    assert not np.shares_memory(c.bare_columns(), c.bare_columns())


@pytest.mark.parametrize('nonfinite', (float('nan'), float('inf')))
def test_complete_spectrum_rejects_nonfinite_vectors_with_finite_energies(monkeypatch, nonfinite):
    original = np.linalg.eigh
    def changed(matrix, *args, **kwargs):
        values, vectors = original(matrix, *args, **kwargs)
        # Eight-dimensional blocks belong to the chain, not the dimer helper.
        if matrix.shape == (8, 8):
            vectors[:, 0] = nonfinite
        return values, vectors
    monkeypatch.setattr(np.linalg, 'eigh', changed)
    with pytest.raises(ArithmeticError, match='complete chain spectral admission'):
        c.setup()


@pytest.mark.parametrize('nonfinite', (float('nan'), float('inf')))
def test_polar_dressing_rejects_nonfinite_vectors_with_positive_gram_spectrum(monkeypatch, nonfinite):
    original = np.linalg.eigh
    def changed(matrix, *args, **kwargs):
        values, vectors = original(matrix, *args, **kwargs)
        # The positive four-column Gram matrix is distinct from the dimer's
        # negative-energy four-dimensional number blocks.
        if matrix.shape == (4, 4) and values[0] > 0:
            vectors[:, 0] = nonfinite
        return values, vectors
    monkeypatch.setattr(np.linalg, 'eigh', changed)
    with pytest.raises(ArithmeticError, match='polar chart is singular or unresolved'):
        c.dressed_columns()


@pytest.mark.parametrize('axis', (True, 0., -1, 3))
def test_axis_validation(axis):
    with pytest.raises(ValueError): c.setup(axis)
    with pytest.raises(ValueError): c.ground_column(axis)


@pytest.mark.parametrize('coupling', (True, .5, -1, 2, '1'))
def test_exact_link_domain(coupling):
    with pytest.raises(ValueError): c.setup(link_strength=coupling)


def test_malformed_rational_and_changed_shared_spin_chart_reject(monkeypatch):
    coupling = F(0); coupling._denominator = 0
    with pytest.raises(ValueError): c.setup(link_strength=coupling)
    changed = [matrix.copy() for matrix in c.SIGMA]
    changed[0][0, 0] = 1
    monkeypatch.setattr(c, 'SIGMA', tuple(changed))
    with pytest.raises(ArithmeticError, match='spin chart'): c.setup()
