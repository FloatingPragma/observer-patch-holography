"""Independent full-CAR and reference tests for physical charged observables."""
from dataclasses import replace
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import chain, observables as o
from ..interacting_field.lattice import SIGMA
from ..particle_response import dimer
from ..full_fock_field.field import FieldGeometry


def literal_car(mode):
    result = np.zeros((256, 256), complex)
    for mask in range(256):
        if mask >> mode & 1:
            result[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return result


def literal_chain(axis=0):
    ladders = [literal_car(mode) for mode in range(8)]
    numbers = [x.conj().T@x for x in ladders]
    result = -sum(numbers)/2
    for site in range(4):
        result += numbers[2*site]@numbers[2*site+1]
    for site in range(3):
        for spin in range(2):
            for other in range(2):
                term = 1j*SIGMA[axis][spin, other]/16*ladders[2*site+spin].conj().T@ladders[2*(site+1)+other]
                result += term+term.conj().T
    return result


@pytest.mark.parametrize('axis', [0, 1, 2])
@pytest.mark.parametrize('kind', ['hole', 'particle'])
def test_same_local_residue_predicts_joining_current_without_fitting(axis, kind):
    result = o.projected_coefficient(kind, axis)
    spectrum = dimer.setup(FieldGeometry(tuple(2 if j == axis else 1 for j in range(3))))
    projector = dimer.candidate_projectors(spectrum)['charge_minus_one_hole']
    removed = spectrum.space.annihilator(0).toarray()@spectrum.vectors[:, 0]
    residue = float(np.vdot(removed, projector@removed).real)
    assert abs(result['coefficient']-residue/16) < 2e-14
    assert abs(residue-(1+1/math.sqrt(17))/4) < 2e-14
    bare = chain.bare_columns(kind, axis)
    number_a = np.diag([(mask & 15).bit_count() for mask in range(256)])
    k = literal_chain(axis)
    independent_current = -1j*(k@number_a-number_a@k)
    np.testing.assert_allclose(bare.conj().T@independent_current@bare, result['current'],
                               atol=2e-14, rtol=0)


@pytest.mark.parametrize('kind', ['hole', 'particle'])
def test_dressed_logical_transfer_has_its_actual_small_charge_contrast(kind):
    data = o.dressed_observables(kind)
    charge = -1 if kind == 'hole' else 1
    assert .006 < data['hopping'] < .007
    assert .07 < data['charge_contrast'] < .09
    assert data['charge_read_offband_operator_norm'] > .5
    assert data['bare_overlap_eigenvalues'][0] < .28
    before = o.history(kind, 'dressed', 0.)
    after = o.history(kind, 'dressed', data['transfer_time'])
    np.testing.assert_allclose(after['dressed_destination_weight'], np.eye(2), atol=2e-12, rtol=0)
    np.testing.assert_allclose(before['QA']-after['QA'], charge*data['charge_contrast']*np.eye(2),
                               atol=2e-12, rtol=0)
    assert np.linalg.norm(after['QA']) > .6  # No complete physical charge transfer.
    np.testing.assert_allclose(data['current'],
        -1j*(data['generator']@data['QA']-data['QA']@data['generator']), atol=2e-12, rtol=0)


@pytest.mark.parametrize('kind', ['hole', 'particle'])
def test_bare_full_history_and_every_record_preserve_an_unknown_reference(kind):
    actual = o.history(kind, 'bare', 4.)
    expected = expm(-4j*literal_chain())@actual['initial']
    coefficients = np.array([[1., 2j, -1j], [2., -1., 3j]])
    coefficients /= np.linalg.norm(coefficients)
    np.testing.assert_allclose(actual['state']@coefficients, expected@coefficients, atol=2e-12, rtol=0)
    assert len(actual['occupation_kraus_rows']) == 256
    gram = np.zeros((3, 3), complex)
    for address, row in enumerate(actual['occupation_kraus_rows']):
        assert row.shape == (1, 2)
        branch = row@coefficients
        gram += branch.conj().T@branch
        np.testing.assert_allclose(row, expected[address:address+1], atol=2e-12, rtol=0)
    np.testing.assert_allclose(gram, coefficients.conj().T@coefficients, atol=2e-12, rtol=0)
    assert actual['block_number_probabilities'].shape == (5, 5, 2)
    np.testing.assert_allclose(actual['block_number_probabilities'].sum(axis=(0, 1)), 1., atol=2e-12)
    assert abs(actual['band_weight'].trace().real/2-.5109883376314) < 2e-12
    assert abs(actual['bare_destination_weight'].trace().real/2-.0058515268725) < 2e-12
    assert actual['right_exterior_to_product_ground'] > .03


def test_bare_charge_and_candidate_transfer_at_full_band_time_are_not_projected_predictions():
    time = o.dressed_observables()['transfer_time']
    actual = o.history('hole', 'bare', time)
    expected = expm(-1j*time*literal_chain())@actual['initial']
    np.testing.assert_allclose(actual['state'], expected, atol=2e-11, rtol=0)
    assert .18 < -actual['QA'].trace().real/2 < .20
    assert .77 < actual['bare_destination_weight'].trace().real/2 < .78
    assert .48 < 1-actual['band_weight'].trace().real/2 < .50


def test_complete_coupled_ground_measure_uses_all_poles_and_no_band_renormalization():
    response = o.charged_response()
    k = literal_chain()
    energies, vectors = np.linalg.eigh(k)
    ground = vectors[:, 0]
    assert response['addition_weights'].shape == (8, 256)
    assert response['removal_weights'].shape == (8, 256)
    assert len(response['frequencies']) == 256
    for mode in range(8):
        a = literal_car(mode)
        for weights, vector in ((response['addition_weights'][mode], a.conj().T@ground),
                                (response['removal_weights'][mode], a@ground)):
            assert abs(weights.sum()-.5) < 2e-13
            for duration in (0., .37, 4., 29.):
                measured = weights@np.exp(-1j*duration*response['frequencies'])
                independent = np.vdot(vector, expm(-1j*duration*(k-energies[0]*np.eye(256)))@vector)
                assert abs(measured-independent) < 2e-12
    assert max(response['addition_band_weights']) < .25
    assert min(response['addition_weights'].sum(axis=1)-response['addition_band_weights']) > .25
    np.testing.assert_allclose(response['addition_band_weights'], response['removal_band_weights'],
                               atol=2e-13, rtol=0)


@pytest.mark.parametrize('kwargs', [dict(axis=True), dict(axis=1.), dict(axis=3), dict(kind=True),
    dict(kind='electron'), dict(preparation='ground'), dict(preparation=True),
    dict(duration=True), dict(duration=float('nan')), dict(duration=float('inf')), dict(duration=513)])
def test_history_rejects_coercion_and_undeclared_domains(kwargs):
    with pytest.raises(ValueError):
        o.history(**kwargs)


@pytest.mark.parametrize('field', ['K', 'energies', 'vectors', 'sectors', 'ground_projector'])
def test_fresh_canonical_validation_rejects_mutated_complete_spectrum(monkeypatch, field):
    original = chain.setup
    def changed(*args, **kwargs):
        data = original(*args, **kwargs)
        value = np.array(getattr(data, field), copy=True)
        if field == 'sectors':
            value = value.astype(float)
        elif field == 'energies':
            value[-1] = np.nan
        else:
            value[0, 0] += .01
        return replace(data, **{field: value})
    o.operators()  # Warm every predecessor before corruption.
    monkeypatch.setattr(chain, 'setup', changed)
    with pytest.raises(ArithmeticError):
        o.charged_response()


def test_mutated_band_cannot_relabel_dark_complement_as_charged_residue(monkeypatch):
    original = chain.band
    def changed(*args, **kwargs):
        band = original(*args, **kwargs)
        return replace(band, projector=np.eye(256))
    o.charged_response()
    monkeypatch.setattr(chain, 'band', changed)
    with pytest.raises(ArithmeticError):
        o.charged_response()


def test_mutated_bare_input_cannot_replace_same_candidate_with_another_charge_state(monkeypatch):
    wrong = np.zeros((256, 4), complex)
    wrong[[7, 11, 13, 14], range(4)] = 1
    o.initial_columns()
    monkeypatch.setattr(chain, 'bare_columns', lambda *args, **kwargs: wrong)
    with pytest.raises(ArithmeticError):
        o.history('hole', 'bare', 4.)


def test_returned_state_arrays_do_not_poison_later_fresh_calls():
    first = o.history()
    with pytest.raises(ValueError):
        first['state'][0, 0] = 1
    first['state'] = np.zeros((256, 2))
    later = o.history()
    np.testing.assert_allclose(later['state'].conj().T@later['state'], np.eye(2), atol=2e-12)
