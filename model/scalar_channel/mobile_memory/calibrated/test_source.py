import itertools

import numpy as np
from scipy import sparse
from scipy.linalg import expm

from ....contrast_memory_spherical_v1 import carrier
from ..mobility import home_seam
from .source import build, leakage_witness


def _literal_source(model):
    """Enumerate the physical transitions without sparse product assembly."""
    coordinates = list(itertools.product(range(2), range(3), range(3)))
    labels = {address: i for i, address in enumerate(coordinates)}
    hamiltonian = np.zeros((len(coordinates), len(coordinates)))
    for i, (load, q, c) in enumerate(coordinates):
        hamiltonian[i, i] = (
            0.5 + model.charge_mass * (1 if q < 2 else -1)
            + model.cage_mass * (1 if c < 2 else -1)
            + model.coupling * (load == 0 and q == c == 2)
        )
        hamiltonian[labels[(1 - load, q, c)], i] = -0.5
        for endpoint, orientation in ((0, 1), (1, -1)):
            if q == endpoint or q == 2:
                other = endpoint + 2 - q
                hamiltonian[labels[(load, other, c)], i] = orientation * model.speed / model.length
            if c == endpoint or c == 2:
                other = endpoint + 2 - c
                hamiltonian[labels[(load, q, other)], i] = orientation * model.speed / model.length
    return hamiltonian


def test_literal_chart_and_spectral_leakage_on_each_native_level():
    for level in (0, 1, 2):
        federation = carrier.federation(level)
        lengths = {home_seam(federation): 1.25 / 2 ** level}
        model = build(federation, lengths, 0.83, selected_seam_only=True)
        assert model.shape == (2, 3, 3)
        assert all(pair in federation.seams for pair in model.pairs)
        assert model.swap_seam in federation.seams
        np.testing.assert_allclose(model.hamiltonian.toarray(), _literal_source(model), atol=2e-10, rtol=0)
        result = leakage_witness(model)
        assert result["negative_charge_probability_lower_bound"] > 0
        assert result["negative_charge_probability"] >= result["negative_charge_probability_lower_bound"]
        control = build(federation, lengths, 0.83, selected_seam_only=True, coupling=0)
        assert leakage_witness(control)["negative_charge_probability"] < 1e-20


def test_same_source_rest_controls_keep_the_tag_phase_during_load_exchange():
    federation = carrier.federation(0)
    model = build(federation, {home_seam(federation): 1.25}, 0.83, selected_seam_only=True)
    wait, split = 0.37, 0.13
    for seed, grade in ((model.charge_rest, model.charge_grading), (model.cage_rest, model.cage_grading)):
        H = model.hamiltonian.toarray()
        lower = expm(-1j * wait * (H - np.diag(grade * split / 2))) @ seed
        upper = expm(-1j * wait * (H + np.diag(grade * split / 2))) @ seed
        np.testing.assert_allclose(np.vdot(lower, upper), np.exp(-1j * split * wait), atol=2e-10, rtol=0)
        assert np.linalg.norm(lower.reshape(model.shape)[1]) > 0.05


def test_full_metric_kinetics_is_retained_by_both_contact_controls():
    federation = carrier.federation(0)
    lengths = {seam: 1.0 for seam in federation.external}
    active = build(federation, lengths, 1.0)
    free = build(federation, lengths, 1.0, coupling=0)
    parked = build(federation, lengths, 1.0, parked_swap=True)
    assert active.shape == (2, 50, 50)
    for control in (free, parked):
        np.testing.assert_array_equal(active.charge_operator.toarray(), control.charge_operator.toarray())
        np.testing.assert_array_equal(active.cage_operator.toarray(), control.cage_operator.toarray())
        np.testing.assert_allclose(active.initial, control.initial, atol=2e-10, rtol=0)
    contact = active.hamiltonian - free.hamiltonian
    contact.eliminate_zeros()
    row, col = contact.nonzero()
    np.testing.assert_array_equal(row, [active.contact_index])
    np.testing.assert_array_equal(col, [active.contact_index])
    np.testing.assert_allclose(contact[active.contact_index, active.contact_index], 1.0, atol=2e-10, rtol=0)
    swap = sparse.kron(sparse.csr_matrix(active.load_operator), sparse.eye(active.shape[1] ** 2))
    assert np.max(abs((swap @ contact - contact @ swap).data)) > 0.45
