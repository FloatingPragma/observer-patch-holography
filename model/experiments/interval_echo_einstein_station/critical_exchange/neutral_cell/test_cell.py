import numpy as np
from scipy.linalg import eigh, expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import cell
from model.experiments.interval_echo_gauge.port_phases.apparatus import path, swap_word
from . import measure as m


def test_exact_native_transport_and_unconditional_selection():
    board = cell()
    contact = candidate.native_law.seam_term(candidate.KAPPA)
    swap = np.eye(9)[[3*b + a for a in range(3) for b in range(3)]]
    np.testing.assert_allclose(expm(-1j * np.pi * contact), swap, atol=5e-14, rtol=0)
    pairs = list(m.combinations(m.STOCK, 2)) + [(m.DONOR, s) for s in m.STATIONS]
    for a, b in pairs:
        steps = path(board, a, b, set())
        labels = list(range(12))
        swap_word(labels, steps + steps[-2::-1])
        intended = list(range(12))
        intended[a], intended[b] = b, a
        assert labels == intended
    assert .5 * abs(m.stock_source() - 1/56).sum() < 1.648e-9


def test_independent_excess_propagation_with_direct_projectors():
    board = cell()
    route = (0, 1)
    material = (2, 4)
    basis = m.words(1, 2)
    h = m.native_generator(board, basis, 8)
    e, v = eigh(h / 2)
    u = (v * np.exp(-1j * m.STEP * e)) @ v.T
    sender, receiver = (np.array(basis)[:, j] == 2 for j in route)
    source = np.zeros(len(e), complex)
    source[basis.index(m.launch(material, route[0], 2))] = 1
    expected = np.zeros(m.TICKS + 2)
    for a in range(1, m.TICKS):
        source = u @ source
        branch = source * receiver
        source[receiver] = 0
        for b in range(1, m.TICKS+1-a):
            branch = u @ branch
            expected[a+b] += np.linalg.norm(branch[sender])**2
            branch[sender] = 0
    expected[-1] = 1 - expected.sum()
    reconstructed, _ = m.flight(board, 8, True)
    actual = reconstructed[route][m.MATERIAL.index(material[::-1])]
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)


def test_held_owner_and_charge_control_with_frozen_counts():
    with np.load(m.HERE / "counts.npz", allow_pickle=False) as bank:
        data, _ = m.run(8, archive=bank)
    assert data["maximum_two_charge_history_TV_on_sign_blind_heralds"] < 1.32e-8
    assert data["measured_local_seam_energy"]["near"] - data["measured_local_seam_energy"]["quiet"] > 45
    assert data["unmonitored_nonunit_at_mark"]["2.0"]["near"] < 0.31
    for g in ("near", "quiet"):
        lo, hi = data["clock_product_intervals_99_percent"][g]
        assert lo < 1 < hi
    bounds = data["clock_mean_intervals_simultaneous_99_percent"]
    for q in ("vacancy", "excess"):
        a, b = (bounds[key][q]["probe"] for key in ("1-0", "1-3"))
        assert a["near"][0] / a["quiet"][1] > b["near"][1] / b["quiet"][0]
    v, e = (bounds["0-1"][q]["probe"]["vacancy_at_mark"] for q in ("vacancy", "excess"))
    assert v[0] - e[1] > 0.10
