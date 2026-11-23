from types import SimpleNamespace

import numpy as np
import pytest

from . import measure as m


def literal_operator(graph, route, charge, strength):
    """Direct embedding of the full 99-state edge term at every defect address."""
    exchange, agreement = m.emission.LAW.local_term_parts()
    native = exchange + strength * agreement
    arms = [tuple(sorted(route[:2])), tuple(sorted(route[2:]))]
    coins = np.column_stack(
        [m.emission.memory_vector([(c, -1, 1), (c, 1, 1)]) for c in (0, 1)]
    )
    blank = np.eye(m.emission.NM)[:, :1]
    n = graph.n
    words = 1 + charge * np.eye(n, dtype=int)
    h = np.kron(
        np.diag(strength * np.sum(words**2 * (6 - graph.degree), axis=1)), np.eye(4)
    )
    for edge in graph.seams:
        frame = coins if edge in arms else blank
        other = [i for i in range(n) if i not in edge]
        embeddings = []
        for word in words:
            embed = np.zeros((m.emission.LAW.DIM, frame.shape[1]), complex)
            embed[m.emission.columns(word[list(edge)])] = frame
            embeddings.append(embed)
        for i, left in enumerate(embeddings):
            for j, right in enumerate(embeddings):
                if not np.array_equal(words[i, other], words[j, other]):
                    continue
                local = (left.conj().T @ native @ right).real
                if edge == arms[0]:
                    term = np.kron(local, np.eye(2))
                elif edge == arms[1]:
                    term = np.kron(np.eye(2), local)
                else:
                    term = local.item() * np.eye(4)
                h[4 * i : 4 * i + 4, 4 * j : 4 * j + 4] += term
    return h - np.eye(4 * n) * h.trace() / (4 * n)


@pytest.mark.parametrize("charge", [-1, 1])
@pytest.mark.parametrize("strength", [0, m.emission.K])
def test_all_free_edges_agree_with_direct_native_columns(charge, strength):
    seams = [(0, 2), (0, 3), (1, 3), (2, 3), (1, 4), (2, 4)]
    graph = SimpleNamespace(
        n=5, seams=seams, degree=np.bincount(np.ravel(seams), minlength=5)
    )
    route = (2, 0, 3, 1)
    h = m.hamiltonian(graph, route, charge, strength).toarray()
    np.testing.assert_allclose(
        h, literal_operator(graph, route, charge, strength), atol=1e-13
    )


def test_operator_calibration_predicts_complex_coin_covariance():
    rng = np.random.default_rng(321)
    v = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    rho = v @ v.conj().T
    rho /= rho.trace()
    product = np.kron(*m.marginals(rho))
    for a, b in zip(m.marginals(product), m.marginals(rho), strict=True):
        np.testing.assert_allclose(a, b)
    _, cov = m.tomography(rho, rng)
    e = np.kron(m.ANALYZERS[1][0], m.ANALYZERS[0][1])
    row = m.read_probe([(1, "test", e)], rho, cov, rng)[0]
    assert (
        abs(row["forecast"] - np.trace(e @ (rho - product)).real)
        < row["forecast_radius_99"]
    )


def test_reference_uses_no_pseudocount_and_matches_symmetric_divergence():
    rng = np.random.default_rng(554)
    singular = m.modular_read({4.0: np.diag([1.0, 0, 0, 0])}, rng)
    assert not singular["faithful_99"] and "held_reads" not in singular
    p = np.array([[0.37, 0.13], [0.18, 0.32]])
    product = np.outer(p.sum(axis=1), p.sum(axis=0))
    prediction = -np.sum((p - product) * np.log(p / product))
    state = np.diag(p.ravel())
    measured = m.modular_read(dict.fromkeys([2.0, 4.0, 8.0], state), rng)
    assert measured["faithful_99"]
    for row in measured["held_reads"]:
        assert abs(prediction - row["measured_nats"]) < row["radius_99_nats"]


def test_blank_and_uniform_controls():
    result = m.uniform_controls(np.random.default_rng(541), [(0, 1)])
    assert (
        result["blank_commutator_nonzero"]
        == result["magnitude_commutator_nonzero"]
        == 0
    )
    assert result["arms"][0]["uniform_seam_energy"] == 2 * m.emission.K
    assert all(row["zz_covariance"] == 0.0 for row in result["arms"][0]["appointments"])
