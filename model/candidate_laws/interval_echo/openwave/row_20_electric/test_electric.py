"""Independent native-word, physical-pulse and endpoint-count checks."""

from itertools import product
import json

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from . import experiment as ex, neutral, schedule


def word_matrices(n, left, right, k):
    local = sum(x * y for x, y in zip(law.local_term_parts(), (1, k), strict=True))
    words = list(product((0, 2), repeat=n))
    lookup = {m: i for i, m in enumerate(words)}
    intra = np.zeros((2**n, 2**n))
    inter = np.zeros_like(intra)
    for col, w in enumerate(words):
        ports = (left + 1, *(value for a in w for value in (a, 2 - a)), right + 1)
        for i, (a, b) in enumerate(zip(ports, ports[1:])):
            idx = law.SRC.index((a, b))
            if i % 2:
                j = i // 2
                intra[col, col] += local[idx, idx]
                row = lookup[w[:j] + (2 - w[j],) + w[j + 1 :]]
                intra[row, col] += local[law.SRC.index((b, a)), idx]
            else:
                inter[col, col] += local[idx, idx]
    return intra, inter


@pytest.mark.parametrize("k", [0, 0.125, 8])
def test_static_three_state_energy_menu(k):
    j, h = ex.couplings(k)
    for n in (2, 3, 5):
        physical = []
        mapped = []
        for a, b in product((-1, 0, 1), repeat=2):
            intra, inter = word_matrices(n, a, b, k)
            physical.append(np.linalg.eigvalsh(intra + inter)[0])
            mapped.append(neutral.minimum(n, a, b, j, h))
        for vals in (physical, mapped):
            x = np.asarray(vals).reshape(3, 3)
            vals[:] = (x + x[1, 1] - x[:, 1, None] - x[None, 1, :]).ravel()
        np.testing.assert_allclose(physical, mapped, atol=3e-12)


def test_actual_contact_pulse_phase_and_critical_bounds():
    j, h = ex.couplings(8)
    for n in (2, 3, 6):
        for factor in (0.5, 1, 2):
            f = factor * h / j
            phase = []
            for r in (-1, 1):
                intra, inter = word_matrices(n, 1, r, 8)
                # Restore common phases only after differencing the two arms.
                a = intra - np.trace(intra) / len(intra) * np.eye(len(intra))
                b = inter - np.trace(inter) / len(inter) * np.eye(len(inter))
                u = expm(-1j * f * schedule.DELTA * b) @ expm(-1j * schedule.DELTA * a)
                phase.append(min(-np.angle(np.linalg.eigvals(u))))
            diff = phase[0] - phase[1]
            np.testing.assert_allclose(
                diff, (-1) ** n * schedule.parity_phase(n, j, h, f), atol=1e-12
            )
    theta = 2 * h * schedule.DELTA
    for r in json.loads((ex.HERE / "schedule.json").read_text())["schedules"][1][
        "rows"
    ]:
        m = r["n"] + 1
        v = 2 * np.tan(theta / 2) * np.sin(np.pi / (4 * m + 2))
        read = np.sin(abs(r["phase"]) / 2)
        assert (
            v * np.cos(theta / 2) ** 2 - 1e-11 <= read <= v * np.cos(theta / 2) + 1e-11
        )


def test_gates_meter_odd_sign_and_held_counts():
    stored = json.loads((ex.HERE / "records.json").read_text())
    for sign, meter in zip((-1, 0, 1), stored["local_meter"], strict=True):
        assert sum(meter) == meter[sign + 2] == ex.SHOTS
    errors = []
    for b in stored["boards"]:
        for row in b["shortest_route_reads"]:
            open_, shut = row["open_and_held"]
            assert sum(open_) == open_[-1] == 4096
            assert sum(shut) == 0
        for trial in b["studies"]:
            for row in trial["rows"]:
                p = row["source_minus_zero_plus"]
                for s, r in zip((-1, 0, 1), p, strict=True):
                    counts = np.array(r["counts_xy"])
                    complex_signal = (2 * counts[0] / ex.SHOTS - 1) + 1j * (
                        2 * counts[1] / ex.SHOTS - 1
                    )
                    # The microscopic energy predicts each quadrature without phase unwrapping.
                    predicted = np.exp(
                        1j * np.array(stored["exposures"]) * s * row["exact"]
                    )
                    assert max(abs(complex_signal - predicted)) < 0.0026
                    errors.append(
                        abs(r["energy"] - s * row["exact"]) / r["simultaneous_error"]
                    )
        for r in b["studies"][0]["guide_held"]:
            assert not r["accepted"]
    assert max(errors) < 1
    neutral_rows = json.loads((ex.HERE / "neutral.json").read_text())["rows"]
    assert neutral_rows[-1]["charge_product"] > 15.74
    assert neutral_rows[-2]["charge_product"] < -15.74


def test_packet_trace_in_a_finite_material_parity_sector():
    from . import flight

    m = 5
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    eye = np.eye(2)

    def tensor(words):
        out = np.array([[1]])
        for w in words:
            out = np.kron(out, w)
        return out

    gamma = np.stack(
        [tensor([z] * a + [s] + [eye] * (m - a - 1)) for a in range(m) for s in (x, y)]
    )
    d = 2**m
    idn = np.eye(d)
    parity = (idn + tensor([z] * m)) / 2
    g = np.random.default_rng(2201)
    packet = np.linalg.qr(g.normal(size=(2 * m, 2)))[0]
    recv = np.linalg.qr(g.normal(size=(2 * m, 2)))[0]
    theta = schedule.DELTA
    H1 = sum(tensor([eye] * a + [x, x] + [eye] * (m - a - 2)) for a in range(m - 1))
    H2 = sum(tensor([eye] * a + [z] + [eye] * (m - a - 1)) for a in range(m))
    u = expm(0.5j * theta * H2) @ expm(0.5j * theta * H1)
    o = np.eye(2 * m)
    schedule.rotate_pairs(o, np.arange(1, 2 * m - 1, 2), theta)
    schedule.rotate_pairs(o, np.arange(0, 2 * m, 2), theta)

    def detector(pair):
        a, b = np.einsum("sa,sjk->ajk", pair, gamma)
        return 1j * a @ b

    rho = parity @ (idn + detector(packet)) / (d / 2)
    advanced = np.linalg.matrix_power(u, 7)
    measured = np.trace(detector(recv) @ advanced @ rho @ advanced.conj().T)
    analytic = np.linalg.det(recv.T @ np.linalg.matrix_power(o, 7) @ packet)
    np.testing.assert_allclose(measured, analytic, atol=1e-12)
    recap = json.loads((ex.HERE / "flight.json").read_text())
    assert len([s for s in recap["scans"] if s["qualified"]]) == 2
    assert 2.96 < recap["coefficient"] < 3.04
    assert len(flight.CENTRES) == 361
