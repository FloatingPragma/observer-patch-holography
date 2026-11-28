"""Native operator certificates and energy measurements on pooled ternary preparations."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from model.executable_model import CELL_SEAMS
from .. import source
from ...closed_region.measure import check

HERE = Path(__file__).parent
K = 8


def certificate(k):
    law = source.LAW
    J = law._write_pattern()
    reflection, d2 = law._label_ops()
    swap, diff = law._load_ops()
    twice_g = np.kron(
        np.eye(10, dtype=int) - reflection + 2 * k * d2, np.eye(9, dtype=int)
    )
    twice_b = np.eye(9, dtype=int) - swap + 2 * k * diff
    h4 = source.formation.emission.four_h(k)
    projector2 = J @ J.T
    exact = {
        "write_isometry": J.T @ J - 2 * np.eye(9, dtype=int),
        "memory_intertwines": twice_g @ J - J @ twice_b,
        "written_energy_complement": h4[9:, 9:] - (2 * twice_g - J @ twice_b @ J.T),
        "blank_energy": h4[:9, :9] - 2 * twice_b,
        "free_memory_ruler": h4[9:, 9:] @ twice_g - twice_g @ h4[9:, 9:],
        "empty_write_energy": h4[9:, 9:] @ J,
        "sharp_write_diagonal": np.diag(J.T @ np.diag(np.diag(h4[9:, 9:])) @ J)
        - 2 * np.diag(twice_b),
    }
    assert all(not x.any() for x in exact.values())
    # In the sharp written basis, compatible rows have half the memory energy,
    # every incompatible row has all of it. The local cost depends only on d.
    valid = projector2.diagonal()
    np.testing.assert_array_equal(
        h4[9:, 9:].diagonal(), (2 - valid) * twice_g.diagonal()
    )
    return (
        J,
        twice_g,
        twice_b,
        h4,
        {k: int(np.count_nonzero(x)) for k, x in exact.items()},
    )


def geometry():
    left, right = np.array(CELL_SEAMS).T
    smallest = None
    witness = {}
    for start in range(0, 3**12, 20000):
        digits = (
            np.arange(start, min(start + 20000, 3**12))[:, None]
            // 3 ** np.arange(12)
            % 3
        )
        top = np.sort(digits * digits, axis=1)[:, -3:].sum(axis=1)
        total = (digits * digits).sum(axis=1)
        internal = ((digits[:, left] - digits[:, right]) ** 2).sum(axis=1)
        cost = total - top + internal
        assert np.all(2 * cost >= 3 * top)
        ties = np.flatnonzero((top > 0) & (2 * cost == 3 * top))
        if smallest is None and ties.size:
            j = ties[0]
            witness = dict(
                ports=digits[j].tolist(),
                internal=int(internal[j]),
                degree=int(total[j] - top[j]),
                gate_square=int(top[j]),
            )
            smallest = Fraction(int(cost[j]), int(top[j]))
    assert smallest == Fraction(3, 2)
    return dict(
        words_checked=3**12,
        minimum_complement_per_three_largest=str(smallest),
        witness=witness,
        coherent_external_bound="3/7",
        pooled_classical_sharp_bound="17/82",
        pooled_coherent_input_sharp_bound="17/80",
    )


def trace_energy(h, state):
    return float(np.sum(h.T * state).real)


def scramble(state):
    tensor = state.reshape(10, 9, 10, 9)
    memory = tensor.trace(axis1=1, axis2=3)
    load = tensor.trace(axis1=0, axis2=2)
    return np.kron(memory, load)


def compare(rho, ops, read=False):
    j, g, b, h4, _ = ops
    sigma = j @ rho @ j.T / 2
    blank = trace_energy(b / 2, rho)
    sharp_input = trace_energy(np.diag(np.diag(b)) / 2, rho)
    out = {}
    for name, formed in (("coherent", sigma), ("sharp", np.diag(sigma.diagonal()))):
        prod = scramble(formed)
        h = h4[9:, 9:] / 4
        before = trace_energy(h, formed)
        after = trace_energy(h, prod)
        changed = prod - formed
        delta = changed.reshape(10, 9, 10, 9)
        np.testing.assert_allclose(delta.trace(axis1=0, axis2=2), 0, atol=2e-14)
        np.testing.assert_allclose(delta.trace(axis1=1, axis2=3), 0, atol=2e-14)
        credit = trace_energy(b / 2, j.T @ prod @ j / 2)
        stored = trace_energy(g / 2, prod)
        np.testing.assert_allclose(
            [after + credit, stored],
            [blank if name == "coherent" else sharp_input] * 2,
            atol=1e-13,
        )
        np.testing.assert_allclose(
            before, 0 if name == "coherent" else sharp_input / 2, atol=1e-13
        )
        assert (
            -1e-12
            <= after - before
            <= (blank if name == "coherent" else 17 / 32 * blank) + 1e-12
        )
        out[name] = dict(
            retained=before,
            erased=after,
            gain=after - before,
            unresolved_input_credit=credit,
            memory_energy=stored,
        )
        if read:
            out[name]["banks"] = (formed, prod)
    return blank, sharp_input, out


def experiment():
    ops = certificate(K)
    _, proj = source.energy_projectors(K)
    rng = np.random.default_rng(191201)
    samples = 2**22
    vacuum = np.eye(9)[0]
    symmetric = (np.eye(9)[1] + np.eye(9)[3]) / np.sqrt(2)
    results = []
    for p in (0.5, 1 / 16):
        rho = (1 - p) * np.outer(vacuum, vacuum) + p * np.outer(symmetric, symmetric)
        inc, d, read = compare(rho, ops, read=True)
        for mode, panel in read.items():
            stats = []
            for state in panel.pop("banks"):
                probs = np.array(
                    [trace_energy(P[9:, 9:] / den, state) for e, P, den in proj]
                )
                np.testing.assert_allclose(probs.sum(), 1.0)
                counts = rng.multinomial(samples, probs.clip(0))
                en = np.array([e for e, _, _ in proj])
                stats.append(
                    dict(
                        counts=counts.tolist(),
                        mean=(counts @ en) / samples,
                        forecast=float(probs @ en),
                        probabilities=probs.tolist(),
                    )
                )
            panel["spectral_counts"] = stats
            panel["measured_gain"] = stats[1]["mean"] - stats[0]["mean"]
        results.append(
            dict(bright_probability=p, incident=inc, diagonal_incident=d, readings=read)
        )
    # Sweep complex, generally charge-coherent density matrices. The analytic
    # trace certificate also permits an entangled exterior by using its reduction.
    ranges = {"coherent": [], "sharp": []}
    for rank in (1, 3, 9):
        for _ in range(80):
            z = rng.normal(size=(9, rank)) + 1j * rng.normal(size=(9, rank))
            rho = z @ z.conj().T
            rho = rho / rho.trace()
            inc, _, data = compare(rho, ops)
            for kind in ranges:
                ranges[kind].append(data[kind]["gain"] / inc)
    return dict(
        operators={k: certificate(k)[-1] for k in (0, 1, 8)},
        geometry=geometry(),
        spectra=[e for e, _, _ in proj],
        trials_per_bank=samples,
        shot_panels=results,
        arbitrary_complex_states=240,
        sampled_fraction_ranges={k: [min(v), max(v)] for k, v in ranges.items()},
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    result = json.loads(json.dumps(experiment()))
    path = HERE / "results.json"
    if p.parse_args().write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), result)
