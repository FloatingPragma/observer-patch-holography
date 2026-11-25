"""Exact cyclicity of the measured charged guide and same-pole local receiver readback."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.stats import beta

from . import measure
from .apparatus import Guide, certify_chart, geometry


def certify_cyclic(g):
    """Exact endpoint Krylov columns of 2H in port order, without a dense determinant."""
    d = [int(a) for a in 2 * g.diag]
    b = [int(a) for a in 2 * g.off]
    np.testing.assert_array_equal(d, 2 * g.diag)
    np.testing.assert_array_equal(b, 2 * g.off)
    assert all(x for x in b)
    col = [1] + [0] * (g.n - 1)
    diagonal = 1
    digest = sha256()
    for k in range(g.n):
        assert col[k] == diagonal and not any(col[k + 1 :])
        digest.update((",".join(map(str, col)) + "\n").encode())
        if k < g.n - 1:
            diagonal *= b[k]
        col = [
            d[i] * v
            + (b[i - 1] * col[i - 1] if i else 0)
            + (b[i] * col[i + 1] if i + 1 < g.n else 0)
            for i, v in enumerate(col)
        ]
    return dict(
        ports=g.n,
        nonzero_integer_subdiagonals=len(b),
        diagonal_entries_nonzero=g.n,
        span_dimension=g.n,
        endpoint_columns_sha256=digest.hexdigest(),
        maximum_exact_eigenvalue_multiplicity=1,
    )


def aperture_certificate(g, receivers):
    """Coprime full and prefix continuants certify nonzero weight at a sharp port."""
    prime = 65537
    d = [int(x) for x in 2 * g.diag]
    b = [int(x) for x in 2 * g.off]
    polys = [[1], [(-d[0]) % prime, 1]]
    for i in range(1, g.n):
        y, z = polys[-1], polys[-2]
        item = [0] + y
        for k, val in enumerate(y):
            item[k] = (item[k] - d[i] * val) % prime
        for k, val in enumerate(z):
            item[k] = (item[k] - b[i - 1] ** 2 * val) % prime
        polys.append(item)

    def degree_gcd(f, h):
        while h:
            while len(f) >= len(h):
                shift = len(f) - len(h)
                mult = f[-1] * pow(h[-1], -1, prime)
                for j, y in enumerate(h):
                    f[j + shift] = (f[j + shift] - mult * y) % prime
                while f and f[-1] == 0:
                    f.pop()
            f, h = h, f
        return len(f) - 1

    degrees = [degree_gcd(polys[-1].copy(), polys[r].copy()) for r in receivers]
    assert degrees == [0] * len(receivers)
    return dict(prime=prime, gcd_degrees=degrees, nonzero_weight_every_finite_pole=True)


def captures(g, r0, last, detectors):
    # Emission A is the phased pair of measured modes. Emission B is the same
    # independently sharp charge as used in the reciprocity experiment.
    emitted = g.prepare((g.columns[:, r0] - g.columns[:, last]) / np.sqrt(2))
    sharp = np.eye(1, g.n, detectors[0]).ravel()
    source = np.column_stack([emitted, sharp - emitted * np.vdot(emitted, sharp)])
    source[:, 1] /= np.linalg.norm(source[:, 1])
    np.testing.assert_allclose(source.conj().T @ source, np.eye(2), atol=2e-12)
    rows = []
    for i in measure.selected(g):
        v = g.columns[:, i]
        # P_E has one column direction for the full charged object. A
        # hypothetical coherent reading across both receivers is stronger
        # than either local charge/head reading separately.
        recv = v[detectors]
        launch = v.conj() @ source
        matrix = np.outer(recv, launch)
        sv = np.linalg.svd(matrix, compute_uv=False)
        np.testing.assert_allclose(
            sv, [np.linalg.norm(recv) * np.linalg.norm(launch), 0], atol=1e-13
        )
        assert all(abs(recv) > 1e-12) and sv[0] > 1e-12
        rows.append(
            dict(
                mode=i,
                energy=float(g.levels[i]),
                actual_receiver_weights=(abs(recv) ** 2).tolist(),
                source_receiver_singular_values=[float(sv[0]), 0],
                exact_receiver_dimension=1,
            )
        )
    return rows, emitted


def spin_reference(g, mode, emitted, detectors, rng):
    """Declared charged spin-1/2 reference on exactly the same orbital path and apertures."""
    i2 = np.eye(2, dtype=complex)
    x = np.array([[0, 1], [1, 0]], complex)
    z = np.diag([1, -1]).astype(complex)
    rx, rz = -1j * x, -1j * z
    # Both half-turns preserve the exact integer H x I on every native bond.
    for hentry in np.unique(np.r_[g.diag, g.off]):
        hblock = hentry * i2
        for turn in (rx, rz):
            np.testing.assert_array_equal(turn @ hblock, hblock @ turn)
    np.testing.assert_array_equal(rx @ rz, -rz @ rx)
    np.testing.assert_array_equal(rx @ rz @ rx.conj().T @ rz.conj().T, -i2)
    v = g.columns[:, mode]
    amp = v.conj() @ emitted
    out = []
    for r in detectors:
        # A spin-sensitive analyzer can transmit any relative spin phase at
        # this exact orbital energy and this same sharp charge address.
        k = v[r] * amp * i2
        sv = np.linalg.svd(k, compute_uv=False)
        np.testing.assert_allclose(sv, abs(v[r] * amp) * np.ones(2), atol=1e-13)
        prob = float(sv[0] ** 2)
        assert prob > 1e-5
        banks = []
        for bits, analyzer in (
            ([1, 0], z),
            ([0, 1], z),
            ([2**-0.5, 2**-0.5], x),
            ([2**-0.5, -(2**-0.5)], x),
        ):
            spin = k @ np.array(bits)
            plus = float(np.vdot(spin, (i2 + analyzer) @ spin).real / 2)
            counts = rng.multinomial(
                2**20, [max(0, plus), max(0, prob - plus), 1 - prob]
            )
            banks.append(counts.tolist())
        assert all(sum(b[:2]) > 20 for b in banks)
        assert banks[0][1] == banks[2][1] == banks[1][0] == banks[3][0] == 0
        lower = float(
            beta.ppf(
                0.01 / 16,
                min(sum(b[:2]) for b in banks),
                2**20 + 1 - min(sum(b[:2]) for b in banks),
            )
        )
        out.append(
            dict(
                index=r,
                success_probability=prob,
                spin_sensitive_singular_values=sv.tolist(),
                charge_only_probability_both_spins=prob,
                starts_per_test=2**20,
                simultaneous99_success_lower=lower,
                z_plus_minus_x_plus_minus=banks,
            )
        )
    # A deliberately split known spin is a counterexample to any inference
    # from a single resolved line to an absence of spin in a broader theory.
    row_bound = float(np.max(abs(g.diag) + np.r_[0, g.w] + np.r_[g.w, 0]))
    splitting = 64
    assert splitting > 2 * row_bound
    p_up = np.diag([1, 0])
    np.testing.assert_array_equal(p_up @ x @ p_up, np.zeros((2, 2)))
    return dict(
        model="supplied H tensor identity on spin",
        mode=mode,
        known_split_control=dict(
            shift=splitting,
            orbital_norm_bound=row_bound,
            individual_fixed_energy_spin_dimension=1,
        ),
        half_turn_commutator_sign=-1,
        analyzers=out,
    )


def run():
    src = json.loads((measure.HERE / "records.json").read_text())
    ports, board = geometry()
    cases = []
    for data in src["observations"]:
        n, p, sign = data["n"], data["period"], data["charge"]
        g = Guide(n, p)
        chart = certify_chart(ports[:n], board, g, sign)
        assert chart == data["charged_chart"]
        apertures = data["source_receiver"]["receivers"]
        indices = [r["index"] for r in apertures]
        assert [int(ports[i]) for i in indices] == [r["port"] for r in apertures]
        reference = data["reference_index"]
        rows, source = captures(g, reference, data["source_receiver"]["mode"], indices)
        cases.append(
            dict(
                length=n,
                charge=sign,
                actual_detector_ports=[int(ports[i]) for i in indices],
                cyclicity=certify_cyclic(g),
                aperture_certificate=aperture_certificate(g, indices),
                pole_readouts=rows,
                declared_spin_control=spin_reference(
                    g,
                    reference,
                    source,
                    indices,
                    np.random.default_rng(n + sign + 37681),
                ),
            )
        )
    return dict(
        registered_law="Interval Echo",
        formation_source="dispersion_search/records.json",
        sign_independent_commuting_centralizer=True,
        fixtures=cases,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2) + "\n"
    if ap.parse_args().write:
        Path(__file__).with_name("spin_fibers.json").write_text(out)
    else:
        print(out, end="")
