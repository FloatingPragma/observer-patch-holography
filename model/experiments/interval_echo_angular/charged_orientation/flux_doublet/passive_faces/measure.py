"""Repeated face transport as a gauge-invariant detector on the guarded charge shell."""

from fractions import Fraction
from functools import cache
import json
from math import prod
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from model.experiments.interval_echo_lepton_register.measure import route
from ...live_trimer.measure import LAW, candidate
from ..measure import CELL, faces, face_edges, links, spatial_order
from ..passive import exact_columns

HERE = Path(__file__).resolve().parent
SHOTS = 2**18
SPECTRUM = (0, 8, 9, 32, 33)


@cache
def twice_native_crossings():
    """Return exact integer 2K for each charge sign and edge orientation."""
    write = LAW._write_pattern()
    q = 2 * np.eye(90, dtype=np.int64) - write @ write.T
    flip, d2 = LAW._label_ops()
    x, c = LAW._load_ops()
    k = candidate.KAPPA
    native = np.zeros((99, 99), dtype=np.int64)
    native[:9, :9] = 4 * (np.eye(9, dtype=np.int64) - x) + 8 * k * c
    generator = np.kron(
        np.eye(10, dtype=np.int64) - flip + 2 * k * d2, np.eye(9, dtype=np.int64)
    )
    assert not np.any(q @ generator @ write)
    native[9:, 9:] = q @ generator @ q
    ex, co = LAW.local_term_parts()
    np.testing.assert_allclose(native, 8 * (ex + k * co), atol=1e-13)
    eye = np.eye(99, dtype=np.int64)
    factors = [native - 8 * e * eye for e in SPECTRUM]
    annihilator = eye.copy()
    for f in factors:
        annihilator = annihilator @ f
    assert not np.any(annihilator)
    # Lagrange interpolation of exp(-i*pi*H) on the certified integer spectrum.
    # Products for each numerator have norm at most 264^4; all fit in int64.
    terms = []
    for i, e in enumerate(SPECTRUM):
        value = eye.copy()
        for j, f in enumerate(factors):
            if j != i:
                value = value @ f
        den = prod(8 * (e - y) for y in SPECTRUM if e != y)
        terms.append(((-1) ** e * value, den))
    direct = expm(-1j * np.pi * (ex + k * co))
    result = {}
    for sign in (-1, 1):
        for reverse in (False, True):
            start = (1, 1 + sign) if reverse else (1 + sign, 1)
            i, j = map(LAW.SRC.index, (start, start[::-1]))
            indices = np.ix_(np.arange(j, 99, 9), np.arange(i, 99, 9))
            projected = [(a[indices], b) for a, b in terms]
            entries = [
                [
                    sum(Fraction(2 * int(a[s, t]), b) for a, b in projected)
                    for t in range(11)
                ]
                for s in range(11)
            ]
            assert all(v.denominator == 1 for row in entries for v in row)
            m = np.array(entries, dtype=np.int64)
            np.testing.assert_allclose(m / 2, direct[indices], atol=5e-13)
            blank = np.zeros((11, 11), dtype=np.int64)
            blank[0, 0] = 1
            np.testing.assert_array_equal(m @ m, 4 * blank)
            hop = native[indices]
            assert hop[0, 0] == -4
            # H is stored in eighths. Two directed hops on the same head
            # coordinate have square (1/4) times its ordinary-blank projector.
            np.testing.assert_array_equal(hop @ hop, 16 * blank)
            np.testing.assert_array_equal(hop @ blank, -4 * blank)
            np.testing.assert_array_equal(blank @ hop, -4 * blank)
            assert not np.any(m[1:, 0]) and not np.any(m[0, 1:])
            result[sign, reverse] = m
    return result


def controls():
    # Conservation of the material head current also covers non-pi waits.
    columns = exact_columns()
    rows = []
    for q in range(20):
        root = (7 * q + 3) % 20
        field = links(q, root=root)
        shift = np.array([(3 * v * v + 2 * q * v) % 20 for v in range(12)])
        gauge = field + shift[None, :] - shift[:, None]
        signs = spatial_order(field)
        assert spatial_order(gauge) == signs
        assert signs == [10 * (q % 2)] * 15
        residue = (2 * q) % 20
        matches_native = residue == 0
        rows.append(
            dict(
                flux=q,
                twice_face_residue=residue,
                negative_half_turn_pairs=15 * (q % 2),
                unit_visibility_two_lap_match=matches_native,
            )
        )
    return columns, rows


def paths():
    shell = faces()
    assert len(shell) == 20
    edges = {tuple(sorted(edge)) for face in shell for edge in face_edges(face)}
    assert edges == set(CELL.seams)
    assert all(
        sum(e in [tuple(sorted(k)) for k in face_edges(f)] for f in shell) == 2
        for e in edges
    )
    # A native dressed source containing a single written odd-parity seam
    # has its two possible directions in a complex superposition.
    return shell, (0, 1)


def write_state(sign, forward, one_lap_phase):
    c = int((sign > 0) ^ (not forward))

    def vector(bit):
        state = np.zeros(11, complex)
        for d in (-1, 1):
            state[1 + LAW.LABELS.index((bit, d))] = d / np.sqrt(2)
        return state

    # Sharp-charge controlled lap and erasure of its balanced path pointer.
    # The two blank links of the triangle give a spatial return with +1.
    initial = vector(c)
    flipped = twice_native_crossings()[sign, not forward] @ initial / 2
    np.testing.assert_allclose(flipped, vector(1 - c), atol=1e-13)
    selected = (initial + np.exp(-1j * one_lap_phase) * flipped) / 2
    np.testing.assert_allclose(np.vdot(selected, selected).real, 0.5)
    return selected * np.sqrt(2)


def quadratures(expected, norm, rng, source_admission=1.0):
    # The receiving beam splitter compares a coherently heralded path with
    # an idle reference. The route with a failed charge hop has a separate bin.
    banks = []
    for phase in (0.0, np.pi / 2):
        v = float((np.exp(-1j * phase) * expected).real)
        p = np.maximum(
            [(1 + norm + 2 * v) / 4, (1 + norm - 2 * v) / 4, (1 - norm) / 2], 0.0
        )
        assert abs(sum(p) - 1) < 1e-12
        p = np.r_[source_admission * p, 1 - source_admission]
        banks.append(rng.multinomial(SHOTS, p / sum(p)).tolist())
    return banks


def supplied_read(sign, rng):
    # A field pulse is measured against the ordinary full edge exchange on
    # both orientations. Infer the site corrections from these separate data.
    field = links(-sign)
    amps = {}
    measured = np.zeros((12, 12), int)
    ex, ct = LAW.local_term_parts()
    states = [LAW.SRC.index((1 + sign, 1)), LAW.SRC.index((1, 1 + sign))]
    contact = expm(-1j * np.pi * (ex + candidate.KAPPA * ct)[np.ix_(states, states)])
    private = ct.diagonal()
    delta = candidate.KAPPA * (
        private[LAW.SRC.index((1 + sign, 0))] - private[LAW.SRC.index((1, 0))]
    )
    np.testing.assert_allclose(delta, -8 if sign < 0 else 24)
    board = federation(0)

    def hold(v, port, angle):
        trail = min(
            (route(port, j) for j in range(12) if board.degree[j] == 5), key=len
        )
        pairs = list(zip(trail, trail[1:]))
        for a, b in pairs:
            v[[a, b]] = contact @ v[[a, b]]
        dwell = (-angle * np.sign(delta) % (2 * np.pi)) / abs(delta)
        v[trail[-1]] *= np.exp(-1j * delta * dwell)
        for a, b in reversed(pairs):
            v[[a, b]] = contact @ v[[a, b]]

    samples = []
    for a, b in CELL.seams:
        amp = 2 * np.pi * field[a, b] / 20
        v = np.eye(12, dtype=complex)
        hold(v, b, -amp)
        v[[a, b]] = contact @ v[[a, b]]
        hold(v, b, amp)
        for start, finish in ((a, b), (b, a)):
            z = v[finish, start]
            assert abs(abs(z) - 1) < 2e-13
            amps[start, finish] = z
        banks = quadratures(amps[a, b], 1, rng)
        means = [(row[0] - row[1]) / SHOTS for row in banks]
        angle = np.arctan2(means[1], means[0])
        n = int(np.rint(20 * angle / (2 * np.pi)))
        assert (n - field[a, b]) % 20 == 0
        measured[a, b] = n
        measured[b, a] = -n
        samples.append(dict(link=[a, b], counts=banks, residue=n % 20))
    residues = [
        int(sum(measured[i, j] for i, j in face_edges(f)) % 20) for f in faces()
    ]
    assert residues == [(-sign) % 20] * 20
    assert spatial_order(measured) == [10] * 15
    return amps, dict(
        charge=sign,
        private_relative_cost=float(delta),
        measured_face_residues=residues,
        negative_order_pairs=15,
        links=samples,
    )


def run():
    transfers = twice_native_crossings()
    columns, classes = controls()
    shell, special = paths()
    rng = np.random.default_rng(909231)
    calibrations = []
    records = []
    basis = np.eye(11)[0].astype(complex)
    for sign in (-1, 1):
        control, calibration = supplied_read(sign, rng)
        calibrations.append(calibration)
        for face in shell:
            # Each trial begins with a sharp local charge at this face.
            # Keep the same material input for all measurements on the shell.
            for stock in ("blank", "one_chiral_head", "supplied_unit_flux"):
                mean = []
                norm = []
                observations = []
                for laps in (1, 2):
                    overlap = 1.0 + 0j
                    prob = 1.0
                    for a, b in face_edges(face):
                        reverse = a > b
                        op = transfers[sign, reverse] / 2
                        phase = 1.0
                        if stock == "supplied_unit_flux":
                            phase = control[a, b]
                        state = basis
                        if (
                            stock == "one_chiral_head"
                            and tuple(sorted((a, b))) == special
                        ):
                            # Orientation convention always follows sorted (0,1).
                            state = write_state(sign, True, np.pi / 2)
                        final = np.linalg.matrix_power(op * phase, laps) @ state
                        overlap *= np.vdot(state, final)
                        prob *= np.vdot(final, final).real
                    mean.append([float(overlap.real), float(overlap.imag)])
                    norm.append(float(prob))
                    observations.append(
                        quadratures(
                            overlap,
                            prob,
                            rng,
                            0.5 if stock == "one_chiral_head" else 1.0,
                        )
                    )
                records.append(
                    dict(
                        charge=sign,
                        face=list(face),
                        source=stock,
                        complex_idle_overlap=mean,
                        route_survival=norm,
                        sampled_X_Y_by_lap=observations,
                    )
                )
    return dict(
        native_columns=columns,
        charge_signs=[-1, 1],
        classes=classes,
        min_complex_distance_to_passive_for_odd_flux=float(np.sin(np.pi / 5)),
        shots_per_axis=SHOTS,
        seed=909231,
        count_order=["bright", "dark", "route_refusal", "source_refusal"],
        independent_calibrations=calibrations,
        phase_control="independently supplied unit outward flux for vacancy, opposite for excess",
        records=records,
    )


if __name__ == "__main__":
    result = run()
    # Human-readable metadata with one row per bank preserves compact receipts.
    arrays = {
        k: result.pop(k) for k in ("classes", "independent_calibrations", "records")
    }
    body = json.dumps(result, indent=2)[:-2]
    for key, values in arrays.items():
        body += (
            ",\n  "
            + json.dumps(key)
            + ": [\n    "
            + ",\n    ".join(json.dumps(v) for v in values)
            + "\n  ]"
        )
    (HERE / "results.json").write_text(body + "\n}\n")
