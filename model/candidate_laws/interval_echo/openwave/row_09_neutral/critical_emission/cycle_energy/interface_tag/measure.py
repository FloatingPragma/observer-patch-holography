"""Compose two local charge records with the original distal vacuum tag."""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import cos, pi, sin, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ........executable_model.carrier import federation
from .......contrast_memory_spherical_v1 import law
from .....row_20_electric import experiment as prior
from ....passive_decay.measure import compare_records
from ...apparatus import Guide, aperture
from ...emitter import absorber, three_exits
from ...release import enroll
from ..measure import DELTA
from model.experiments.interval_echo_consensus.signed_effects import (
    fixed_pointer as pointer,
)
from model.experiments.material_inertia.pair_spectrum.native_source_probe import (
    physical,
)
from model.experiments.material_inertia.pair_spectrum.native_source_probe.measure import (
    routes,
)


HERE = Path(__file__).resolve().parent
SHOTS = 2**42
SEED = 90133
LABELS = (
    "parent_refused",
    "mode_refused",
    "charge_integrity_refused",
    "minus_two_tag",
    "minus_two_other",
    "zero_tag",
    "zero_other",
    "plus_two_tag",
    "plus_two_other",
)


def pointer_columns():
    # Exact relative phase on each data word: V(a,2)-V(a,0)=65(1-a)/2.
    def potential(a, b):
        return 8 * (a - b) ** 2 + Fraction(int(a != b), 2)

    assert all(
        potential(a, 2) - potential(a, 0) == Fraction(65, 2) * (1 - a) for a in range(3)
    )
    ids = np.array(
        [
            [
                int(np.array((a, b, 1, 2, 0)) @ [81, 27, 9, 3, 1]),
                int(np.array((a, b, 1, 0, 2)) @ [81, 27, 9, 3, 1]),
            ]
            for a, b in product(range(3), repeat=2)
        ]
    )
    columns = np.eye(243)[:, ids[:, 0]]
    basis = pointer.neutral_basis(3, 4)
    probabilities = {}
    for name, ratio, sine in (
        ("odd_magnitude", Fraction(2, 65), False),
        ("odd_sign", Fraction(1, 65), True),
        ("even_magnitude", Fraction(1, 65), False),
        ("even_sign", Fraction(1, 130), True),
    ):
        word = (
            basis
            + pointer.diagonal_word(0, 4, ratio)
            + pointer.diagonal_word(1, 4, ratio)
        )
        if sine:
            word += [pointer.pulse("reference", 4, 4, Fraction(1, 64))]
        out = physical.primitive(word + pointer.inverse(basis), columns)
        pair = out[ids, np.arange(9)[:, None]]
        projected = np.zeros_like(out)
        projected[ids, np.arange(9)[:, None]] = pair
        np.testing.assert_allclose(out, projected, atol=8e-13, rtol=0)
        np.testing.assert_allclose(pair[2], pair[6], atol=8e-13, rtol=0)
        expected = [
            (
                1
                + (sin if sine else cos)(
                    pi * float(ratio * Fraction(65, 2)) * (a + b - 2)
                )
            )
            / 2
            for a, b in product(range(3), repeat=2)
        ]
        np.testing.assert_allclose(abs(pair[:, 0]) ** 2, expected, atol=8e-13, rtol=0)
        probabilities[name] = expected
    for i, (a, b) in enumerate(product(range(3), repeat=2)):
        charge = a + b - 2
        first = probabilities["odd_magnitude"][i]
        assert first == float(charge % 2 == 0)
        if first:
            quiet = probabilities["even_magnitude"][i]
            decoded = 0 if quiet else 2 if probabilities["even_sign"][i] else -2
        else:
            decoded = 1 if probabilities["odd_sign"][i] else -1
        assert decoded == charge
    return probabilities


def routing(board, ports):
    protected = {p // 12 for p in ports}
    result = []
    for pair in (ports[-4:-2], ports[-2:]):
        m = routes(board, protected)
        source = [m["source_park"], m["source_cross"]]
        paths = [
            m["endpoint_transposition_paths"][m["target_ports"].index(p)] for p in pair
        ]
        rounds = [2 * physical.validate_route(board, path, *source) for path in paths]
        result.append(
            {
                "data_ports": pair,
                "pointer_ports": source,
                "routes": paths,
                "round_trip_swaps": rounds,
            }
        )
        protected.add(source[0] // 12)
    return result


def unlike_column(g):
    # gamma_(2j+1) gamma_(2j+2) represents the two neighboring rail signs.
    # D=(I+Z_j Z_(j+1))/2 on the encoded dimers. Its one-mode column
    # follows from contracting i gamma_a gamma_b between creator states.
    a, b = 2 * g.n - 3, 2 * g.n - 2
    w = g.modes
    result = 0.5j * (w[a, :].conj() * w[b, 0] - w[b, :].conj() * w[a, 0])
    result[0] += (1 + g.c[a, b]) / 2
    np.testing.assert_allclose(result[0].imag, 0, atol=1e-14, rtol=0)
    return result


def independent_full_words():
    """Check the charge branches on four dimers using full rail words."""
    n = 4
    g = Guide(n, "native", 0.5)
    identity = np.eye(2)
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    gamma = []
    for j, spin in product(range(n), (z, y)):
        term = np.ones((1, 1))
        for p in range(n):
            term = np.kron(term, x if p < j else spin if p == j else identity)
        gamma.append(term)
    gamma = np.asarray(gamma)
    destroy = np.einsum("ri,rxy->ixy", g.modes.conj() / 2, gamma)
    vacuum_projector = np.eye(2**n, dtype=complex)
    for gate in destroy:
        vacuum_projector = vacuum_projector @ gate @ gate.conj().T
    j = int(np.argmax(vacuum_projector.diagonal().real))
    vac = vacuum_projector[:, j] / sqrt(vacuum_projector[j, j].real)
    frame = np.column_stack([a.conj().T @ vac for a in destroy])
    psi = frame[:, 0]
    beta = np.array([0.31 + 0.11j, -0.17j, 0.24, 0.08 - 0.21j])
    tag = beta @ frame.conj().T
    orientations = list(product((0, 1), repeat=n))
    rails = [tuple(u for s in row for u in (s, 1 - s)) for row in orientations]
    index = {row: j for j, row in enumerate(product((0, 1), repeat=2 * n))}
    dest = [index[row] for row in rails]
    tag_full = np.zeros(2 ** (2 * n), complex)
    tag_full[dest] = tag
    changed = np.zeros_like(tag_full)
    s, c = sin(DELTA / 65), cos(DELTA / 65)
    swap = (*range(2 * n - 3), 2 * n - 2, 2 * n - 3, 2 * n - 1)
    for i, row in enumerate(rails):
        if row[-3] != row[-2]:
            changed[dest[i]] += c * psi[i]
            changed[index[tuple(row[k] for k in swap)]] += 1j * s * psi[i]
        else:
            changed[dest[i]] += psi[i]
    projection = beta @ unlike_column(g)
    charged_weights = []
    for q in (-2, 0, 2):
        mask = np.array([2 * (sum(row[-2:]) - 1) == q for row in index])
        word = changed * mask
        direct = tag_full @ word
        if q:
            reversed_word = word[[index[tuple(row[k] for k in swap)] for row in index]]
            restored = tag_full @ reversed_word
            np.testing.assert_allclose(direct, 0, atol=2e-14, rtol=0)
            np.testing.assert_allclose(
                restored, 1j * s * projection / 2, atol=2e-14, rtol=0
            )
            charged_weights.append(float(np.vdot(word, word).real))
        else:
            np.testing.assert_allclose(
                direct, beta[0] - (1 - c) * projection, atol=2e-14, rtol=0
            )
    np.testing.assert_allclose(
        charged_weights, [s**2 * unlike_column(g)[0].real / 2] * 2, atol=2e-14, rtol=0
    )
    return {"logical_dimers": n, "full_extreme_rail_words": len(index)}


def panel(n, archive, emitter_row, board, path, rng):
    g = Guide(n, "native", 0.5)
    one, admission = enroll(g)
    target = emitter_row["target"]
    beta = g.receiver * aperture((target * g.period) - g.phi, 2048)
    nominal = abs(beta[0]) ** 2
    np.testing.assert_allclose(
        nominal, three_exits(one, absorber(g, target), False)[0], atol=3e-13, rtol=0
    )
    d = unlike_column(g)
    b = beta @ d
    s, c = sin(DELTA / 65), cos(DELTA / 65)
    indices = [law.SRC.index(w) for w in ((0, 2), (2, 0))]
    h = law.seam_term(8)
    # Removing only the guard replaces the diagonal phase on 02/20 by
    # that same phase times c I + i s SWAP.
    native = expm(-1j * (2 * DELTA / 65) * h)[np.ix_(indices, indices)]
    expected = np.exp(-1j * (2 * DELTA / 65) * 32.5) * np.array(
        [[c, 1j * s], [1j * s, c]]
    )
    np.testing.assert_allclose(native, expected, atol=2e-15, rtol=0)
    swap = np.eye(9)[[law.SRC.index((y, x)) for x, y in law.SRC]]
    np.testing.assert_allclose(expm(-1j * pi * h), swap, atol=4e-14, rtol=0)
    p_signed = s * s * d[0].real / 2
    q0tag = abs(beta[0] - (1 - c) * b) ** 2
    charged_return = s * s * abs(b / 2) ** 2
    old = archive["guard_only_terminal_charge_control"]
    np.testing.assert_allclose(
        2 * p_signed, old["motion_per_admitted_exposed_dimer"], rtol=1e-12, atol=0
    )
    np.testing.assert_allclose(
        admission, archive["parent_refused_mode_refused_admitted"], rtol=0, atol=1e-13
    )
    assert abs(b) > 0.02 and charged_return > 0 and 0 < q0tag < 1 - 2 * p_signed
    ready = float(admission[2])
    # A weak independent bound from the archived bright row suffices to
    # certify a positive tag even if its remaining phases are adversarial.
    # The inherited mode certificate bounds eigenphase error below 2^-38.
    # Four 2048-label delays give an aperture derivative at most 8188.
    tol = archive["conservative_covariance_allowance"] + 3e-8
    bright = sqrt(emitter_row["bright_only_tag_upper"])
    tagged_lower = (
        (ready - tol)
        * s
        * s
        * max(
            0,
            (sqrt(nominal) - tol) * (d[0].real - tol) / 2
            - (bright + tol) * sqrt((d[0].real + tol) / 2),
        )
        ** 2
    )
    assert 0 < tagged_lower <= ready * charged_return
    histories = {}
    for name, tag in (("direct", 0), ("charge_conditioned_return", charged_return)):
        inclusive = ready * np.array(
            [tag, p_signed - tag, q0tag, 1 - 2 * p_signed - q0tag, tag, p_signed - tag]
        )
        raw = np.r_[admission[:2], 0.0, inclusive]
        assert np.min(raw) >= 0 and abs(raw.sum() - 1) < 8e-15
        counts = rng.multinomial(SHOTS, raw)
        assert counts.sum() == SHOTS
        assert (
            (counts[3] > 0 and counts[7] > 0) if tag else (counts[3] == counts[7] == 0)
        )
        histories[name] = {
            "probabilities": raw.tolist(),
            "seeded_draws": counts.tolist(),
        }
    return {
        "dimers": n,
        "charge_meter_routes": routing(board, path[1 : 1 + 2 * n]),
        "original_low_mode_vacuum_tag_probability_given_admission": float(nominal),
        "original_low_mode_amplitude_real_imaginary": [beta[0].real, beta[0].imag],
        "one_mode_matrix_element_of_total_unlike_projector": d[0].real,
        "guard_and_aperture_allowance": tol,
        "coherent_end_contraction_real_imaginary": [b.real, b.imag],
        "charged_tag_probability_per_original_lower_after_deliberate_return": float(
            tagged_lower
        ),
        "guarded_and_return_arms": histories,
    }


def run():
    a = HERE.parent / "results.json"
    e = HERE.parents[1] / "emitter.json"
    old = json.loads(a.read_text())
    emitter = json.loads(e.read_text())
    rng = np.random.default_rng(SEED)
    board = federation(2)
    _, path, _ = prior.guide(board)
    results = []
    for p in old["panels"]:
        if p["schedule"] != "native":
            continue
        em = next(
            x
            for x in emitter["panels"]
            if x["dimers"] == p["dimers"] and x["source_link_fraction"] == 0.5
        )
        results.append(panel(p["dimers"], p, em, board, path, rng))
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "cycle_energy_sha256": sha256(a.read_bytes()).hexdigest(),
        "emission_sha256": sha256(e.read_bytes()).hexdigest(),
        "calibration_input_word_order": [
            str(a) + str(b) for a, b in product(range(3), repeat=2)
        ],
        "native_meter_head20_probabilities": pointer_columns(),
        "full_rail_oracle": independent_full_words(),
        "shots_per_arm": SHOTS,
        "seed": SEED,
        "original_outcome_order": list(LABELS),
        "receiving_and_neighboring_charge_per_valid_record": [[-2, 2], [0, 0], [2, -2]],
        "cases": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    if args.acquire:
        with args.acquire.open("x") as out:
            out.write(json.dumps(data, indent=1) + "\n")
    else:
        compare_records(json.loads((HERE / "results.json").read_text()), data)
        print(
            "Verified both pre-absorption charge records, the direct tag and the charged-arm return."
        )


if __name__ == "__main__":
    main()
