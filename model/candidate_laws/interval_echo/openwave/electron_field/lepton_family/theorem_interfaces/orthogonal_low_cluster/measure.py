"""Certify corner spin insertion and extraction for the isolated cap-three band."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path

import sympy as s

from ..charged_endpoints.measure import SOURCE_REVISION, read_archive


HERE = Path(__file__).resolve().parent
MINMAX = "model/native_charged_cluster_minmax/results.json"
BRIDGE = (
    "model/field_record_addenda/native_orthogonal_wilson_bridge/research_reference.json"
)
BRIDGE_ALGEBRA = (
    "model/field_record_addenda/native_orthogonal_wilson_bridge/retained/algebra.py.raw"
)
CHARGED_BAND = (
    "model/field_record_addenda/native_orthogonal_wilson_charged_band/results.json"
)
NATIVE = "model/native_wilson_maxwell/results.json"
SCALE = 2**52


def root_lo(x):
    assert x >= 0
    return F(isqrt(x.numerator * SCALE**2 // x.denominator), SCALE)


def root_hi(x):
    return root_lo(x) + F(1, SCALE)


def exact_matrix_records():
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    kron = s.kronecker_product
    alpha = [kron(pauli[0], p) for p in pauli]
    beta = s.diag(1, 1, -1, -1)
    gamma = alpha[0] * alpha[1] * alpha[2] * beta
    h = kron(s.eye(8), beta)
    for j in range(3):
        x = kron(*[pauli[0] if k == j else s.eye(2) for k in range(3)])
        y = kron(*[pauli[1] if k == j else s.eye(2) for k in range(3)])
        h += kron((s.eye(8) - x) / 2, gamma) + kron(y / 2, alpha[j])

    identity = s.eye(32)
    z, z0 = h * h, s.Rational(5, 4)

    def r(matrix):
        return (
            (4 * matrix - 13 * identity)
            * (16 * matrix**2 - 136 * matrix + 241 * identity)
            * (16 * matrix**2 - 264 * matrix + 977 * identity)
        )

    squares = (r(z) / r(z0 * identity)[0, 0]).applyfunc(s.expand)
    free = (squares / 2 + squares * h / s.sqrt(5)).applyfunc(s.expand)
    assert (free**2 - free).applyfunc(s.expand) == s.zeros(32)
    assert (h * free - s.sqrt(5) * free / 2).applyfunc(s.expand) == s.zeros(32)
    assert free.H == free and s.trace(free) == 2

    vertices = tuple(product((0, 1), repeat=3))

    def rotate(v):
        return v[2], v[0], v[1]

    vertex = s.zeros(8)
    for j, v in enumerate(vertices):
        vertex[vertices.index(rotate(v)), j] = 1
    axis = sum(pauli, s.zeros(2)) / s.sqrt(3)
    spin = (s.eye(2) - s.I * s.sqrt(3) * axis) / 2
    rotation = kron(vertex, s.diag(spin, spin))
    assert (rotation**3 + identity).applyfunc(s.expand) == s.zeros(32)
    assert (rotation * h * rotation.H - h).applyfunc(s.expand) == s.zeros(32)
    assert s.simplify(s.trace(rotation * free)) == 1
    assert s.simplify(s.trace(rotation**2 * free)) == -1
    assert vertex[:, 0] == s.eye(8)[:, 0] and vertex[:, 7] == s.eye(8)[:, 7]
    assert (spin.H * spin).applyfunc(s.expand) == s.eye(2)

    # C3 acts on the physical cochains by preserving the incidence and the
    # plaquette quadratic form; the same lift preserves the total photon cap.
    edges = [
        (v, tuple(w + int(k == j) for k, w in enumerate(v)))
        for v in vertices
        for j in range(3)
        if v[j] == 0
    ]
    d = s.zeros(12, 8)
    permutation = s.zeros(12)
    for row, (v, w) in enumerate(edges):
        d[row, vertices.index(v)] = -1
        d[row, vertices.index(w)] = 1
        permutation[edges.index((rotate(v), rotate(w))), row] = 1
    assert permutation * d == d * vertex
    c = s.zeros(6, 12)
    row = 0
    for v in vertices:
        for a, b in ((0, 1), (0, 2), (1, 2)):
            if v[a] or v[b]:
                continue
            first = tuple(w + int(k == a) for k, w in enumerate(v))
            second = tuple(w + int(k == b) for k, w in enumerate(v))
            end = tuple(w + int(k == b) for k, w in enumerate(first))
            for edge, sign in (
                ((v, first), 1),
                ((first, end), 1),
                ((second, end), -1),
                ((v, second), -1),
            ):
                c[row, edges.index(edge)] = sign
            row += 1
    assert row == 6 and c * d == s.zeros(6, 8)
    k2 = c.T * c
    assert permutation * k2 * permutation.T == k2

    local = {}
    for site in (0, 7):
        at = 4 * site
        b = free[at : at + 2, at : at + 2]
        for sign in (-1, 1):
            q = (s.eye(2) + sign * axis) / 2
            assert (q**2 - q).applyfunc(s.simplify) == s.zeros(2)
            assert (q * b - b * q).applyfunc(s.simplify) == s.zeros(2)
            residue = s.simplify(s.trace(q * b))
            weight_sign = sign * (1 if site == 0 else -1)
            assert (
                s.simplify(
                    residue
                    - (
                        s.Rational(1, 16)
                        + s.sqrt(5) / 40
                        + weight_sign * (s.sqrt(3) / 32 + s.sqrt(15) / 80)
                    )
                )
                == 0
            )
            local[(site, sign)] = str(residue)
    assert (
        free[28:30, 0:2] - (s.Rational(1, 32) + s.sqrt(5) / 80) * s.eye(2)
    ).applyfunc(s.simplify) == s.zeros(2)
    conjugation = s.I * beta * alpha[1]
    assert conjugation * beta.conjugate() * conjugation.H == -beta
    assert (
        conjugation * s.diag(spin, spin).conjugate() * conjugation.H
        - s.diag(spin, spin)
    ).applyfunc(s.simplify) == s.zeros(4)
    return {
        "free_corner_beta_weights": [
            {
                "site_index": site,
                "relative_spin_axis_sign": chi,
                "exact_weight": local[(site, chi)],
            }
            for site in (0, 7)
            for chi in (-1, 1)
        ],
        "free_root_to_opposite_compressed_projector_scalar": "1/32 + sqrt(5)/80",
        "corner_rotation_fixes_actual_root_and_opposite_site": True,
        "spin_lift_cubed": "-I",
        "physical_incidence_and_photon_quadratic_covariance_checked": True,
        "local_charge_conjugate_hole_orbital": "i beta alpha_2 conjugate(v)",
    }


def atan_small(x, n=48):
    assert 0 <= x < F(1, 2)
    t = x
    accum = F(0)
    for j in range(n):
        accum += t / (2 * j + 1)
        t *= -x * x
    error = abs(t) / (2 * n + 1)
    return accum - error, accum + error


def phases(x):
    p, q = atan_small(F(1, 5)), atan_small(F(1, 239))
    pi0, pi1 = 16 * p[0] - 4 * q[1], 16 * p[1] - 4 * q[0]
    if x < F(1, 2):
        a, b = atan_small(x)
        return a / pi1, b / pi0
    assert x >= 1
    a, b = atan_small((x - 1) / (x + 1))
    return F(1, 4) + a / pi1, F(1, 4) + b / pi0


def angle_error(eps):
    return root_hi(2 - 2 * root_lo(1 - eps**2))


def phase_error(size, delta):
    x = F(1, size) / delta
    return (x + x * x) / 2


def run():
    source = {
        name: read_archive(name)
        for name in (MINMAX, BRIDGE, BRIDGE_ALGEBRA, CHARGED_BAND, NATIVE)
    }
    intervals = json.loads(source[MINMAX])
    bridge = json.loads(source[BRIDGE])
    assert sha256(source[BRIDGE]).hexdigest() == intervals["frozen_bridge_receipt"]
    assert (
        sha256(source[BRIDGE_ALGEBRA]).hexdigest()
        == bridge["owned_source_sha256"]["algebra.py"]
    )
    assert F(
        bridge["energy_alignment"]["robust_declared_positive_complement_coefficient"]
    ) == F(21, 2)
    assert F(bridge["energy_alignment"]["g"]) == F(1, 64)
    native = json.loads(source[NATIVE])
    pcode = F(native["actual_preparation"]["global_code_probability_lower"])
    exact = exact_matrix_records()
    energy = intervals["conservative_actual_energy"]
    n0lo, n0hi = [
        F(energy[f"neutral_code_shift_{d}_over_m"]) for d in ("lower", "upper")
    ]
    b1lo, b1hi = [
        F(energy[f"charged_reference_corrected_shift_{d}_over_m"])
        for d in ("lower", "upper")
    ]
    gap = F(intervals["minmax"]["full_charged_rank2_gap_by_completion"]["21/2"])
    code_gap = F(intervals["minmax"]["charged_code_gap_lower_over_m"])
    charged_band = json.loads(source[CHARGED_BAND])
    assert (
        charged_band["source"]["frozen_bridge_sha256"]
        == intervals["frozen_bridge_receipt"]
    )
    assert charged_band["source"]["neutral_active_ground_unchanged_from_bridge"]
    nnu = F(charged_band["bounds_exact"]["neutral_perturbation_norm_over_m_upper"])
    e0 = root_hi(F(energy["neutral_trial_residual_squared_upper"])) / (2 - nnu)
    e1 = root_hi(F(energy["charged_trial_residual_squared_upper"])) / code_gap
    eta = angle_error(e0) + angle_error(e1)
    w = {
        -1: F(1, 16) + root_lo(F(5)) / 40 - root_hi(F(3)) / 32 - root_hi(F(15)) / 80,
        1: F(1, 16) + root_lo(F(5)) / 40 + root_lo(F(3)) / 32 + root_lo(F(15)) / 80,
    }
    x = {chi: root_lo(weight) for chi, weight in w.items()}
    assert min(x.values()) > eta
    residue = {chi: (value - eta) ** 2 for chi, value in x.items()}
    product_amplitude = (x[-1] - eta) * (x[1] - eta)
    wait = F(1)
    sin_error = {chi: eta / value for chi, value in x.items()}
    cosine = (
        root_lo(1 - sin_error[-1] ** 2) * root_lo(1 - sin_error[1] ** 2)
        - sin_error[-1] * sin_error[1]
        - wait * (b1hi - b1lo) / 2
    )
    assert cosine > 0

    margin0, margin1 = F(1, 64), F(1, 128)
    assert phases(max(-n0lo, n0hi))[1] < F(1, 8) - margin0
    assert phases(2 - nnu)[0] > F(1, 8) + margin0
    assert phases(b1lo)[0] > F(1, 4) + margin1
    assert phases(b1hi)[1] < F(9, 32) - margin1
    assert phases(root_lo(F(5)) / 2 + gap)[0] > F(9, 32) + margin1
    regs = 2**24
    compiler = F(1, 2**20)
    beta0 = phase_error(regs, margin0) + compiler
    beta1 = phase_error(regs, margin1) + compiler
    error = 2 * beta0 + beta1
    g_overlap = root_lo(1 - e0**2)
    sig_min = g_overlap * product_amplitude / 2 - error
    bell = g_overlap * product_amplitude * cosine / 2 - error
    assert bell > 0 and sig_min > 0
    return {
        "archive_revision": SOURCE_REVISION,
        "archive_sha256": {
            name: sha256(data).hexdigest() for name, data in source.items()
        },
        "source": {
            "name": "orthogonal Wilson–Maxwell finite open cube",
            "completion": F(21, 2),
            "coupling": F(1, 64),
            "total_photon_cap": 3,
            "root_and_receiver_sites": [[0, 0, 0], [1, 1, 1]],
            "isolated_rank_in_each_signed_charge_sector": 2,
            "gap_above_full_charged_cluster_over_m_lower": gap,
            "charged_minus_interacting_neutral_minus_reference_cost_four_over_m": energy[
                "charged_minus_actual_neutral_minus4_over_m"
            ],
            "reference_fixed_root": True,
        },
        "exact_free_orbital_checks": exact,
        "neutral_ground_free_subspace_angle_upper": e0,
        "charged_code_free_subspace_angle_upper": e1,
        "combined_wavefunction_difference_upper": eta,
        "actual_ground_insertion_band_residue_lower": {
            "root_spin_minus_and_opposite_spin_plus": residue[-1],
            "root_spin_plus_and_opposite_spin_minus": residue[1],
        },
        "ideal_opposite_corner_on_band_minimum_singular_value_lower_for_any_wait": product_amplitude,
        "actual_held_wait_in_m_inverse": wait,
        "neutral_clock": {
            "dimension": regs,
            "acceptance_radius_in_turns": F(1, 8),
            "alias_safe_margin_in_turns": margin0,
            "coherent_projection_error_upper": beta0,
        },
        "charged_clock": {
            "dimension": regs,
            "acceptance_in_turns": [F(1, 4), F(9, 32)],
            "alias_safe_margin_in_turns": margin1,
            "coherent_projection_error_upper": beta1,
        },
        "each_clock_implementation_operator_error_allowance": compiler,
        "includes_explicit_half_amplitude_factor_from_complete_pump_and_extraction": True,
        "all_complete_record_failures_included": True,
        "original_preparation_code_admission_lower": pcode,
        "complete_finite_pump_clock_extraction_minimum_singular_value_lower": sig_min,
        "complete_signed_bell_rearrival_probability_lower_per_original_start": pcode
        * bell**2
        / 2,
        "complete_receiver_spin_unread_probability_lower_per_original_start_either_sign": pcode
        * sig_min**2,
        "complete_orthogonal_source_matrix_or_clock_tape_executed": False,
        "autonomous_packet_speed_or_endpoint_mass_assigned": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = HERE / "results.json"
    actual = json.dumps(run(), default=str, sort_keys=True, indent=2) + "\n"
    if args.write:
        receipt.write_text(actual)
    else:
        assert receipt.read_text() == actual
        print(
            "Verified cap-three isolated charged band and finite corner spin extraction bounds."
        )


if __name__ == "__main__":
    main()
