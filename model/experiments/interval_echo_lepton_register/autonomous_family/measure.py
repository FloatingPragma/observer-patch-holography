"""Apply a static mass term to the full twelve-port analyzer and read its field."""

import argparse
import json
from fractions import Fraction
from math import factorial, isqrt
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame import (
    field,
    measure as face,
)
from model.executable_model import cell, rotations

from ..phase_feedback import circuit
from .decay import weak_records
from .face_locality import reencode, two_current_records
from .reconstruct import (
    RECOIL_SETTINGS,
    geometry_record,
    native_frequency_check,
    native_phase_banks,
    ramsey_records,
    recoil_record,
    reconstruct,
    terminal_record,
)

PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
BETA = np.kron(PAULI[0], np.eye(2))
ALPHA = np.stack([np.kron(-PAULI[2], sigma) for sigma in PAULI])


def port_current(port_rotation, source, p):
    """Intertwine the analyzer rails and the neutral A5 band in the same orientation."""
    bands, tables = face.theorem_projectors()
    band = bands[1][0]
    oriented_face = next(
        f for f in face.faces() if list(face.face_turn(f)) == port_rotation
    )
    frames = face.character_frames(oriented_face, band)
    # The analyzer Fourier columns carry eigenvalues 1, omega^-1, omega^-2.
    axes = []
    for projector in frames[[0, 2, 1]]:
        col = projector[:, np.argmax(projector.diagonal().real)]
        axes.append(col / np.linalg.norm(col))
    neutral = np.column_stack(axes)
    encoding = neutral @ source.conj().T
    turn = face.permutation(port_rotation)
    assert_allclose(turn @ encoding, encoding @ turn, atol=2e-14)
    assert_allclose(encoding.conj().T @ encoding, p, atol=2e-14)
    assert_allclose(encoding @ encoding.conj().T, band, atol=2e-14)

    raising = field.transfer(3, 4)
    occupation = np.stack([np.diag(field.transfer(a, a)) for a in range(5)])
    electric = -occupation[:3].sum(axis=0) / 3 + occupation[3]
    assert electric[[23, 15, 24]].tolist() == [-1, 0, 1]
    assert np.flatnonzero(raising[:, 23]).tolist() == [15]
    assert abs(raising[15, 23]) == 1 and not raising[:, 24].any()
    color = [field.transfer(a, b) for a in range(3) for b in range(3) if a != b]
    color += [
        field.transfer(0, 0) - field.transfer(1, 1),
        field.transfer(1, 1) - field.transfer(2, 2),
    ]
    assert not np.any(np.stack(color)[:, :, [23, 15, 24]])
    assert np.linalg.norm(np.stack(color)[:, :, 1]) > 0
    current = raising[15, 23] * encoding

    # Check all rotated realizations of the same rank-one family projectors.
    for rotation in rotations():
        g = face.permutation(rotation)
        rotated = tuple(rotation[i] for i in oriented_face)
        assert_allclose(
            g @ frames @ g.T, face.character_frames(rotated, band), atol=2e-14
        )
    evidence = dict(
        oriented_face=oriented_face,
        rotations_checked=60,
        charge=electric[[23, 15, 24]].tolist(),
        current_sign=raising[15, 23],
        family_band=tables,
    )
    return current, neutral, evidence


def field_symbol(mass, momentum, spacing):
    """The squared Wilson term lifts the corners without changing the rest curvature."""
    if spacing:
        angle = spacing * momentum
        kinetic = np.sin(angle) / spacing
        wilson = (2 * np.sin(angle / 2) ** 2).sum() ** 2 / spacing
    else:
        kinetic = momentum
        wilson = 0.0
    identity = np.eye(len(mass))
    shifted = mass + wilson * identity
    symbol = np.kron(np.einsum("i,ijk->jk", kinetic, ALPHA), identity) + np.kron(
        BETA, shifted
    )
    return symbol, kinetic, wilson


def same_field_observations(mass, p, modes, energies):
    """Read the full spinor spectrum and the same-mode rest and recoil response."""
    all_masses = np.linalg.eigvalsh(mass)
    banks = []
    spin_z = np.kron(np.eye(2), PAULI[2]) / 2
    spin_state = np.array([1, 0, 1, 0]) / np.sqrt(2)
    states = np.kron(spin_state[:, None], modes)
    assert_allclose(np.kron(BETA, mass) @ states, states * energies, atol=2e-14)
    assert_allclose(np.kron(spin_z, p) @ states, states / 2, atol=2e-14)
    ramsey = ramsey_records(np.kron(BETA, mass), states)
    for energy, row in zip(energies, ramsey["preparations"], strict=True):
        low, high = row["energy_interval_99"]
        assert low < energy < high
    rng = np.random.default_rng(726044)
    recoil = [
        recoil_record(field_symbol(mass, k, 0)[0], states, k, rng)
        for k in np.array(RECOIL_SETTINGS)
    ]

    # Full turn applied to the same positive-energy, charged rest preparations.
    turn = np.cos(np.pi) * np.eye(4) - 2j * np.sin(np.pi) * spin_z
    assert_allclose(np.kron(turn, p) @ states, -states, atol=2e-14)
    # Exterior lift on the same two modes, preserving both orderings.
    for state in states.T:
        first, second = np.kron(np.eye(2), state).reshape(2, 2 * len(state))
        antisymmetric = (np.outer(first, second) - np.outer(second, first)) / np.sqrt(2)
        assert_allclose(np.linalg.norm(antisymmetric), 1, atol=2e-14)
        assert_allclose(antisymmetric.T, -antisymmetric, atol=2e-14)

    for spacing in (0.0, 1.0, 1 / 8, 1 / 64):
        for momentum in np.array(
            [[0, 0, 0], [0.01, 0, 0], [0, 0.1, 0], [0.2, 0.3, 0.4]]
        ):
            h, k, w = field_symbol(mass, momentum, spacing)
            conjugate, _, _ = field_symbol(mass.conj(), momentum, spacing)
            levels = np.linalg.eigvalsh(h)
            assert_allclose(np.linalg.eigvalsh(conjugate), levels, atol=3e-14)
            positive = levels[len(mass) * 2 :]
            expected = np.hypot(np.linalg.norm(k), all_masses + w)
            assert_allclose(positive, np.repeat(expected, 2), atol=3e-14)
            restriction = np.kron(np.eye(4), p)
            assert_allclose(h @ restriction, restriction @ h, atol=2e-14)
            banks.append(
                dict(
                    spacing=spacing,
                    momentum=momentum.tolist(),
                    positive_energies=positive[::2].tolist(),
                )
            )
    return dict(
        full_turn_phase=-1,
        aligned_exchange_phase=-1,
        records=banks,
        rest_interferometer=ramsey,
        recoil_interferometers=recoil,
    )


def ratio_enclosures():
    """Rational enclosures for the exact angle and character ratios."""
    angle = Fraction(2, 9)

    def trig(shift):
        terms = [
            (-1) ** j * angle ** (2 * j + shift) / factorial(2 * j + shift)
            for j in range(9)
        ]
        return sum(terms[:8]), sum(terms)

    def radical(value):
        denominator = 10**25
        root = isqrt(value * denominator**2)
        return Fraction(root, denominator), Fraction(root + 1, denominator)

    c0, c1 = trig(0)
    s0, s1 = trig(1)
    r0, r1 = radical(2)
    t0, t1 = radical(3)
    light = r0 - c1 - t1 * s1, r1 - c0 - t0 * s0
    middle = r0 - c1 + t0 * s0, r1 - c0 + t1 * s1
    heavy = r0 + 2 * c0, r1 + 2 * c1
    assert 0 < light[0] < light[1] < middle[0] < middle[1] < heavy[0]
    return [((lo / light[1]) ** 2, (hi / light[0]) ** 2) for lo, hi in (middle, heavy)]


def comparisons(source, p, q, mass, unitary, gap):
    """Compare two other static attachments and eliminate the common virtual gap."""
    raw = np.linalg.eigvalsh(source.conj().T @ q @ source)
    regular = np.linalg.eigvalsh(mass)[:3]
    assert_allclose(regular, raw, atol=3e-14)
    raw_ratio = raw[1:] / raw[0]

    single = np.linalg.eigvalsh(q + gap * (np.eye(len(p)) - p))[:3]
    expected_single = (
        2 * gap * raw / (gap + 1 + np.sqrt((gap + 1) ** 2 - 4 * gap * raw))
    )
    assert_allclose(single, expected_single, atol=3e-14)
    return_amplitude = source.conj().T @ unitary @ source
    virtual_h = np.block(
        [
            [np.zeros((3, 3)), return_amplitude.conj().T],
            [return_amplitude, gap * np.eye(3)],
        ]
    )
    virtual = -np.linalg.eigvalsh(virtual_h)[:3][::-1]
    expected_virtual = 2 * raw / (gap + np.sqrt(gap**2 + 4 * raw))
    assert_allclose(virtual, expected_virtual, atol=3e-14)

    # Only this comparison uses laboratory values. Use the direct CODATA ratio.
    electron, electron_u = map(Fraction, ("0.51099895069", "0.00000000016"))
    tau, tau_u = map(Fraction, ("1776.93", "0.09"))
    mu_ratio, mu_ratio_u = map(Fraction, ("206.7682827", "0.0000046"))
    mu_bounds = [mu_ratio - 2 * mu_ratio_u, mu_ratio + 2 * mu_ratio_u]
    tau_bounds = [
        (tau - 2 * tau_u) / (electron + 2 * electron_u),
        (tau + 2 * tau_u) / (electron - 2 * electron_u),
    ]
    mu_interval, tau_interval = (
        np.array(row, dtype=float) for row in (mu_bounds, tau_bounds)
    )

    def required_eta(q_ref, interval):
        # The inverse decreases on 1 < r <= q_ref. Negative values are inadmissible.
        low, high = interval[::-1]
        return [
            max(0.0, (q_ref - low) / (low**2 - q_ref)),
            (q_ref - high) / (high**2 - q_ref),
        ]

    mu_eta = required_eta(raw_ratio[0], mu_interval)
    tau_eta = required_eta(raw_ratio[1], tau_interval)
    assert mu_eta[0] > tau_eta[1]
    # Conservative rational bounds used in the independent Lean comparison.
    exact = ratio_enclosures()
    assert exact[0][0] >= Fraction("206.770315") and mu_bounds[1] <= Fraction(
        "206.768292"
    )
    assert exact[1][1] <= Fraction("3477.473") and tau_bounds[0] >= Fraction("3477.012")
    assert mu_bounds[0] >= 206
    assert_allclose(
        raw_ratio, [float((a + b) / 2) for a, b in exact], atol=1e-7, rtol=0
    )

    endpoint_eta = np.array(mu_eta)
    rhs = raw_ratio[1] * (1 + endpoint_eta)
    tau_at_mu = 2 * rhs / (1 + np.sqrt(1 + 4 * endpoint_eta * rhs))
    assert max(tau_at_mu) < tau_interval[0]
    return dict(
        quotient=float(regular.sum() / np.sqrt(regular).sum() ** 2),
        rest_energies=regular.tolist(),
        all_rest_energies=np.linalg.eigvalsh(mass).tolist(),
        paired_ratios=raw_ratio.tolist(),
        single_projector_ratios=(single[1:] / single[0]).tolist(),
        virtual_ratios=(virtual[1:] / virtual[0]).tolist(),
        middle_ratio_error_ppm=float(1e6 * (raw_ratio[0] / float(mu_ratio) - 1)),
        lab_windows=[mu_interval.tolist(), tau_interval.tolist()],
        virtual_eta_windows=[mu_eta, tau_eta],
        virtual_tau_at_allowed_mu=sorted(tau_at_mu.tolist()),
    )


def run():
    sites, rotation, theta, unitary, _, fourier, _, _ = circuit.build()
    address = np.eye(12)[:, sites[:, 0]]
    p = address @ address.T
    raw_q = unitary.conj().T @ p @ unitary
    reflection = 2 * p - np.eye(12)
    gap = circuit.marked.KAPPA
    mass = (raw_q + reflection @ raw_q @ reflection) / 2 + gap * (np.eye(12) - p)
    assert_allclose(mass, mass.conj().T, atol=2e-14)
    assert_allclose(p @ mass, mass @ p, atol=2e-14)
    assert_allclose(p @ mass @ p, p @ raw_q @ p, atol=2e-14)
    assert (
        np.linalg.eigvalsh(mass)[3] > 1 and np.linalg.eigvalsh(mass)[3] >= gap - 2e-14
    )

    source = address @ fourier
    energies = np.einsum("ij,ij->j", source.conj(), mass @ source).real
    assert_allclose(mass @ source, source * energies, atol=2e-14)
    current, neutral, gauge = port_current(rotation, source, p)
    active = neutral.conj().T @ current @ source
    assert_allclose(active, gauge["current_sign"] * np.eye(3), atol=2e-14)

    native = circuit.release()
    assert native["minimum_leakage_squared"] > 0.74
    support = np.triu(np.abs(mass) > 1e-12, 1)
    beyond = support & ~cell().adjacency().toarray().astype(bool)
    local_mass, local_p, local_modes, local_current, locality = reencode(
        address,
        face.permutation(rotation),
        gauge["oriented_face"],
        mass,
        source,
        current,
        gap,
    )
    locality["original_off_seam_pairs"] = int(beyond.sum())
    local_active = neutral.conj().T @ local_current @ local_modes
    assert_allclose(local_active, active, atol=2e-14)
    locality["weak_current_controls"] = two_current_records(
        locality["ports"], local_modes, neutral, local_current
    )
    observations = same_field_observations(local_mass, local_p, local_modes, energies)
    geometric_counts = geometry_record(unitary, source[:, [2, 0, 1]], address[:, 0])
    registered = terminal_record(native["survival_on_three_ports_at_one"])
    native_generator = circuit.free_generator()
    registered_phases = native_phase_banks(native_generator, source)
    diagonal = native_generator.diagonal()
    row_radius = np.asarray(abs(native_generator).sum(axis=1)).ravel() - abs(diagonal)
    native["filled_host_energy_support"] = [
        float(min(diagonal - row_radius)),
        float(max(diagonal + row_radius)),
    ]
    reconstructed = reconstruct(
        geometric_counts,
        observations["rest_interferometer"],
        observations["recoil_interferometers"],
        registered,
    )
    frequency_check = native_frequency_check(registered_phases)
    assert max(bank["lower_bound"] for bank in frequency_check) > 0
    order = reconstructed["formation_tags_light_to_heavy"]
    assert reconstructed["mapping"]["geometry_channel_for_each_formation_tag"] == (
        1,
        2,
        0,
    )
    report = comparisons(source, p, raw_q, mass, unitary, gap)
    assert_allclose(report["rest_energies"], energies[order], atol=2e-14)
    return dict(
        theta=theta,
        spectator_gap=gap,
        central_ports=sites[:, 0].tolist(),
        spectrum=report,
        face_current=gauge,
        field=observations,
        observer_records=dict(
            geometry=geometric_counts,
            native_free_return=registered,
            native_phase_counts=registered_phases,
            native_frequency_check=frequency_check,
            reconstruction=reconstructed,
        ),
        weak=weak_records(energies[order], local_active[np.ix_(order, order)], PAULI),
        locality=locality,
        native_release_control=native,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    else:
        print(json.dumps(result["observer_records"]["reconstruction"], indent=2))
