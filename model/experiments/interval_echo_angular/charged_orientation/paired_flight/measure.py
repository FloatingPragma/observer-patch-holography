"""Sharp charge arrivals, spatial turns and energy reads on a paired material guide."""

import argparse
from functools import cache
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.experiments.interval_echo_angular.charged_orientation.live_trimer import (
    detectors,
)
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet.measure import (
    half_pairs,
    compile_turn,
)
from model.experiments.interval_echo_gauge.flag_measure.measure import meter
from model.candidate_laws.interval_echo.openwave.row_11_spectrum.conjugate_flight.measure import (
    binomial_interval,
)
from .material import (
    RAILS,
    exact_lines,
    native_matrix,
    pair,
    preparation,
    spectral_word,
)

HERE = Path(__file__).resolve().parent
SHOTS = 2**22
TIMES = list(range(1, 25))
N = len(RAILS[0])
START = [1, N + 1]
END = [N - 2, 2 * N - 2]


def logical_sources(k, sign):
    """Matrices for local Pauli eigenstates before the two buffer-to-gateway swaps."""
    x, y, z = detectors.X, detectors.Y, detectors.Z
    blank = pair(k, sign)
    half = expm(-0.5j * np.pi * blank)
    rate = -k * ((1 + sign) ** 2 - 1)
    angle = np.pi / 2
    time = (angle * np.sign(rate) % (2 * np.pi)) / abs(rate)
    phase = np.diag([np.exp(1j * time * rate), 1])
    native = []
    for basis, pulse in ((z, np.eye(2)), (x, phase @ half), (y, half)):
        columns = [pulse[:, i] for i in (0, 1)]
        columns.sort(key=lambda c: -np.vdot(c, basis @ c).real)
        for direction, c in zip((1, -1), columns, strict=True):
            np.testing.assert_allclose(basis @ c, direction * c, atol=2e-13)
        native.extend(columns)
    return native


def sample(probabilities, rng):
    raw = np.asarray(probabilities).real
    if min(raw) < -5e-12:
        raise ValueError("negative outcome probability")
    p = np.maximum(raw, 0)
    np.testing.assert_allclose(sum(p), 1, atol=5e-11)
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def read(kraus, sources, rng):
    raw = []
    for src in sources:
        arrived = kraus @ src
        fail = 1 - np.vdot(arrived, arrived).real
        raw.append(
            [
                sample(
                    [
                        abs(np.vdot(sources[j], arrived)) ** 2,
                        abs(np.vdot(sources[j + 1], arrived)) ** 2,
                        fail,
                    ],
                    rng,
                )
                for j in (0, 2, 4)
            ]
        )
    pt = detectors.measured_partial(raw, SHOTS)
    p = np.array(raw)[:, 0, :2].sum() / (6 * SHOTS)
    return dict(
        sharp_and_coherent=raw,
        measured_mixed_arrival=float(p),
        unnormalized_partial_transpose_min=float(pt[0]),
        analytic_singular_squared=(
            np.linalg.svd(kraus, compute_uv=False) ** 2
        ).tolist(),
    )


@cache
def turn_pair():
    r, s = half_pairs()[0]
    native = []
    words = []
    for perm in (r, s):
        q = np.eye(12)
        word = compile_turn(perm)
        for a, b in word:
            q[[a, b]] = q[[b, a]]
        np.testing.assert_array_equal(q, np.eye(12)[:, perm])
        native.append(q)
        words.append([list(map(int, x)) for x in word])
    return (*native, words)


def spatial_read(arrival, sources, rng):
    r, s, _ = turn_pair()
    frame = np.eye(12)[:, [a[-1] % 12 for a in RAILS]]
    x, z = detectors.X, detectors.Z
    counts = {}
    for label, left, right in [
        ("perpendicular_half_turn_orders", r @ s @ frame, s @ r @ frame),
        ("one_full_spatial_turn", r @ r @ frame, frame),
        ("buffer_swap_and_local_phase_orders", frame @ x @ z, frame @ z @ x),
    ]:
        left, right = left @ arrival, right @ arrival
        rows = []
        for v in sources:
            bright = np.linalg.norm((left + right) @ v) ** 2 / 4
            dark = np.linalg.norm((left - right) @ v) ** 2 / 4
            rows.append(sample([bright, dark, 1 - bright - dark], rng))
        counts[label] = rows
    return counts


def spectral_read(k, sign, h, rng):
    path = h[:N, :N]
    values, modes = np.linalg.eigh(path)
    diagonal = path[1, 1]
    i = int(np.argmin(abs(values - diagonal)))
    mode = modes[:, i]
    prepared, contacts = spectral_word(mode, k, sign)
    err = float(np.linalg.norm(prepared - mode))
    if err > 1e-11:
        raise AssertionError("Native exposure synthesis missed the selected line")
    # Require the complete one-defect arrival code, over both rails. A
    # physical terminal word has one nonunit load on exactly one rail.
    tags = meter()
    assert [tags[x] for x in (0, 1, 2)] == [0, 1, 2]
    charge_energy = []
    for launch in range(2):
        psi = np.zeros(2 * N, complex)
        psi[launch * N : (launch + 1) * N] = prepared
        at = expm(-13j * h) @ psi
        site = sample(abs(at) ** 2, rng)
        reshaped = np.array(site).reshape(2, N)
        at_ends = int(reshaped[:, [0, -1]].sum())
        static = k * (2 - at_ends / SHOTS + ((1 + sign) ** 2 - 1) * at_ends / SHOTS)
        # A uniformly sampled one of 24 physical bonds has odd-occupation
        # exchange parity probability given by its antisymmetric projection.
        odd = sum(
            abs(at[r * N + j] - at[r * N + j + 1]) ** 2
            for r in range(2)
            for j in range(N - 1)
        ) / (4 * (N - 1))
        minus = int(rng.binomial(SHOTS, odd))
        energy = static + 2 * (N - 1) * minus / SHOTS
        interval_bond = 2 * (N - 1) * binomial_interval(minus, SHOTS, 0.01 / 12)
        interval_load = (
            k * ((1 + sign) ** 2 - 2) * binomial_interval(at_ends, SHOTS, 0.01 / 12)
        )
        interval = 2 * k + np.sort(interval_load) + interval_bond
        charge_energy.append(
            dict(
                port_counts=site,
                antisymmetric_bonds=minus,
                simultaneous_99_interval=interval.tolist(),
                measured_energy=float(energy),
                receiver_probability=float(abs(at[END[launch]]) ** 2),
            )
        )
    # Change the relative phase on one rail with a physical pi phase at its
    # blank private buffer. Positive swaps address each chosen port and return.
    # Run the inverse of the entire eigenmode source to read the two original
    # orbital amplitudes and their orthogonal complement in fresh copies.
    reversal = {}
    starts = np.kron(np.ones(2) / np.sqrt(2), prepared)
    for label, sites in [
        ("one_receiver_gateway", [N - 2]),
        ("whole_second_rail", range(N)),
    ]:
        final = starts.copy()
        for site in sites:
            phase_clock = -k * ((1 + sign) ** 2 - 1)
            dwell = np.pi / abs(phase_clock)
            final[N + site] *= np.exp(1j * phase_clock * dwell)
        projected = np.array(
            [np.vdot(prepared, final[j * N : (j + 1) * N]) for j in (0, 1)]
        )
        plus = abs(sum(projected)) ** 2 / 2
        minus = abs(projected[0] - projected[1]) ** 2 / 2
        reversal[label] = dict(
            plus_minus_other_line=sample([plus, minus, 1 - plus - minus], rng),
            theory_plus_minus=[float(plus), float(minus)],
        )
    # A finite sequence of controlled free waits against an energy reference
    # implements the power-of-two phase filter. It is a supplied coherent
    # timing apparatus; the source/receiver and their rail loads are native.
    filters = np.prod(
        (1 + np.exp(-1j * np.outer(2.0 ** np.arange(13) / 8, values - values[i]))) / 2,
        axis=0,
    )
    filter_map = (modes * filters) @ modes.T
    filter_only = abs(filter_map[-2, 1]) ** 2
    before = int(
        np.random.default_rng(9200 + 2 * k + sign).binomial(SHOTS, filter_only)
    )
    post = expm(-13j * path) @ filter_map
    probability = float(abs(post[-2, 1]) ** 2)
    native_selected_weight = float(abs(mode[1]) ** 2)
    out = read(post[-2, 1] * np.eye(2), logical_sources(k, sign), rng)
    out["simultaneous_filter_refusal_charge_miss_hit"] = sample(
        [
            1 - float(np.sum(abs(filters * modes[1]) ** 2)),
            float(np.sum(abs(filters * modes[1]) ** 2)) - probability,
            probability,
        ],
        rng,
    )
    return dict(
        filtered_receiver_before_extra_flight=before,
        unfiltered_receiver_before_flight=0,
        level=float(values[i]),
        gap=float(min(values[i] - values[i - 1], values[i + 1] - values[i])),
        target_rule="closest to the diagonal energy at the launch gateway",
        eigenspace_multiplicity=2,
        full_path_spectrum=values.tolist(),
        two_rail_gateway_weight=native_selected_weight,
        pulse_state_error=err,
        addressed_contacts_per_rail=contacts,
        independent_energy_banks=charge_energy,
        in_line_logical_phase=reversal,
        filtered_spatial_orders=spatial_read(
            post[-2, 1] * np.eye(2), logical_sources(k, sign), rng
        ),
        filter_conditional_line_weight=float(
            abs(mode[1]) ** 2 / np.sum(abs(filters * modes[1]) ** 2)
        ),
        filtered_gateway_read=out,
        twice_energy_endpoint_determinant=exact_lines(k, sign),
    )


def run():
    geometry = preparation()
    rng = np.random.default_rng(42897)
    cases = []
    for k, sign in ((8, -1), (8, 1), (4, -1)):
        h = native_matrix(k, sign)
        np.testing.assert_allclose(h, np.kron(np.eye(2), h[:N, :N]), atol=2e-14)
        vals, v = np.linalg.eigh(h[:N, :N])
        calibration = [
            int(
                rng.binomial(
                    2**16, abs(sum(v[1] * v[-2] * np.exp(-1j * t * vals))) ** 2
                )
            )
            for t in TIMES
        ]
        if k == 8 and sign == -1:
            flight = TIMES[int(np.argmax(calibration))]
        arrival = expm(-1j * flight * h)[np.ix_(END, START)]
        sources = logical_sources(k, sign)
        cut = expm(-1j * flight * native_matrix(k, sign, cut=True))[np.ix_(END, START)]
        cases.append(
            dict(
                kappa=k,
                charge=sign,
                flight=flight,
                training_gateway_hits_out_of_65536=calibration,
                gateway=read(arrival, sources, rng),
                one_cut=read(cut, sources, rng),
                orders=spatial_read(arrival, sources, rng),
                spectral=spectral_read(k, sign, h, rng),
            )
        )
    return dict(
        shots=SHOTS,
        training_times=TIMES,
        rails=RAILS,
        marks=geometry,
        turns=turn_pair()[2],
        pt_simultaneous_99_radius=detectors.simultaneous_bound(SHOTS, 9),
        cases=cases,
    )


def verify(full, stored):
    """Compare every archived apparatus setting, raw count and numeric reduction."""
    if isinstance(stored, dict):
        assert set(full) == set(stored)
        for name, value in stored.items():
            verify(full[name], value)
    elif isinstance(stored, list):
        for observed, expected in zip(full, stored, strict=True):
            verify(observed, expected)
    elif isinstance(stored, float):
        np.testing.assert_allclose(full, stored, rtol=2e-10, atol=2e-11)
    else:
        assert full == stored


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        old = json.loads((HERE / "results.json").read_text())
        verify(json.loads(json.dumps(data)), old)
        print("Paired charge flight, spectral lines and phase-bank checks passed.")


if __name__ == "__main__":
    main()
