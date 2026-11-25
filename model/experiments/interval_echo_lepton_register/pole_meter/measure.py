"""Same-charge Ramsey spectroscopy against a remote compact vacancy."""

import argparse
from functools import cache
import json
from math import asin, lcm, log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, eigh
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo.candidate import KAPPA
from model.candidate_laws.interval_echo.openwave.row_04_compact import (
    material_reference,
)
from model.candidate_laws.interval_echo.openwave.row_04_compact.gateway_bound import (
    dark_projector,
)
from model.executable_model import federation
from model.experiments.interval_echo_lepton_register import measure as marked

HERE = Path(__file__).resolve().parent
STEPS = 512
STEP = 0.25
SHOTS = marked.SHOTS
SMALL_TIME = 1 / 64


def compact_preparation(target, sharp):
    operator = np.eye(12, dtype=complex)
    remaining = 1.0
    for port, coefficient in enumerate(target):
        if port != sharp and coefficient:
            marked.rotate(operator, sharp, port, asin(coefficient / remaining), 1)
            remaining = sqrt(remaining**2 - coefficient**2)
    np.testing.assert_allclose(operator[:, sharp], target, atol=2e-13)
    return operator


@cache
def storage_record():
    return json.loads(
        (
            Path(material_reference.__file__).parent / "material_reference.json"
        ).read_text()
    )["native"][0]


def storage_calibration():
    bank = storage_record()["panels"]["6"]["vacancy"]
    counts = np.asarray(bank["detector_counts"])[0]
    signal = (counts[:, 0] - counts[:, 1]) / material_reference.SHOTS
    return material_reference.fit(signal, bank["timer_swaps"])


@cache
def setup():
    board = federation(0)
    charts = [
        tuple(map(int, np.flatnonzero(local == 6)))
        for local in board.degree.reshape(-1, 12)
    ]
    projections = [dark_projector(chart) for chart in charts]
    denominator = lcm(*(int(x.q) for p in projections for x in p))
    projectors = [(np.array(p * denominator, dtype=np.int64)) for p in projections]
    numerator = block_diag(*projectors)
    assert np.array_equal(numerator @ numerator, denominator * numerator)

    stored = storage_record()
    home, remote = stored["home_cells"]
    sharp, reference_tip = stored["launch_ports"]
    assert (
        home == 0
        and charts[home] == charts[remote]
        and reference_tip == 12 * remote + sharp
    )
    column = projectors[home][:, sharp].astype(float)
    target = column / np.linalg.norm(column)
    operator = compact_preparation(target, sharp)
    raw = marked.addresses(marked.CELL_SEAMS[0])
    a, _, plus, minus = raw
    character_ports = [a, minus, plus]
    analyzer = marked.reader(raw, 2 / 9, 1)
    states = np.zeros((board.n, 4), complex)
    states[:12, :3] = analyzer.conj().T[:, character_ports]
    states[:12, 3] = operator[:, sharp]
    ref = np.zeros(board.n, complex)
    ref[12 * remote : 12 * (remote + 1)] = operator[:, sharp]
    parents = shortest_path(
        board.adjacency(), directed=False, unweighted=True, return_predecessors=True
    )[1]
    routes = []
    for source in [*character_ports, sharp]:
        trail = [reference_tip]
        while trail[-1] != source:
            trail.append(int(parents[source, trail[-1]]))
        routes.append(trail[::-1])
    return board, numerator, denominator, states, ref, routes


def geometry_certificate(h, projection, denominator):
    """Integer positivity witness for the spectral gap on the compact complement."""
    shifted = np.rint(2 * (h - 35 * np.eye(len(h)))).astype(np.int64)
    np.testing.assert_array_equal(shifted / 2, h - 35 * np.eye(len(h)))
    assert not np.any(shifted @ projection)
    # Certify A = (H-35)^2 + P - (1/5)^2 I strictly positive, from
    # the exact dyadic factor plus a strictly diagonally dominant remainder.
    scale = lcm(4, 25, denominator)
    eye = np.eye(len(h), dtype=np.int64)
    integer = (
        (scale // 4) * (shifted @ shifted)
        + (scale // denominator) * projection
        - (scale // 25) * eye
    )
    dyadic = 2**16
    factor = np.rint(dyadic * np.linalg.cholesky(integer - 0.5 * eye)).astype(np.int64)
    assert len(h) * int(np.max(abs(factor))) ** 2 < 2**62
    assert int(np.max(abs(integer))) * dyadic**2 < 2**62
    residual = integer * dyadic**2 - factor @ factor.T
    margin = np.diag(residual) - (
        np.abs(residual).sum(axis=1) - np.abs(np.diag(residual))
    )
    if min(margin) <= 0:
        raise ArithmeticError("Dyadic certificate does not establish the gap")
    # Gershgorin gives this coarse enclosure on every one-vacancy eigenvalue.
    radius = np.abs(h).sum(axis=1) - np.diag(h)
    bounds = [float(min(np.diag(h) - radius)), float(max(np.diag(h) + radius))]
    return dict(
        compact_dimension=int(np.trace(projection) // denominator),
        proved_pole=35,
        proved_separation=0.2,
        gershgorin=bounds,
        dyadic_denominator=dyadic,
        certificate_scale=scale,
        minimum_integer_dominance_margin=int(min(margin)),
    )


@cache
def spectrum(speed):
    board, numerator, denominator, states, ref, _ = setup()
    h = speed * board.laplacian().toarray() / 2 + np.diag(
        KAPPA * (2 * board.degree - 6)
    )
    ref_energy = float(np.vdot(ref, h @ ref).real)
    np.testing.assert_allclose(h @ ref, ref_energy * ref, atol=5e-12)
    np.testing.assert_allclose(states.conj().T @ ref, 0, atol=1e-14)
    values, vectors = eigh(h)
    weights = abs(vectors.T @ states) ** 2
    exact_weight = (states.conj() * (numerator @ states)).sum(axis=0).real / denominator
    return h, ref_energy, values, weights, exact_weight


def outcomes(amplitude, phase=0, coherent=True):
    visibility = np.real(amplitude * np.exp(1j * phase)) if coherent else 0
    return np.stack(
        (
            (1 + abs(amplitude) ** 2 + 2 * visibility) / 4,
            (1 + abs(amplitude) ** 2 - 2 * visibility) / 4,
            (1 - abs(amplitude) ** 2) / 2,
        ),
        axis=-1,
    )


def counts(rng, p):
    # Structural dark bins are exact zeros on this fixed detector grid.
    return marked.count_bank(rng, p)


def quadrature(counts_):
    return float((counts_[0] - counts_[1]) / SHOTS)


def measure():
    board, proj, denom, _, _, routes = setup()
    h, ref, energies, weights, direct = spectrum(1)
    cert = geometry_certificate(h, proj, denom)
    delays = np.arange(STEPS) * STEP
    frequencies = energies - ref
    amp = np.exp(-1j * np.outer(delays, frequencies)) @ weights
    law = np.array([STEPS, *(2 * (STEPS - j) for j in range(1, STEPS))]) / STEPS**2
    assert law.sum() == 1
    min_gap = cert["proved_separation"]
    end = max(abs(x - ref) for x in cert["gershgorin"])
    assert end * STEP / 2 < pi and np.sin(end * STEP / 2) > np.sin(min_gap * STEP / 2)
    tail_bound = 1 / (STEPS * np.sin(STEP * min_gap / 2)) ** 2
    mass_read = []
    rng = np.random.default_rng(155121)
    # Each trial independently draws one absolute difference of two uniform
    # integer waits. The observed categories are bright, dark, or neither.
    for label, phase, coherent in [
        ("compact_reference", np.zeros(STEPS), True),
        ("detuned_by_nine", 9 * delays, True),
        ("unread_arm", np.zeros(STEPS), False),
    ]:
        channels = []
        for k in range(4):
            p = outcomes(amp[:, k], phase, coherent)
            marginal = law @ p
            bank = counts(rng, marginal)
            channels.append(
                dict(
                    port_counts=bank.tolist(),
                    mean_from_counts=quadrature(bank),
                    mean_from_law=float(marginal[0] - marginal[1]),
                )
            )
        mass_read.append(dict(menu=label, source_channels=channels))

    # A separate material dimer checks the short exposure used for the centroid.
    # For short waits, two-clock quadrature data obey a uniform cubic bias bound.
    time = SMALL_TIME
    timer_steps = 16
    timer_probability = np.sin(timer_steps * time / 2) ** 2
    exchange_timer = counts(rng, np.array([timer_probability, 1 - timer_probability]))
    estimated_time = 2 * asin(sqrt(exchange_timer[0] / SHOTS)) / timer_steps
    moments = []
    for speed in (1, 3):
        _, anchor, e, w, _ = spectrum(speed)
        z = np.exp(-1j * time * (e - anchor)) @ w
        delta = (e - anchor) @ w
        second = (e - anchor) ** 2 @ w
        abs_third = abs(e - anchor) ** 3 @ w
        banks = [counts(rng, p) for p in outcomes(z, pi / 2)]
        moments.append(
            dict(
                exchange_speed=speed,
                reference_energy=round(anchor, 9),
                generator_centroid=(anchor + delta).tolist(),
                generator_spread=np.sqrt(np.maximum(0, second - delta**2)).tolist(),
                short_wait_port_counts=[b.tolist() for b in banks],
                relative_slope_from_nominal_wait=[quadrature(b) / time for b in banks],
                relative_slope_from_dimer_time=[
                    quadrature(b) / estimated_time for b in banks
                ],
                finite_wait_bias_bounds=(time**2 * abs_third / 6).tolist(),
            )
        )

    # This optional ideal-meter calculation has finer energy resolution than the
    # finite Ramsey bank. Merge numerical degeneracies before any inference.
    groups = np.r_[0, np.flatnonzero(np.diff(energies) > 1e-8) + 1]
    bands = np.add.reduceat(weights[:, :3], groups, axis=0)
    support = [0, 1, 5, 7]
    np.testing.assert_array_equal(
        proj[np.ix_(support, support)] * (6 // denom),
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 2, -1], [0, 0, -1, 2]],
    )
    predicted = (2 - np.cos(4 / 9 + 4 * pi * np.arange(3) / 3)) / 9
    np.testing.assert_allclose(direct[:3], predicted, atol=1e-13)
    eps = sqrt(2 * log(2 * 12 / 0.01) / SHOTS)
    return dict(
        reference_local_chart=[int(p) for p in np.flatnonzero(board.degree[:12] == 6)],
        pulse_routes=routes,
        shots_per_bank=SHOTS,
        step=STEP,
        steps=STEPS,
        source_order=[
            "character_0",
            "character_1",
            "character_2",
            "home_compact_control",
        ],
        nominal_reference_energy=35,
        independent_six_seam_store_calibration=storage_calibration(),
        short_wait=SMALL_TIME,
        short_wait_clock_swaps_and_stays=exchange_timer.tolist(),
        short_wait_clock_estimate=estimated_time,
        clock_repetitions=timer_steps,
        certificate=cert,
        off_pole_filter_bound=tail_bound,
        simultaneous_99_percent_sampling_radius=eps,
        reference_pole_weights=direct.tolist(),
        randomized_wait_banks=mass_read,
        finite_wait_moments=moments,
        ideal_energy_reader=dict(
            numerical_grouping_tolerance=1e-8,
            resolved_bands=len(groups),
            visible_per_source=(bands > 1e-10).sum(axis=0).tolist(),
            equal_prior_guess=float(bands.max(axis=1).sum() / 3),
            rigorous_common_pole_ceiling=float(
                1 - (direct[:3].sum() - max(direct[:3])) / 3
            ),
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = measure()
    dest = HERE / "results.json"
    if parser.parse_args().write:
        dest.write_text(json.dumps(result, indent=2) + "\n")
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
