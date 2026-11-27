"""Decode complete receiver and clock counts without source amplitudes."""

from math import sqrt, pi
import numpy as np

from .sensitivity import cut_enclosure


def frequencies(n, j, h):
    """Independent open-bond spectrum, including both signs for real records."""
    b = np.diag(np.full(n + 1, 2 * j)) + np.diag(np.full(n, -2 * h), -1)
    energies = np.linalg.svd(b, compute_uv=False)
    signed = np.r_[energies, -energies]
    return np.unique(np.round((signed[:, None] - signed).ravel(), 13))


def leakage(delta, length):
    """Absolute response of the unit-sum periodic Hann window."""

    def geometric(x):
        x = (x + pi) % (2 * pi) - pi
        return (
            np.exp(1j * (length - 1) * x / 2)
            * np.sinc(length * x / (2 * pi))
            / np.sinc(x / (2 * pi))
        )

    offset = 2 * pi / length
    return abs(
        geometric(delta) - (geometric(delta + offset) + geometric(delta - offset)) / 2
    )


def attribute(signal, weights, dt, radius, n, head):
    sectors = {"intact": n, "cut_right": n - head}
    catalog = {k: frequencies(length, 8.0, 0.5) for k, length in sectors.items()}
    length = len(weights)
    t = np.arange(length) * dt
    rows = {}
    for label, other in (("intact", "cut_right"), ("cut_right", "intact")):
        tests = catalog[label]
        tests = tests[(tests >= 0.1) & (tests <= 2)]
        coeff = abs(np.exp(-1j * tests[:, None] * t) @ (weights[:, None] * signal))
        # For a normalized mixture of either sector alone the all-start difference
        # has coefficient l1 norm <= sqrt(d) (Frobenius covariance-norm bound).
        spill = np.sqrt(2 * (sectors[other] + 1)) * leakage(
            (tests[:, None] - catalog[other]) * dt, length
        ).max(axis=1)
        strength = coeff.max(axis=1)
        excess = strength - spill - radius
        best = int(np.argmax(excess))
        rows[label] = dict(
            tested=len(tests),
            witnesses=int((excess > 0).sum()),
            best_frequency=float(tests[best]),
            best_receiver_index=int(coeff[best].argmax()),
            observed=float(strength[best]),
            other_sector_bound=float(spill[best]),
            margin=float(excess[best]),
        )
    return rows


def wavelength(profile, radius):
    """All cos(k) consistent with the five complex spatial recurrences."""
    lo, hi = -1.0, 1.0
    for a, b in zip(profile[2:] + profile[:-2], 2 * profile[1:-1], strict=True):
        if b == 0:
            if abs(a) > 4 * radius:
                return []
            continue
        center = (a * b.conjugate()).real / abs(b) ** 2
        slack = (4 * radius) ** 2 - abs(a - b * center) ** 2
        if slack < 0:
            return []
        dx = sqrt(slack) / abs(b)
        lo, hi = max(lo, center - dx), min(hi, center + dx)
        if lo > hi:
            return []
    return [float(np.arccos(hi)), float(np.arccos(lo))]


def decode(counts, clocks, sites, route, sources, shots, hand_shots):
    shots, hand_shots = int(shots), int(hand_shots)
    if (
        counts.shape != (2, 2, 4096, 7, 2, 8)
        or shots != 2**19
        or hand_shots != 2**24
        or clocks.shape != (4098, 2)
        or not np.array_equal(sites, np.arange(7, 14))
        or route.shape != (14, 2)
        or not np.array_equal(sources, [4, 1])
        or counts.dtype.kind not in "iu"
        or clocks.dtype.kind not in "iu"
        or counts.min() < 0
        or clocks.min() < 0
        or clocks.max() > hand_shots
        or not np.all(counts.sum(axis=-1) == shots)
    ):
        raise ValueError(
            "incomplete receiver enrollment, source pair, clock or outcome totals"
        )
    length = counts.shape[2]
    quadratures = 2 * clocks.astype(float) / hand_shots - 1
    clock = np.unwrap(np.arctan2(quadratures[:, 1], quadratures[:, 0]))
    hand_err = sqrt(2 * np.log(2 * clocks.size / 0.002) / hand_shots) + 8 / 2**26
    angle_error = float(np.arcsin(sqrt(2) * hand_err))
    dt = float((clock[-1] - clock[2]) / (length - 1))
    dt_error = 2 * angle_error / (length - 1)
    nominal = np.r_[0, 1, 2 + np.arange(length) / 8]
    max_clock = float(np.max(abs(clock - nominal)))
    if (
        max_clock > angle_error
        or np.max(np.diff(clock)) >= pi - 2 * angle_error
        or np.max(np.abs(quadratures)) > 1
    ):
        raise ValueError("clock outside registered program or unwrap bound")
    profiles = counts[..., ::2].sum(axis=-1) / shots
    contrast = profiles[..., 0] - profiles[..., 1]
    w = (1 - np.cos(2 * pi * np.arange(length) / length)) / length
    # Each real and imaginary component of a normalized contrast is a sum of
    # two independent Bernoulli means. All tested frequencies are enumerated.
    # The eight weight-roundings give <= 16 / GRID in a source difference.
    hypotheses = 2 * 2 * 2 * 7 * (length + 4 * 30**2)
    coord = np.sqrt((w @ w) * np.log(2 * hypotheses / 0.004) / shots)
    radius = float(sqrt(2) * coord + 16 / 2**26 + 2e-12)
    raw_error = float(np.sqrt(np.log(2 * contrast.size / 0.004) / shots) + 16 / 2**26)
    transformed = np.fft.rfft(w[None, None, :, None] * contrast, axis=2)
    f = 2 * pi * np.fft.rfftfreq(length, dt)
    eligible = (f >= 0.1) & (f <= 2.0)
    discovery = abs(transformed[0, 0]).max(axis=1)
    selected = []
    for _ in range(2):
        i = int(np.argmax(np.where(eligible, discovery, -1)))
        selected.append(i)
        eligible[max(0, i - 3) : i + 4] = False
    held = transformed[0, 1, selected, :]
    singular = np.linalg.svd(held, compute_uv=False)
    rows = []
    for i, h in zip(selected, held):
        rows.append(
            dict(
                bin=i,
                angular_frequency=float(f[i]),
                clock_half_width=float(f[i] * dt_error / (dt - dt_error)),
                bin_spacing=float(f[1]),
                strongest_bond=int(sites[np.argmax(abs(h))]),
                held_peak=float(np.max(abs(h))),
                allowed_wave_numbers=wavelength(h, radius),
                profile_real=h.real.tolist(),
                profile_imag=h.imag.tolist(),
            )
        )
    # Conditional on the supplied fixed cadence. Measured dt and its interval
    # are reported for the independent observer conversion above.
    sector = attribute(contrast[0, 1], w, 1 / 8, radius, 14, 5)
    box = cut_enclosure(sector["intact"]["best_frequency"], dt, dt_error)
    # Absence of every difference leaves only DC and pair sums: 2*(2J-2h)
    # through 2*(2J+2h). Bound using exact frequencies, including DC.
    high = np.r_[frequencies(14, 8.0, 0.5), frequencies(9, 8.0, 0.5)]
    high = high[(abs(high) > 15) | (high == 0)]
    high_spill = [
        float(sqrt(30) * leakage((2 * pi * i / length - high / 8), length).max())
        for i in selected
    ]
    addresses = [[int(route[i - 1, 1]), int(route[i, 0])] for i in sites]
    if addresses != [
        [42, 43],
        [44, 39],
        [104, 105],
        [101, 97],
        [103, 96],
        [106, 107],
        [100, 99],
    ]:
        raise ValueError("receiver words belong to a different physical route")
    decisions = {
        "both_frequency_bins_resolved": bool(
            min(np.max(abs(h)) for h in held) > radius
        ),
        "profiles_require_more_than_one_spatial_vector": bool(
            singular[1] > sqrt(14) * radius
        ),
        "single_real_wavenumber_necessary_test": [
            "compatible" if row["allowed_wave_numbers"] else "excluded" for row in rows
        ],
        "changed_law_zero_response": "compatible"
        if np.max(abs(contrast[1])) <= raw_error
        else "excluded",
        "nominal_cut_only_explanation": "excluded"
        if sector["intact"]["witnesses"]
        else "unresolved",
        "nominal_intact_only_explanation": "excluded"
        if sector["cut_right"]["witnesses"]
        else "unresolved",
        "cut_only_on_meter_box": "excluded"
        if sector["intact"]["observed"] > radius + box["cut_only_ceiling"]
        else "unresolved",
    }
    return dict(
        summary=dict(
            peaks=[r["angular_frequency"] for r in rows],
            held_peaks=[r["held_peak"] for r in rows],
            selected_bins=selected,
            complex_radius=radius,
            second_spatial_singular=float(singular[1]),
            spatial_bound=float(sqrt(14) * radius),
            zero_contrast_max=float(np.max(abs(contrast[1]))),
            zero_contrast_bound=raw_error,
            sector_witnesses={k: v["witnesses"] for k, v in sector.items()},
        ),
        decisions=decisions,
        clocks=dict(
            spacing=dt,
            spacing_error=dt_error,
            angle_error=angle_error,
            max_discrepancy=max_clock,
        ),
        bands=rows,
        sector=sector,
        cut_contact_sensitivity=box,
        dc_and_sum_only_bound=high_spill,
        frequency_unit="per exchange radian",
        sources=sources.tolist(),
        sites=sites.tolist(),
        receiver_ports=[[int(route[i - 1, 1]), int(route[i, 0])] for i in sites],
        outcomes="00,01,10,11 (earlier coin then later coin), wall then no wall for each word",
        shots_per_source_site_menu=shots,
        total_wall_starts=int(counts.sum()),
        sampling="Seeded conditional binomial quantiles of integer-grid joint outcome probabilities.",
    )


def verify(now, expected):
    if isinstance(expected, dict):
        if now.keys() != expected.keys():
            raise ValueError("receipt fields changed")
        for key in expected:
            verify(now[key], expected[key])
    elif isinstance(expected, (list, tuple)):
        if len(now) != len(expected):
            raise ValueError("receipt enrollment differs")
        for item, target in zip(now, expected, strict=True):
            verify(item, target)
    elif isinstance(expected, float):
        np.testing.assert_allclose(now, expected, rtol=0, atol=3e-11)
    elif now != expected:
        raise ValueError("discrete receipt mismatch")
