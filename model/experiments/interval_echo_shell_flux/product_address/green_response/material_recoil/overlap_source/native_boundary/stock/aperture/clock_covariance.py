"""Two independent ordinary Ramsey readers of a single retained source pair."""

from fractions import Fraction

import numpy as np
from scipy.linalg import expm

from . import occupied_release as flight
from . import regional_tomography as meter

N = meter.N
EPS = meter.HOEFFDING
DWELL = 2 * np.pi / 65
LINE = np.array([0, 32, 33])


def bank(rng, p):
    assert np.min(p) >= -3e-13
    p = np.maximum(p, 0).ravel()
    np.testing.assert_allclose(p.sum(), 1, atol=3e-12, rtol=0)
    return rng.multinomial(N, p / p.sum())


def signed_read(rng, z):
    return rng.binomial(N, (1 + np.array([z.real, z.imag])) / 2)


def source_receiver(load):
    """Balanced rails, with a future quadrature sign reporting +sin(a t)."""
    columns = {}
    blocks = flight.lab.contact.blocks()
    full = flight.lab.contact.stock.kernels()[2]
    for magnitude, (baseline, shift) in blocks["full"].items():
        heads = 9 * np.array(
            [
                i + 1
                for i, (_, m) in enumerate(flight.lab.contact.stock.LAW.LABELS)
                if abs(m) == magnitude
            ]
        )
        absent = heads + 3 * load + 1
        filled = heads + 4
        base = full[np.ix_(absent, absent)]
        extra = full[np.ix_(filled, filled)] - base
        np.testing.assert_array_equal(base, baseline)
        np.testing.assert_array_equal(extra, shift)
        first = expm(-1j * DWELL * base / 4)
        second = expm(-1j * DWELL * (base + extra) / 4)
        values, vecs = np.linalg.eigh(-extra / 4)
        effects = []
        for phase in (1, -1, -1j, 1j):
            k = (first + phase * second) / 2
            np.testing.assert_allclose(k @ extra, extra @ k, atol=2e-11)
            effects.append(k.conj().T @ k / 2)
        np.testing.assert_allclose(sum(effects), np.eye(len(base)), atol=1e-13)
        # Verify nondemolition through both full native sectors and the waiting interval.
        for h in (base, base + extra):
            wait = expm(-1j * 64 * h / 4)
            np.testing.assert_allclose(wait @ extra, extra @ wait, atol=3e-10)
        for value, vec in zip(values, vecs.T, strict=True):
            a = int(round(value))
            assert abs(a - value) < 1e-12
            observed = np.array([np.vdot(vec, m @ vec).real for m in effects])
            expected = (
                np.array(
                    [
                        1 + np.cos(a * DWELL),
                        1 - np.cos(a * DWELL),
                        1 + np.sin(a * DWELL),
                        1 - np.sin(a * DWELL),
                    ]
                )
                / 4
            )
            np.testing.assert_allclose(observed, expected, atol=5e-14)
            columns[a] = observed
    return np.array([columns[a] for a in LINE]).T


def clock_calibration(rng, matrix):
    counts = np.column_stack([bank(rng, p) for p in matrix.T])
    measured = counts / N
    assert np.linalg.matrix_rank(measured) == 3
    return counts, np.linalg.pinv(measured)


def native_sequence():
    contact = flight.lab.contact
    _, weights, _, _ = contact.source_counts(full=True)
    P = weights / weights.sum()
    priors = np.array([P, np.outer(P.sum(1), P.sum(0))])
    sectors = np.flatnonzero((contact.stock.MAG != 1) & (np.arange(11) != 0)) * 9
    full = contact.stock.kernels()[2]
    reference = full[np.ix_(sectors + 4, sectors + 4)] / 4
    operations = {}
    for name, load in (("vacancy", 0), ("excess", 2)):
        loc = sectors + 3 * load + 1
        probe = expm(-1j * DWELL * full[np.ix_(loc, loc)] / 4)
        ref = expm(-1j * DWELL * reference)
        operations[name] = [
            (probe + phase * ref) / (2 * np.sqrt(2)) for phase in (1, -1, -1j, 1j)
        ]
    return priors, reference, operations


def two_clocks(native, first, second, age, rng, same_axis=None):
    priors, reference, operations = native
    wait = expm(-1j * age * reference)
    a, b = operations[first], operations[second]
    if same_axis is None:
        # Both heads start in the actual joint law of the six sharp written
        # labels. Operations on the unobserved first head after its read
        # preserve its trace. The second head evolves throughout the wait.
        p = np.array([(abs(m) ** 2).sum(0) for m in a])
        q = np.array([(abs(m @ wait) ** 2).sum(0) for m in b])
        laws = [p @ prior @ q.T for prior in priors]
    else:
        # Same physical head, including the disturbance of the first clock.
        pairs = np.array(
            [[(abs(last @ wait @ initial) ** 2).sum(0) for last in b] for initial in a]
        )
        laws = [pairs @ prior.sum(axis=same_axis) for prior in priors]
    return np.array([bank(rng, p).reshape(4, 4) for p in laws]), np.array(laws)


def frequency_laws(counts, left, right):
    out = []
    for samples in counts:
        rows = []
        for i in range(3):
            row = []
            for j in range(3):
                dual = np.outer(left[i], right[j])
                middle = np.sum(samples / N * dual)
                error = 5 * EPS * np.ptp(dual) + 1e-12
                row.append([middle - error, middle + error])
            rows.append(row)
        out.append(rows)
    return out


def estimate(counts, left, right, reading):
    # A constant score has exactly zero arm difference. Remove it before
    # computing dual weights to keep statistical errors small.
    offset = (reading.min() + reading.max()) / 2
    dual = left.T @ (reading - offset) @ right
    middle = np.sum((counts[1] - counts[0]) / N * dual)
    span = np.ptp(dual)
    # TV per calibrated detector column <= 2 eps. Two readers, two arms;
    # add two independent bounded terminal means (one for each arm).
    width = 10 * EPS * span + 1e-13 * (1 + np.max(abs(dual)))
    return middle + np.array([-width, width])


def gain(load, rng):
    # Two active alternatives launch on the same native blank contact.
    # Exchange out of the monitored word contributes zero fringe. A balanced
    # coin assigns misses equal signed weights on the two terminal outputs.
    h = flight.lab.contact.stock.kernels()[2][:9, :9] / 4
    dwell = 1 / 64
    u = expm(-1j * dwell * h)
    records = []
    for c in (1, 2):
        src, ref = 3 * load + c, 3 + c
        z = u[src, src] * u[ref, ref].conjugate()
        rate = h[src, src] - h[ref, ref]
        assert abs(rate * dwell) < np.pi / 2
        np.testing.assert_allclose(-np.angle(z) / dwell, rate, atol=1e-11)
        count = signed_read(rng, z)
        measured = (2 * count / N - 1) @ np.array([1, 1j])
        phase_width = np.arcsin(2 * np.sqrt(2) * EPS / abs(measured))
        band = (-np.angle(measured) + np.array([-phase_width, phase_width])) / dwell
        assert band[0] <= rate <= band[1]
        records.append(dict(positive=count.tolist(), frequency=band.tolist()))
    response = meter.difference_interval(
        records[1]["frequency"], records[0]["frequency"]
    )
    # The native blank word has a half-integer diagonal. Its difference
    # and hence this difference of differences also lie on that lattice.
    np.testing.assert_array_equal(2 * h.diagonal(), np.rint(2 * h.diagonal()))
    low, high = int(np.ceil(2 * response[0])), int(np.floor(2 * response[1]))
    assert low == high
    return dict(
        dwell=dwell,
        returns=records,
        gain=response.tolist(),
        certified_gain=str(Fraction(low, 2)),
    )


def experiment(board, old, edges, rng):
    adj = board.adjacency().toarray().astype(bool)
    for clamp in range(board.n):
        ends = [[s for s in edge if adj[clamp, s]] for edge in edges]
        if clamp not in np.ravel(edges) and all(ends):
            stations = [min(x) for x in ends]
            break
    else:
        raise ValueError("no common ordinary calibration site")
    head = np.array(flight.sector_prior(old))
    law = (head @ flight.CONDITIONAL).reshape(2, 3, 3)
    delta = flight.lab.probe.log_factor(old["source"])
    matrices = {name: source_receiver(t) for name, t in (("vacancy", 0), ("excess", 2))}
    gains = {}
    counts, inverses = {}, {}
    for species in ("vacancy", "excess"):
        gains[species] = gain(0 if species == "vacancy" else 2, rng)
        counts[species], inverses[species] = clock_calibration(rng, matrices[species])
    sequence = native_sequence()
    both = {}
    clock = {}
    measured = {}
    freq = LINE - 16
    for age in (0, 64):
        for first, second in (
            ("vacancy", "vacancy"),
            ("vacancy", "excess"),
            ("excess", "excess"),
        ):
            pair = (age, first, second)
            clock[pair], records = two_clocks(sequence, first, second, age, rng)
            np.testing.assert_allclose(
                records, matrices[first] @ law @ matrices[second].T, atol=2e-12, rtol=0
            )
            measured[pair] = estimate(
                clock[pair], inverses[first], inverses[second], np.outer(freq, freq)
            )
            both[f"{age}/{first}/{second}"] = dict(
                counts=clock[pair].tolist(), covariance=measured[pair].tolist()
            )
    ratio = meter.product_interval(
        gains["excess"]["gain"], 1 / np.array(gains["vacancy"]["gain"])
    )
    vv = measured[0, "vacancy", "vacancy"]
    ordinary_cross = meter.product_interval(vv, ratio)
    ordinary_ee = meter.product_interval(ordinary_cross, ratio)
    source_covariance = delta * (65 / 4) ** 2
    for band in measured.values():
        assert max(band[0], source_covariance[0]) <= min(band[1], source_covariance[1])
    certified = Fraction(gains["excess"]["certified_gain"]) / Fraction(
        gains["vacancy"]["certified_gain"]
    )
    shared = {}
    for axis, station in ((2, stations[0]), (1, stations[1])):
        # Two different defects interrogate the SAME retained head in order.
        marginal = law.sum(axis=axis)
        source = np.array([np.diag(p) for p in marginal])
        same, records = two_clocks(
            sequence, "vacancy", "excess", 64, rng, same_axis=axis - 1
        )
        np.testing.assert_allclose(
            records,
            matrices["vacancy"] @ source @ matrices["excess"].T,
            atol=2e-12,
            rtol=0,
        )
        dual = inverses["vacancy"].T @ (LINE[:, None] - LINE) ** 2 @ inverses["excess"]
        span = np.ptp(dual)
        widths = 5 * EPS * span + 1e-13 * (1 + np.max(abs(dual)))
        center = np.sum(same / N * dual, axis=(1, 2))
        lo, hi = center - widths, center + widths
        assert np.all(lo <= 0) and np.all(hi >= 0)
        difference = meter.difference_interval(
            gains["vacancy"]["gain"], gains["excess"]["gain"]
        )
        shared[str(station)] = dict(
            counts=same.tolist(),
            squared_line_difference=np.stack((lo, hi), axis=1).tolist(),
            scalar_mean_square_upper=(hi / min(abs(difference)) ** 2).tolist(),
        )
    shared["definition"] = (
        "(vacancy frequency shift - excess frequency shift) / (vacancy ordinary gain - excess ordinary gain)"
    )
    return (
        clamp,
        law,
        delta,
        clock[0, "vacancy", "vacancy"],
        inverses["vacancy"],
        certified,
        dict(
            stations=stations,
            ordinary_clamp=clamp,
            gains=gains,
            local_frequency_calibration={k: v.tolist() for k, v in counts.items()},
            source_covariance=source_covariance.tolist(),
            clock_pairs=both,
            inferred_frequency_laws=frequency_laws(
                clock[0, "vacancy", "vacancy"], inverses["vacancy"], inverses["vacancy"]
            ),
            ordinary_gain_ratio=ratio.tolist(),
            certified_ratio=str(certified),
            same_head=shared,
            matched_field_cross=ordinary_cross.tolist(),
            matched_field_excess=ordinary_ee.tolist(),
        ),
    )
