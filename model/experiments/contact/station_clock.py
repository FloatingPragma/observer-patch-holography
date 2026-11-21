"""Read a receiver's freely evolving written collar when its local scalar detector clicks."""

import math

import numpy as np
from scipy.linalg import expm
from scipy.stats import beta

from model.candidate_laws.scalar_channel.carrier import require
from .law import arm_energies


def local_pointers(spec, rate, elapsed):
    """Read plus-X and minus-Y, so the free written doublet advances with positive phase."""
    native = np.array([arm_energies(spec["kappa"], 0, 0, sign)[1] for sign in (1, -1)])
    filled = np.array([arm_energies(spec["kappa"], 1, 1, sign)[1] for sign in (1, -1)])
    require(np.allclose(filled - native, 1, atol=1e-13, rtol=0), "occupied scalar splits the endpoint clock")
    gap = native[1] - native[0]
    require(abs(gap - 1) < 1e-13, "unexpected written clock gap")
    # The control adds a local shift to the antisymmetric written state only.
    h = np.diag([0, gap + rate - 1])
    plus = np.ones(2) / math.sqrt(2)
    analyzers = np.array([plus, [1, -1j]], dtype=complex)
    analyzers[1] /= math.sqrt(2)
    output = []
    for t in elapsed:
        arms = expm(-1j * t * h) @ plus
        output.append(abs(analyzers.conj() @ arms) ** 2)
    return np.array(output)


def acquire(spec, legs):
    """The receiving station reads its clock on first hit. A miss produces a third letter."""
    endpoints = len(spec["endpoints"])
    slots = spec["clock_slots"]
    rates = spec["clock_rates"]
    shape = (2, len(rates), endpoints, endpoints, len(slots), 2, 3)
    counts = np.zeros(shape, dtype=np.int64)
    rng = np.random.Generator(np.random.PCG64([spec["seed"], 309]))
    for r, rate in enumerate(rates):
        for s, slot in enumerate(slots):
            elapsed = (slot + np.arange(1, spec["ticks"] + 1)) * spec["dt"]
            pointers = local_pointers(spec, rate, elapsed)
            for epoch, leg in enumerate(legs):
                for a in range(endpoints):
                    for b in range(endpoints):
                        if a == b:
                            continue
                        weights = leg[a, b, :-1]
                        p_plus = weights @ pointers
                        p_minus = weights @ (1 - pointers)
                        for axis in range(2):
                            p = [p_plus[axis], p_minus[axis], leg[a, b, -1]]
                            require(min(p) >= -1e-13 and abs(sum(p) - 1) < 2e-12,
                                    "invalid receiving clock law")
                            # Fixed rational acquisition keeps seeded draws portable.
                            first = np.rint(np.array(p[:2]) * 10**12) / 10**12
                            counts[epoch, r, a, b, s, axis] = rng.multinomial(
                                spec["shots"], [*first, 1 - sum(first)])
    return counts


def emission_counts(spec):
    """Fresh source-station cohorts read their own local written clock at emission."""
    probabilities = local_pointers(spec, 1, np.array(spec["clock_slots"]) * spec["dt"])
    rng = np.random.Generator(np.random.PCG64([spec["seed"], 310]))
    p = np.clip(np.rint(probabilities * 10**12) / 10**12, 0, 1)
    plus = rng.binomial(spec["shots"], p, size=(len(spec["endpoints"]), *p.shape))
    return np.stack([plus, spec["shots"]-plus], axis=-1).astype(np.int64)


def fringe(selected, shots, tail):
    """Simultaneous circles for two signed orthogonal pointer means."""
    lower = np.zeros_like(selected, dtype=float)
    upper = np.ones_like(selected, dtype=float)
    positive = selected > 0
    below = selected < shots
    lower[positive] = beta.ppf(tail, selected[positive], shots - selected[positive] + 1)
    upper[below] = beta.isf(tail, selected[below] + 1, shots - selected[below])
    center = (selected[..., 0] - selected[..., 1]) / shots
    width = np.maximum(center - (lower[..., 0] - upper[..., 1]),
                       upper[..., 0] - lower[..., 1] - center) + 1e-12
    return center[..., 0] + 1j * center[..., 1], np.hypot(width[..., 0], width[..., 1])


def score(spec, counts, emission):
    """Use only terminal clock records and launch slots. No receive-tick or carrier data."""
    e = len(spec["endpoints"])
    shots = spec["shots"]
    mask = ~np.eye(e, dtype=bool)
    require(counts.shape == (2, 2, e, e, 3, 2, 3) and counts.dtype == np.int64
            and np.all((0 <= counts) & (counts <= shots)) and np.all(counts[:, :, ~mask] == 0)
            and np.all(counts[:, :, mask].sum(axis=-1) == shots), "invalid clock counts")
    kept = counts[:, :, mask]
    # Marginal plus/minus counts are binomials even though the two counts within a trial
    # are dependent. Bonferroni covers all reported complex records, including controls.
    selected = kept[..., :2]
    require(emission.shape == (e, 3, 2, 2) and emission.dtype == np.int64 and
            np.all((0 <= emission) & (emission <= shots)) and np.all(emission.sum(axis=-1) == shots),
            "invalid sending clock counts")
    tail = spec["clock_alpha"] / (2 * (selected.size + emission.size))
    z, radius = fringe(selected, shots, tail)
    source, source_error = fringe(emission, shots, tail)
    sender = np.nonzero(mask)[0]
    report = []
    for epoch in range(2):
        rset = []
        for r in range(2):
            series = z[epoch, r]
            uncertainty = radius[epoch, r]
            enough = abs(series[:, 0]) > uncertainty[:, 0]
            results = []
            for i, j in ((0, 1), (0, 2), (1, 2)):
                delay = (spec["clock_slots"][j] - spec["clock_slots"][i]) * spec["dt"]
                bound = uncertainty[:, i] + uncertainty[:, j]
                def violations(rate):
                    residual = abs(series[:, j] - np.exp(1j * delay * rate) * series[:, i])
                    return int(np.count_nonzero(residual > bound))
                rotation = source[sender, j] * source[sender, i].conj()
                error = source_error[sender, j] + abs(source[sender, j]) * source_error[sender, i]
                observed_residual = abs(series[:, j] - rotation * series[:, i])
                results.append({"launch_slots": [spec["clock_slots"][i], spec["clock_slots"][j]],
                                "sender_clock_exclusions": int(np.sum(
                                    observed_residual > bound + abs(series[:, i]) * error)),
                                "unit_rate_exclusions": violations(1),
                                "prepared_rate_exclusions": violations(spec["clock_rates"][r])})
            rset.append({"prepared_rate": spec["clock_rates"][r],
                         "resolved_first_slot_fringes": int(enough.sum()), "comparisons": results})
        report.append(rset)
    return report
