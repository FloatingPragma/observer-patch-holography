"""Original-denominator heralds, pair coincidences and Ramsey records on one source."""

import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from model.candidate_laws.interval_echo.openwave.photon_identity.native_emission import (
    measure as emission,
)
from model.executable_model import federation
from model.experiments.material_inertia import neutral_defect as pair
from model.experiments.material_inertia.pair_recoil.measure import PLACEMENTS

HERE = Path(__file__).resolve().parent
STARTS = 2**40
STANDARDS = (1, 32, 33)
DWELLS = (1 / 8, 1 / 4, 1 / 2)
FORMATION_TIME = 1 / 8


def native_clocks():
    """Verify all incoming ternary words, the invariant columns and their exact energies."""
    H = emission.four_h(8)
    for parity, g in ((1, 32), (-1, 33)):
        for coin in (-1, 1):
            states = np.zeros((99, 9), dtype=np.int64)
            for j, word in enumerate(emission.LAW.SRC):
                for c in (0, 1):
                    for d in (-2, 2):
                        index = 1 + emission.LAW.LABELS.index((c, d))
                        states[emission.columns(word)[index], j] = coin**c * (
                            1 if d > 0 else parity
                        )
            levels = np.array(
                [g if coin == -1 or word != (1, 1) else 0 for word in emission.LAW.SRC]
            )
            np.testing.assert_array_equal(H @ states, 4 * states * levels)


def layout(board):
    p, _ = board.external[PLACEMENTS[0]]
    u = min(
        v for v in board.neighbours[p] if board.degree[v] == 5 and v // 12 == p // 12
    )
    start = min(
        v for v in board.neighbours[u] if v != p and v not in board.neighbours[p]
    )
    # Nearest clock illuminates a different pair of terminals from the herald;
    # the second uses a seam in a disjoint host with no one-hop source access.
    available = [
        e
        for e in board.seams
        if e[0] // 12 == e[1] // 12 and u not in e and p not in e and start not in e
    ]
    near = min(e for e in available if any(z in board.neighbours[u] for z in e))
    outside = min(
        e
        for e in available
        if e[0] // 12 not in (u // 12, p // 12)
        and all(
            z not in board.neighbours[u] and z not in board.neighbours[p] for z in e
        )
    )
    return (start, p), (u, p), {"first": near, "held": outside}


def enrolled_clock(board, x, y, h, seam):
    """Replace exactly one blank seam by the measured magnitude-two clock sector."""
    a, b = seam
    left = 1 + (x == a).astype(int) - (y == a)
    right = 1 + (x == b).astype(int) - (y == b)
    hit = (left != 1) | (right != 1)
    different = np.flatnonzero(left != right)

    def exchange(address):
        return np.where(address == a, b, np.where(address == b, a, address))

    target = pair.ordered_index(exchange(x[different]), exchange(y[different]), board.n)
    diag = 8 * (left - right) ** 2 + (left != right) / 2
    blank = sparse.diags(diag) + sparse.csr_matrix(
        (np.full(len(different), -0.5), (target, different)), shape=h.shape
    )
    base = (h - blank).tocsr()
    assert (base != base.T).nnz == 0
    return base, hit


def fly(h, state, dwell):
    return expm_multiply(
        (-1j * dwell) * h, state, traceA=-1j * dwell * h.diagonal().sum()
    )


def certified_flight(h, state, dwell):
    """Scaled Taylor replay with a uniform absolute l2 error (binary64 sparse operations)."""
    steps = round(128 * dwell)
    dt = 1 / 128
    assert steps >= 1 and steps * dt == dwell
    radius = float(abs(h).sum(axis=1).max()) * dt
    # Symmetry bounds both ||H||_2 and ||abs(H)||_2 by the maximum row sum.
    max_terms = int(np.diff(h.indptr).max())
    assert radius < 1 and max_terms < 16
    degree = 16
    tail = radius ** (degree + 1) / factorial(degree + 1)
    # Each sparse dot has <16 real dyadic coefficients per row. In the
    # real/imaginary components a step of the series multiplies by i.
    # The gamma_32 allowance includes dt/k and multiplication roundoff.
    # Errors on terms sum to <=32*u*sum k*r^k/k!; adding the 17 terms
    # contributes <=32*u*exp(r). The bound 4e-14 covers both for r<1.
    per_step = tail + 4e-14
    err = 2 * steps * per_step
    result = state.astype(complex)
    for _ in range(steps):
        term = result
        advance = result.copy()
        for k in range(1, degree + 1):
            term = (-1j * dt / k) * (h @ term)
            advance += term
        result = advance
    assert err < 2e-11
    return result, err


def bounds(k, N, outcomes):
    tail = 0.01 / (2 * outcomes)
    lo = 0.0 if k == 0 else float(beta.ppf(tail, k, N - k + 1))
    hi = 1.0 if k == N else float(beta.isf(tail, k + 1, N - k))
    return [max(0.0, lo - 1e-11), min(1.0, hi + 1e-11)]


def formation_energy(board, ports, rng):
    """Sharp loads and independent swap-parity counts on every unequal seam."""
    occupation = np.ones(board.n, dtype=int)
    occupation[ports[0]], occupation[ports[1]] = 2, 0
    counts = []
    contrast = 8 * int((6 - board.degree) @ (occupation**2 - 1))
    for a, b in board.seams:
        discrepancy = int(occupation[a] - occupation[b])
        contrast += 8 * discrepancy**2
        if discrepancy:
            counts.append([int(a), int(b), int(rng.binomial(STARTS, 0.5))])
    odd = sum(row[2] / STARTS for row in counts)
    width = len(counts) * np.sqrt(np.log(2 * 100 / 0.01) / (2 * STARTS))
    return dict(
        contrast_energy=contrast,
        odd_counts=counts,
        mean_above_uniform=[contrast + odd - width, contrast + odd + width],
    )


def unconditioned_read(base, hit, initial, adjacent, rng):
    """Every incoming preparation continues directly to a hand. No formation read."""
    t = 0.25
    dark, _ = certified_flight(base, initial, t)
    result = []
    for g in (1, 32):
        bright = (
            dark
            if g == 1
            else certified_flight(
                (base + sparse.diags(g * hit, dtype=float)).tocsr(), initial, t
            )[0]
        )
        good, bad = (bright + dark) / 2, (bright - dark) / 2
        rates = np.array(
            [
                np.linalg.norm(arm[mask]) ** 2
                for arm in (good, bad)
                for mask in (adjacent, ~adjacent)
            ]
        )
        assert abs(rates.sum() - 1) < 1e-11
        counts = list(map(int, rng.multinomial(STARTS, rates / rates.sum())))
        lo, hi = bounds(sum(counts[2:]), STARTS, 100)
        result.append(
            dict(
                clock_line=g,
                dwell=t,
                bins=counts,
                wrong_clock_interval=[max(0.0, lo - 1e-10), min(1.0, hi + 1e-10)],
                pair_contact_interval=bounds(counts[0] + counts[2], STARTS, 100),
            )
        )
    return result


def run():
    native_clocks()
    board = federation(0)
    before, formed, clocks = layout(board)
    x, y, adjacent, _, h, cost = pair.generators(board)
    origin = int(pair.ordered_index(*before, board.n))
    herald = int(pair.ordered_index(*formed, board.n))
    assert (
        not adjacent[origin]
        and adjacent[herald]
        and cost[herald] == 32
        and h[herald, origin] == -0.5
    )
    state = np.zeros(len(x))
    state[origin] = 1
    launched, err = certified_flight(h, state, FORMATION_TIME)
    formed_probability = float(abs(launched[herald]) ** 2)
    assert 2 * err < 1e-10 and formed_probability > 0.0001
    np.testing.assert_allclose(
        launched[herald], fly(h, state, FORMATION_TIME)[herald], atol=1e-11, rtol=0
    )
    state[origin], state[herald] = 0, 1
    rng = np.random.default_rng(261127)
    panels = []
    # 18 clock settings, each with a herald, phase bit and pair coincidence;
    # 18 separate vacuum comparison banks are sufficient. Use 100 intervals.
    tests = 100
    for name, seam in clocks.items():
        base, hit = enrolled_clock(board, x, y, h, seam)
        assert not hit[herald]
        q = float(np.linalg.norm((base @ state)[hit]) ** 2)
        for t in DWELLS:
            dark, err = certified_flight(base, state, t)
            for g in STANDARDS:
                # For the unit clock both parities have the same base flight.
                if g == 1:
                    bright = dark
                else:
                    bright, err2 = certified_flight(
                        (base + sparse.diags(g * hit, dtype=float)).tocsr(), state, t
                    )
                    assert err2 < 2e-11
                # Relative to its separately calibrated blank phase, the wrong port
                # of the second Ramsey pulse has positive norm squared / 4.
                difference = (bright - dark) / 2
                bad = abs(difference) ** 2
                good = abs((bright + dark) / 2) ** 2
                rates = np.array(
                    [
                        v[mask].sum()
                        for v in (good, bad)
                        for mask in (adjacent, ~adjacent)
                    ]
                )
                assert abs(rates.sum() - 1) < 1e-11
                rates /= rates.sum()
                # One multinomial records nonheralds and all four coincidence bins.
                all_counts = rng.multinomial(
                    STARTS, np.r_[1 - formed_probability, formed_probability * rates]
                )
                missed, *outcomes = map(int, all_counts)
                accepted = STARTS - missed
                clock_error = sum(outcomes[2:])
                wrong = bounds(clock_error, accepted, tests)
                # The numerical nudge covers both independent vector errors and
                # accumulation of the positive sums in the wrong-port probability.
                wrong = [max(0.0, wrong[0] - 1e-10), min(1.0, wrong[1] + 1e-10)]
                panels.append(
                    dict(
                        mount=name,
                        seam=list(map(int, seam)),
                        dwell=t,
                        clock_line=g,
                        approach_rate=q,
                        nonheralds=missed,
                        accepted_bins=outcomes,
                        herald_interval=bounds(accepted, STARTS, tests),
                        wrong_clock_interval=wrong,
                        pair_contact_interval=bounds(
                            outcomes[0] + outcomes[2], accepted, tests
                        ),
                        reference_no_source_wrong=0,
                    )
                )
    base, hit = enrolled_clock(board, x, y, h, clocks["first"])
    all_starts = unconditioned_read(base, hit, launched, adjacent, rng)
    energy = {
        name: formation_energy(board, ports, rng)
        for name, ports in (("incoming", before), ("after_encounter_read", formed))
    }
    for key, idx in (("incoming", origin), ("after_encounter_read", herald)):
        lo, hi = energy[key]["mean_above_uniform"]
        assert lo < h[idx, idx] + 12 * 8 < hi
    return dict(
        unconditioned=all_starts,
        source_energy=energy,
        initial_pair=list(map(int, before)),
        detected_pair=list(map(int, formed)),
        formation_dwell=FORMATION_TIME,
        formation_probability=formed_probability,
        starts_per_setting=STARTS,
        bin_order=[
            "right hand + contact",
            "right hand + separate",
            "wrong hand + contact",
            "wrong hand + separate",
        ],
        panels=panels,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
