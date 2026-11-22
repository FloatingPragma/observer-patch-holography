"""Finite signed-offer reference on the Interval Echo carrier.

Two disjoint chains carry independent integer loads through the finite states
0, 1 and 2. A controller supplies idle, downward and upward offer multiplicities
and fresh donor loads. A saturated offer leaves a failure record. Each witness
buffers the first pair whose two chains finish. The observed word is the signed
difference of their loads. An incident cohort uses every other seam meeting the
selected seam. It accepts unanimous words. The hand's exchange-parity projection
is applied only once. The common reader compares its sharp result with one
additional independent witness: zero for even and a chart-oriented unit for odd.
The fair reader accepts even with certainty and odd with probability one half.

The banks, replenished resources, equality rule, chart and arbitration are
apparatus choices. The common comparator was proposed with the working coupling
already known. Its finite calibration returns the number of incident witnesses
for any resolved source odds. The two held populations can still fail this test.
The module returns mathematical probabilities, including exceptional events that
would need very large numbers of preparations to measure by random sampling.
"""

from collections import deque
from fractions import Fraction
from math import asin, ceil, log, log1p, sin, sqrt
from numbers import Integral

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

from model.candidate_laws.interval_echo.candidate import native_law as law
from model.executable_model import federation


def _history_counts(offers, counts, word):
    if isinstance(offers, bool) or not isinstance(offers, Integral) or offers < 0:
        raise ValueError("offers must be a nonnegative integer")
    if any(
        len(weights) != 3
        or not all(isinstance(n, Integral) and n >= 0 for n in weights)
        or not sum(weights)
        for weights in (counts, word)
    ):
        raise ValueError(
            "counts and word each require three nonnegative integers with positive sum"
        )
    port_counts = tuple(map(int, counts))
    idle, down, up = map(int, word)
    denominator = sum(port_counts)
    first_failure = []
    for _ in range(offers):
        low, middle, high = port_counts
        denominator *= idle + down + up
        first_failure.append(Fraction(down * low + up * high, denominator))
        port_counts = (
            idle * low + down * middle,
            up * low + idle * middle + down * high,
            up * middle + idle * high,
        )
    survive = Fraction(sum(port_counts), denominator)
    if sum(first_failure) + survive != 1 or not survive:
        raise ValueError("Invalid or empty surviving source")
    return port_counts, survive, first_failure


def signed_probabilities(offers, counts=(1, 1, 1), word=(2, 1, 1)):
    """Return exact (P(d=0), P(d=+1), P(d=+2)) for two independent surviving chains.

    Each negative orientation has its positive partner's probability even for
    asymmetric offers, since the endpoints use the same distribution.
    """
    loads, _, _ = _history_counts(offers, counts, word)
    a, b, c = loads
    denom = sum(loads) ** 2
    return (
        Fraction(a * a + b * b + c * c, denom),
        Fraction(a * b + b * c, denom),
        Fraction(a * c, denom),
    )


def incident_witnesses(board, seam):
    a, b = sorted(seam)
    if b not in board.neighbours[a]:
        raise ValueError("The selected endpoints have no native seam")
    return sorted(
        {tuple(sorted((u, v))) for u in (a, b) for v in board.neighbours[u]} - {(a, b)}
    )


def finite_service(
    offers, witnesses, *, reader="common", cap=None, counts=(1, 1, 1), word=(2, 1, 1)
):
    """One local epoch, starting with empty witness buffers and a protected current mark.

    A capped epoch stops at its first overflowing buffer. The cost includes all
    failed and successful pair preparations up to that outcome.
    """
    if not isinstance(witnesses, Integral) or witnesses < 1:
        raise ValueError("At least one witness is required")
    if cap is not None and (not isinstance(cap, Integral) or cap < 1):
        raise ValueError("cap must be a positive integer or None")
    p0, p1, p2 = signed_probabilities(offers, counts, word)
    admissions = {
        "common": (p0, p1),
        "common_reversed_chart": (p0, p1),
        "common_sign_blind": (p0, 2 * p1),
        "common_inverted": (p1, p0),
        "fair": (Fraction(1), Fraction(1, 2)),
    }
    if reader not in admissions:
        raise ValueError("Unknown parity reader")
    even, odd = admissions[reader]
    if min(even, odd, p0, p1, p2) <= 0 or even == odd:
        raise ValueError(
            "The five-mode calibration needs full support and a resolved parity split"
        )
    m = int(witnesses)
    buffers = m + int(reader != "fair")
    agreement = p0**m + 2 * p1**m + 2 * p2**m
    weights = (p0**m * even, p1**m * even, p1**m * odd, p2**m * even, p2**m * odd)
    pointer_when_filled = sum(weights)
    modes = tuple(p / pointer_when_filled for p in weights)
    z, one_even, one_odd, two_even, two_odd = modes
    _, survives, first_failure = _history_counts(offers, counts, word)
    pair_success = survives**2
    # Floating complements avoid huge rational powers at practical retry caps.
    fill_one = (
        1.0
        if cap is None or pair_success == 1
        else -np.expm1(cap * log1p(-float(pair_success)))
    )
    fill_all = fill_one**buffers
    chance = fill_all * float(pointer_when_filled)
    pairs_per_epoch = (
        fill_one / float(pair_success) * sum(fill_one**j for j in range(buffers))
    )
    beta = log(even / odd)
    held = (two_even * z**3 / one_even**4, two_odd * z**3 / (one_even**3 * one_odd))
    if (
        reader in ("common", "common_reversed_chart")
        and z / one_even != (one_even / one_odd) ** m
    ):
        raise ValueError("The common odds calibration failed")
    return dict(
        witnesses=m,
        buffers=buffers,
        signed_probabilities=(p0, p1, p2),
        first_failure_by_offer=first_failure,
        pair_source_success=pair_success,
        signed_agreement_given_sources=agreement,
        modes=modes,
        beta=beta,
        inferred_kappa=log(z / one_even) / beta,
        held_ratios=held,
        conditional_on_filling_pointer_flip=pointer_when_filled,
        overflow=1 - fill_all,
        disagreement=fill_all * (1 - float(agreement)),
        parity_reject=fill_all * float(agreement - pointer_when_filled),
        pointer_flip=chance,
        contraction=1 - chance,
        trace_distance_1pct_epochs=ceil(log(0.01) / log1p(-chance)),
        mean_pairs_per_epoch=pairs_per_epoch,
        mean_pairs_per_flip=pairs_per_epoch / chance,
        mean_pairs_per_held_odd=pairs_per_epoch / (chance * float(two_odd)),
        mean_pairs_per_either_held=pairs_per_epoch
        / (chance * float(two_even + two_odd)),
    )


def check_native(offers=17, witnesses=8):
    """Verify the material source gates, signed comparison, parity, and bank pointer."""
    blank = np.eye(law.DIM)[:, : law.N_SRC]
    echo = law.echo()
    written = echo @ blank
    exchange, contrast = law.local_term_parts()
    swap = np.eye(law.N_SRC)[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    signed_label = np.repeat([-9, *(d for _, d in law.LABELS)], law.N_SRC)
    for d in range(-2, 3):
        expected = np.diag([a - b == d for a, b in law.SRC])
        projected = (signed_label == d)[:, None] * written
        np.testing.assert_allclose(echo @ projected, blank @ expected, atol=1e-13)
    equality = np.rint(
        np.diag(written.T @ ((signed_label == 0)[:, None] * written))
    ).astype(int)
    for a in range(-2, 3):
        for b in range(-2, 3):
            equal = (
                equality[law.SRC.index((abs(a), abs(b)))]
                * equality[law.SRC.index((int(a < 0), int(b < 0)))]
            )
            if bool(equal) != (a == b):
                raise ValueError("The sharp signed-record comparison failed")
    loads, _, _ = _history_counts(offers, (1, 1, 1), (2, 1, 1))
    p0, p1, p2 = signed_probabilities(offers)
    sign_p = (p2, p1, p0, p1, p2)
    joint = np.array(
        [
            float(
                Fraction(loads[a] * loads[b], sum(loads) ** 2)
                * sign_p[a - b + 2] ** (witnesses - 1)
            )
            for a, b in law.SRC
        ]
    )
    rho = np.diag(joint / joint.sum())
    oriented_read = np.tile([a <= b for a, b in law.SRC], 1 + law.N_LAB)
    magnitudes = np.tile([abs(a - b) for a, b in law.SRC], 1 + law.N_LAB)
    reads = [magnitudes == 0]
    for m in (1, 2):
        reads.extend(
            [(magnitudes == m) & oriented_read, (magnitudes == m) & ~oriented_read]
        )
    theta = np.diag([1j if a < b else 1 for a, b in law.SRC])
    sign_phase = np.where(signed_label < 0, 1j, 1)
    sign_phase[: law.N_SRC] = 1
    for kappa in (0, 0.5, 1.3, 4.5, 8):
        half = expm(-0.5j * np.pi * (exchange + kappa * contrast))
        full = half @ half
        # Valid non-idle offers use a native repair at magnitude two and a
        # nondestructive comparison at one. Zero records a saturated offer.
        for n in range(3):
            for donor, step in ((0, -1), (2, 1)):
                incoming = blank[:, law.SRC.index((n, donor))]
                mag = abs(n - donor)
                measured = (abs(signed_label) == mag) * (echo @ incoming)
                np.testing.assert_allclose(measured, echo @ incoming, atol=1e-13)
                next_state = measured if mag == 2 else echo @ measured
                if 0 <= n + step <= 2:
                    wrong_next = np.tile(
                        [b != n + step for a, b in law.SRC], 1 + law.N_LAB
                    )
                    np.testing.assert_allclose(next_state[wrong_next], 0, atol=1e-13)
                elif mag != 0:
                    raise ValueError(
                        "Failed boundary offer has a nonzero written magnitude"
                    )
            # Idle transport, disposal, sharp donor loading, and the one-unit
            # current-bank pointer all use this full-contact identity.
            left, right = law.SRC.index((n, 0)), law.SRC.index((0, n))
            phase = np.exp(-1j * np.pi * kappa * n * n)
            np.testing.assert_allclose(
                full @ blank[:, [left, right]],
                phase * blank[:, [right, left]],
                atol=1e-12,
            )
        forward = half @ echo @ (sign_phase[:, None] * written)
        rewind = theta.conj().T @ blank.T @ np.linalg.matrix_power(half, 3)
        maps = []
        for orientation, sign in ((oriented_read, 1), (~oriented_read, -1)):
            op = rewind @ (orientation[:, None] * forward)
            parity = (np.eye(law.N_SRC) + sign * swap) / 2
            np.testing.assert_allclose(op.conj().T @ op, parity, atol=1e-12)
            np.testing.assert_allclose(swap @ op, sign * op, atol=1e-12)
            for a, b in ((0, 1), (1, 2), (0, 2)):
                state = (
                    np.eye(law.N_SRC)[:, law.SRC.index((a, b))]
                    + sign * np.eye(law.N_SRC)[:, law.SRC.index((b, a))]
                ) / np.sqrt(2)
                projected = op @ state
                np.testing.assert_allclose(
                    np.outer(projected, projected.conj()),
                    np.outer(state, state),
                    atol=1e-12,
                )
            maps.append(op)
        # For a reversible proposal, prepare the chosen orthogonal eigenmode
        # on a scratch bank with a half contact and a sign-selective phase.
        negative = signed_label < 0
        negative[: law.N_SRC] = False
        for m in (1, 2):
            left, right = (blank[:, law.SRC.index(pair)] for pair in ((m, 0), (0, m)))
            for odd in (False, True):
                pulse = np.where(negative, 1j if odd else -1j, 1)
                actual = pulse * (echo @ half @ left)
                expected = echo @ ((left + (-1 if odd else 1) * right) / np.sqrt(2))
                np.testing.assert_allclose(
                    np.outer(actual, actual.conj()),
                    np.outer(expected, expected.conj()),
                    atol=1e-12,
                )
        for reader, admit in (
            ("common", (p0, p1)),
            ("fair", (Fraction(1), Fraction(1, 2))),
        ):
            record = finite_service(offers, witnesses, reader=reader)
            after = sum(
                float(p) * (op @ rho @ op.conj().T)
                for p, op in zip(admit, maps, strict=True)
            )
            np.testing.assert_allclose(
                after.trace().real,
                float(
                    record["conditional_on_filling_pointer_flip"]
                    / record["signed_agreement_given_sources"]
                ),
                atol=1e-12,
            )
            after /= after.trace()
            accepted = np.array([float(q) for q in record["modes"]])
            material = written[law.N_SRC :] @ after @ written[law.N_SRC :].T
            memory = np.trace(
                material.reshape(law.N_LAB, law.N_SRC, law.N_LAB, law.N_SRC),
                axis1=1,
                axis2=3,
            )
            for coin in (0, 1):
                ids = [law.LABELS.index((coin, d)) for d in range(-2, 3)]
                block = 2 * memory[np.ix_(ids, ids)]
                for observable in (np.eye(5)[::-1], np.diag([4, 1, 0, 1, 4])):
                    np.testing.assert_allclose(
                        block @ observable, observable @ block, atol=1e-12
                    )
                np.testing.assert_allclose(
                    np.linalg.eigvalsh(block), sorted(accepted), atol=1e-12
                )
            final = forward @ after @ forward.conj().T
            np.testing.assert_allclose(
                [np.diag(final)[r].sum().real for r in reads], accepted, atol=1e-12
            )
            # Dephasing the hand's signed write erases the parity distinction.
            hand = written @ after @ written.T
            hand[signed_label[:, None] != signed_label[None, :]] = 0
            decoder = half @ echo @ np.diag(sign_phase)
            washed = np.diag(decoder @ hand @ decoder.conj().T).real
            bins = [washed[r].sum() for r in reads]
            np.testing.assert_allclose(
                [bins[1], bins[3]], [bins[2], bins[4]], atol=1e-12
            )
            # With independent replacement, tracing the unselected bank gives
            # (1-success)*old + success*new for either initial address.
            target = np.diag(accepted)
            old = np.ones((5, 5)) / 5
            success = record["pointer_flip"]
            for old_address in (0, 1):
                tensor = (
                    np.kron(old, target) if old_address == 0 else np.kron(target, old)
                ).reshape(5, 5, 5, 5)
                banks = (
                    np.trace(tensor, axis1=1, axis2=3),
                    np.trace(tensor, axis1=0, axis2=2),
                )
                result = (1 - success) * banks[old_address] + success * banks[
                    1 - old_address
                ]
                np.testing.assert_allclose(
                    result, (1 - success) * old + success * target, atol=1e-13
                )


def _cell_path(board, start, size):
    def extend(path):
        if len(path) == size:
            return path
        for v in board.neighbours[path[-1]]:
            if v // 12 == start // 12 and v not in path:
                result = extend([*path, v])
                if result:
                    return result

    return extend([start])


def _source_chains(board, count, offers, excluded_cells):
    """Allocate simple, mutually disjoint paths of up to 24 vertices."""
    chains = []
    if offers < 12:
        for cell in range(board.cells):
            if cell not in excluded_cells:
                chains.append(_cell_path(board, 12 * cell, offers + 1))
                if len(chains) == count:
                    break
    else:
        used = set(excluded_cells)
        first = (offers + 1) // 2
        for a, b in board.external:
            if a // 12 not in used and b // 12 not in used:
                left = _cell_path(board, a, first)
                right = _cell_path(board, b, offers + 1 - first)
                if left and right:
                    used.update((a // 12, b // 12))
                    chains.append([*reversed(left), *right])
                    if len(chains) == count:
                        break
    if len(chains) != count or not all(chains):
        raise ValueError("The board lacks the requested disjoint source chains")
    if len({v for chain in chains for v in chain}) != (offers + 1) * count:
        raise ValueError("The source preparations overlap")
    return chains


def _transport(board, start, end, protected, pending):
    parent = {start: None}
    queue = deque([start])
    while queue:
        port = queue.popleft()
        if port == end:
            path = []
            while port is not None:
                path.append(port)
                port = parent[port]
            return path[::-1]
        for other in board.neighbours[port]:
            if (
                other not in parent
                and other not in pending
                and tuple(sorted((port, other))) not in protected
            ):
                parent[other] = port
                queue.append(other)
    raise ValueError("No blank path connects the source to its witness")


def physical_layout(board, seam, offers, *, reader="common"):
    """Return a sufficient native layout for one cohort, including protected banks.

    Source chains and written collars are held out of the transfer service.
    Remote full swaps on blank paths move sharp inputs and restore spectators.
    Fresh preparation/disposal trits and feedback are supplied on later epochs.
    """
    if not 0 <= offers < 24 or reader not in ("common", "fair"):
        raise ValueError(
            "This layout supports 0–23 offers and either the common or fair reader"
        )
    target = tuple(sorted(seam))
    witnesses = incident_witnesses(board, target)
    stations = witnesses.copy()
    if reader == "common":
        touched = {p for edge in [target, *witnesses] for p in edge}
        reporter = next(
            (
                e
                for e in board.seams
                if not touched.intersection(e) and e[0] // 12 == target[0] // 12
            ),
            None,
        )
        if reporter is None:
            raise ValueError("No independent local parity reporter")
        stations.append(reporter)
    used_cells = {v // 12 for edge in [target, *stations] for v in edge}
    chains = _source_chains(board, 2 * len(stations), offers, used_cells)
    protected = {
        tuple(sorted((a, b))) for chain in chains for a, b in zip(chain, chain[1:])
    }
    pending = {chain[-1] for chain in chains}
    occupied = {v for chain in chains for v in chain}
    moves = []
    for i, edge in enumerate(stations):
        legs = []
        for k, port in enumerate(edge):
            source = chains[2 * i + k][-1]
            pending.remove(source)
            path = _transport(board, source, port, protected, pending)
            # The reverse leg omits the last swap. Any spectator is restored.
            tokens = path.copy()
            for j in (*range(len(path) - 1), *range(len(path) - 3, -1, -1)):
                tokens[j], tokens[j + 1] = tokens[j + 1], tokens[j]
            if tokens != [path[-1], *path[1:-1], path[0]]:
                raise ValueError("Remote transfer failed to restore a spectator")
            legs.append(path)
            occupied.update(path)
        protected.add(edge)
        moves.append(legs)
    comparison = []
    for edge in board.seams:
        if not occupied.intersection(edge):
            comparison.append(edge)
            occupied.update(edge)
            if len(comparison) == 3:
                break
    banks = []
    for a, b in board.seams:
        if a in occupied or b in occupied:
            continue
        first = next(
            (v for v in board.neighbours[a] if v != b and v not in occupied), None
        )
        second = next(
            (
                v
                for v in board.neighbours[b]
                if v not in (a, first) and v not in occupied
            ),
            None,
        )
        if first is not None and second is not None:
            banks.append(((a, b), (first, second)))
            occupied.update((a, b, first, second))
            if len(banks) == 2:
                break
    pointers = []
    for edge in board.seams:
        if not occupied.intersection(edge):
            pointers.append(edge)
            occupied.update(edge)
            if len(pointers) == 2:
                break
    if tuple(map(len, (comparison, banks, pointers))) != (3, 2, 2):
        raise ValueError("Missing comparator, bank/disposal, or pointer rails")
    return dict(
        seam=target,
        witnesses=witnesses,
        parity_reporter=stations[-1] if reader == "common" else None,
        chains=chains,
        transfer_paths=moves,
        comparator_rails=comparison,
        banks=banks,
        pointers=pointers,
        source_ports=sum(map(len, chains)),
        laboratory_cells=sorted({v // 12 for v in occupied}),
    )


def exact_stationary_source(offers=17, witnesses=8):
    """Calibrate the ternary preparation by one valid offer, before reading collars.

    Let (a,b,a) be the endpoint load probabilities and r=b/a. The unstopped
    valid-offer transfer has weights (2a+b,2a+2b,2a+b)/4. Equality of the
    normalized ternary distributions requires r=(2+2r)/(2+r). Its unique
    positive solution is sqrt(2). On the monotone interval [0,pi], the native
    isolated one-load contact supplies b=sin(theta/2)^2. A separate fair
    binary decision assigns each remaining preparation to load 0 or 2.

    The analog calibration uses only pre- and post-offer load frequencies. It
    assumes a continuously tunable exposure, separate donor initialization and
    protected measurements. Digital finite-depth preparations in finite_service
    need no such exposure. The conditional five-mode spectrum belongs to an
    isolated written collar after evacuating both load ports into vacancies.
    The native echo extension on all 99 states is still part of the law.
    """
    if (
        not isinstance(offers, Integral)
        or offers < 0
        or not isinstance(witnesses, Integral)
        or witnesses < 1
    ):
        raise ValueError("Require a nonnegative offer depth and positive witness count")
    import sympy as sp

    root = sp.sqrt(2)
    a = 1 / (2 + root)
    b = root * a
    survival = (2 + root) / 4
    zero, unit, double = 2 * a * a + b * b, 2 * a * b, a * a
    transformed = (
        sp.Matrix([[2, 1, 0], [1, 2, 1], [0, 1, 2]]) * sp.Matrix([a, b, a]) / 4
    )
    if any(
        sp.simplify(x - y) != 0
        for x, y in zip(transformed, survival * sp.Matrix([a, b, a]))
    ):
        raise ValueError("The ternary preparation changed after a valid offer")
    identities = (
        2 * a + b - 1,
        b - (root - 1),
        zero + 2 * unit + 2 * double - 1,
        zero - root * unit,
        zero - 4 * double,
    )
    if any(sp.simplify(x) != 0 for x in identities):
        raise ValueError("Exact oriented source weights failed")
    m = int(witnesses)
    # Write up/down = r². The killed transfer below acts on column masses.
    r, down = sp.symbols("r down", positive=True)
    idle = sp.symbols("idle", nonnegative=True)
    up = down * r**2
    killed = sp.Matrix([[idle, down, 0], [up, idle, down], [0, up, idle]])
    masses = sp.Matrix([1, root * r, r**2])
    residual = killed * masses - (idle + root * down * r) * masses
    if any(sp.simplify(value) != 0 for value in residual):
        raise ValueError("The proposed biased source is not conditionally stationary")
    same_weight = masses.dot(masses)
    unit_weight = masses[1] * (masses[0] + masses[2])
    double_weight = masses[0] * masses[2]
    # Parity cancels for either held mode. This checks the whole family
    # against the event frequencies before substituting any numerical rates.
    multiplier = sp.factor(double_weight * same_weight**3 / unit_weight**4)
    if sp.simplify(multiplier - 1 - (up - down)**2 / (4 * up * down)) != 0:
        raise ValueError("The stationary held-ratio factorization failed")
    ports = np.array([float(a), float(b), float(a)])
    probs = np.array([float(zero), float(unit), float(double)])
    beta = log(2) / 2
    theta = 2 * asin(sqrt(float(b)))
    levels = [0, m, m + 1, 4 * m, 4 * m + 1]
    boltzmann = np.exp(-beta * np.array(levels))
    # Read the modes from the physical sign consensus and parity reporter;
    # evaluate the quadratic Gibbs spectrum separately as a check.
    p0, p1, p2 = probs
    raw_modes = np.array(
        [p0 ** (m + 1), p1**m * p0, p1 ** (m + 1), p2**m * p0, p2**m * p1]
    )
    modes = raw_modes / raw_modes.sum()
    np.testing.assert_allclose(modes, boltzmann / boltzmann.sum(), atol=2e-14)

    def transfer(q, word):
        idle, down, up = word
        low, middle, high = q
        raw = np.array(
            [
                idle * low + down * middle,
                up * low + idle * middle + down * high,
                up * middle + idle * high,
            ]
        ) / sum(word)
        return raw.sum(), raw / raw.sum()

    def calibration_residual(angle):
        middle = sin(angle / 2) ** 2
        q = np.array([(1 - middle) / 2, middle, (1 - middle) / 2])
        _, after = transfer(q, (2, 1, 1))
        return after[1] - middle

    calibrated = brentq(calibration_residual, 0, np.pi, xtol=1e-14)
    np.testing.assert_allclose(calibrated, theta, atol=2e-14, rtol=0)

    def control(q, word):
        rate, after = transfer(q, word)
        # Compute the unspent double-mode prediction from a fresh pair of endpoints.
        low, middle, high = after
        z = low * low + middle * middle + high * high
        u = middle * (low + high)
        v = low * high
        return dict(
            start_loads=list(map(float, q)),
            after_one_valid_offer=list(map(float, after)),
            max_load_probability_change=float(max(abs(after - q))),
            one_offer_survival=float(rate),
            held_ratio_after_one_offer=float((v * z**3 / u**4) ** m),
        )

    blank = np.eye(law.DIM)[:, : law.N_SRC]
    echo = law.echo()
    h0, h1 = law.local_term_parts()
    # Construct the actually written state by treating the last witness's two
    # distinct loads separately, then trace the original loads after evacuation.
    counts = np.array(
        [ports[i] * ports[j] * probs[abs(i - j)] ** (m - 1) for i, j in law.SRC]
    )
    source = np.diag(counts / counts.sum())
    swap = np.eye(law.N_SRC)[:, [law.SRC.index((j, i)) for i, j in law.SRC]]
    even = (np.eye(law.N_SRC) + swap) / 2
    odd = (np.eye(law.N_SRC) - swap) / 2
    admitted = probs[0] * (even @ source @ even) + probs[1] * (odd @ source @ odd)
    admitted /= admitted.trace()
    W = echo[law.N_SRC :, : law.N_SRC]
    written = W @ admitted @ W.T
    memory = np.trace(
        written.reshape(law.N_LAB, law.N_SRC, law.N_LAB, law.N_SRC), axis1=1, axis2=3
    )
    unit_load = blank[:, law.SRC.index((1, 0))]
    moved = blank[:, law.SRC.index((0, 1))]
    pulse_controls = []
    for k in dict.fromkeys((-8, 0, 1.5, 8, m)):
        h = h0 + k * h1
        prepared = expm(-1j * theta * h) @ unit_load
        expected = (
            np.exp(-1j * k * theta)
            * (
                (1 + np.exp(-1j * theta)) * unit_load
                + (1 - np.exp(-1j * theta)) * moved
            )
            / 2
        )
        np.testing.assert_allclose(prepared, expected, atol=2e-13)
        pulse_controls.append(
            dict(
                generator_kappa=k,
                middle_load_probability=float(abs(np.vdot(moved, prepared)) ** 2),
            )
        )
        # A full native contact to a vacant disposal port evacuates any load;
        # its phase travels with the discarded load. W itself is a null sector.
        clear = expm(-1j * np.pi * h)
        for n in range(3):
            np.testing.assert_allclose(
                clear @ blank[:, law.SRC.index((n, 0))],
                np.exp(-1j * np.pi * k * n * n) * blank[:, law.SRC.index((0, n))],
                atol=2e-12,
            )
        np.testing.assert_allclose(h @ echo @ blank, 0, atol=2e-13)
        for coin in (0, 1):
            ids = [law.LABELS.index((coin, d)) for d in range(-2, 3)]
            empty_columns = np.eye(law.DIM)[
                :, [(1 + i) * law.N_SRC + law.SRC.index((0, 0)) for i in ids]
            ]
            projected = (np.eye(5) - np.eye(5)[::-1]) / 2 + k * np.diag([4, 1, 0, 1, 4])
            np.testing.assert_allclose(
                h @ empty_columns, empty_columns @ projected, atol=2e-13
            )
            if k == m:
                gibbs = expm(-beta * projected)
                np.testing.assert_allclose(
                    2 * memory[np.ix_(ids, ids)], gibbs / gibbs.trace(), atol=2e-13
                )
    wrong_angle = []
    for delta in (-0.05, 0.05):
        middle = sin((theta + delta) / 2) ** 2
        wrong_angle.append(
            dict(
                offset_radians=delta,
                **control(
                    np.array([(1 - middle) / 2, middle, (1 - middle) / 2]), (2, 1, 1)
                ),
            )
        )
    one = float(survival)
    pair = one ** (2 * int(offers))
    first_failure = [(1 - one) * one**n for n in range(offers)]
    np.testing.assert_allclose(1, one**offers + sum(first_failure), atol=2e-14)
    raw_odd = raw_modes[4]
    agreement = p0**m + 2 * p1**m + 2 * p2**m
    accepted = raw_modes.sum()
    filled = pair ** (m + 1)
    batch = dict(
        any_source_failed=1 - filled,
        unequal_signed_records=filled * (1 - agreement),
        parity_reject=filled * (agreement - accepted),
        accepted=filled * accepted,
    )
    np.testing.assert_allclose(sum(batch.values()), 1, atol=2e-14)
    return dict(
        calibrated_pulse_angle=theta,
        independent_one_offer_root=calibrated,
        endpoint_probabilities=ports.tolist(),
        signed_probabilities=probs.tolist(),
        per_window_survival=one,
        first_failure_one_chain=first_failure,
        pair_survival=pair,
        beta=beta,
        native_energy_levels=levels,
        held_gap_units=[4 * m, 4 * m + 1],
        mode_probabilities=modes.tolist(),
        even_magnitude_gaps=[m, 4 * m],
        held_ratios=[1, 1],
        stationary_offer_family=dict(
            rates="Positive state-independent up and down probabilities p and q. Invalid boundary offers reject.",
            conditional_endpoint_weights="1, sqrt(2*p/q), p/q",
            held_even_and_odd_multiplier="(1 + (p-q)^2/(4*p*q))^m",
            one_to_one=float(multiplier.subs(r, 1)**m),
            two_to_one=float(multiplier.subs(r, root)**m),
        ),
        mean_raw_pairs_per_held_odd=(m + 1) / (pair * raw_odd),
        one_batch_branches=batch,
        buffered_epoch_branches=dict(
            unequal_signed_records=1 - agreement,
            parity_reject=agreement - accepted,
            accepted=accepted,
        ),
        pulse_controls=pulse_controls,
        wrong_angle_controls=wrong_angle,
        biased_offer_control=control(ports, (1, 1, 2)),
    )


def reversible_diagnostics(
    offers, witnesses, cap=None, counts=(1, 1, 1), word=(2, 1, 1)
):
    """A sharp five-mode bank updated through local signed-record competitions.

    Three fair material bits propose an idle with probability one half, or one
    of the four other eigenmodes with probability one eighth each. For each
    coordinate that differs, compare fresh pair histories until either requested
    signed word appears or a local cap expires. A proposed winner advances to
    the next comparison. A current winner or timeout leaves the bank alone.
    Magnitude uses every incident witness; parity uses one identical service.
    Equivalent coordinates require no comparison. No transition consults a
    target probability. Sharp proposals and freshly initialized trits are
    supplied resources. Costs count raw pair preparations before resetting
    sources or retiring old banks; they exclude reset work and routing gates.
    """
    reference = finite_service(offers, witnesses, counts=counts, word=word)
    if cap is not None and (not isinstance(cap, Integral) or cap < 1):
        raise ValueError("The per-comparison cap must be a positive integer")
    probs = reference["signed_probabilities"]
    pair_success = float(reference["pair_source_success"])
    target = np.array(list(map(float, reference["modes"])))
    modes = ((0, 0), (1, 0), (1, 1), (2, 0), (2, 1))
    transitions = np.zeros((5, 5))
    proposal_cost = np.zeros((5, 5))
    failures = np.zeros((5, 2))  # lost to the current word, or timed out
    comparisons = {}
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            hit = pair_success * float(probs[i] + probs[j])
            resolved = 1.0 if cap is None or hit == 1 else -np.expm1(cap * log1p(-hit))
            wins = resolved * float(probs[j] / (probs[i] + probs[j]))
            loses = resolved * float(probs[i] / (probs[i] + probs[j]))
            comparisons[i, j] = (wins, loses, 1 - resolved, resolved / hit)
    for i, (mi, pi) in enumerate(modes):
        for j, (mj, pj) in enumerate(modes):
            if i == j:
                continue
            # For a cap, the chance of resolving by the deadline is symmetric
            # in the two competing words. It cancels in forward/reverse odds.
            exact_odds = (probs[mj] / probs[mi]) ** witnesses * probs[pj] / probs[pi]
            if exact_odds != reference["modes"][j] / reference["modes"][i]:
                raise ValueError("Reversible target ratio failed")
            gates = ([(mi, mj)] * witnesses if mi != mj else []) + (
                [(pi, pj)] if pi != pj else []
            )
            pending = 1 / 8
            for pair in gates:
                wins, loses, timeout, mean_attempts = comparisons[pair]
                proposal_cost[i, j] += pending * mean_attempts
                failures[i] += pending * np.array([loses, timeout])
                pending *= wins
            transitions[i, j] = pending
        transitions[i, i] = 0.5 + failures[i].sum()
    np.testing.assert_allclose(transitions.sum(axis=1), 1, atol=2e-14)
    flux = target[:, None] * transitions
    np.testing.assert_allclose(flux, flux.T, atol=2e-16, rtol=1e-12)
    eigenvalues = np.linalg.eigvalsh(
        transitions * np.sqrt(target[:, None] / target[None, :])
    )
    gap = 1 - eigenvalues[-2]
    upper = first = None
    # Below floating point resolution, omit an unreliable numerical mixing time.
    if gap > 1e-14:
        upper = ceil(log(0.02 / np.sqrt(1 / target.min() - 1)) / log1p(-gap))
        low, high = 0, upper
        while high - low > 1:
            mid = (high + low) // 2
            distance = (
                np.abs(np.linalg.matrix_power(transitions, mid) - target)
                .sum(axis=1)
                .max()
                / 2
            )
            if distance <= 0.01:
                high = mid
            else:
                low = mid
        first = high
    changes = 1 - transitions.diagonal()
    branches = np.column_stack((np.full(5, 0.5), failures, changes))
    avg = float(target @ proposal_cost.sum(axis=1))
    flips = float(target @ changes)
    entries = float(flux[:3, 3:].sum())
    odd_entries = float(flux[:4, 4].sum())
    return dict(
        witnesses=witnesses,
        offers=offers,
        comparison_cap=cap,
        modes=reference["modes"],
        transition=transitions.tolist(),
        branch_order=[
            "idle",
            "current_word_wins",
            "timeout",
            "proposed_mode_committed",
        ],
        branches_by_current_mode=branches.tolist(),
        mean_pairs_by_current_mode=proposal_cost.sum(axis=1).tolist(),
        detailed_balance_error=float(abs(flux - flux.T).max()),
        gap=float(gap) if upper is not None else None,
        tv_1pct_spectral_ceiling=upper,
        tv_1pct_all_point_starts_first_epoch=first,
        mean_raw_pair_attempts_per_epoch=avg,
        mean_raw_pair_attempts_per_change=avg / flips,
        mean_raw_pair_attempts_per_entry_into_magnitude_two=avg / entries,
        mean_raw_pair_attempts_per_entry_into_odd_magnitude_two=avg / odd_entries,
    )


def replay():
    """Run the native controls and return compact independently replayable statistics."""
    board = federation(1)
    check_native()
    layouts = []
    targets = ((0, 5), (0, 7), (0, 1), (1, 5), board.external[0], board.external[1])
    for seam in targets:
        for reader in ("common", "fair"):
            layout = physical_layout(board, seam, 17, reader=reader)
            layouts.append(
                dict(
                    seam=seam,
                    reader=reader,
                    witnesses=len(layout["witnesses"]),
                    source_ports=layout["source_ports"],
                    laboratory_cells=len(layout["laboratory_cells"]),
                )
            )
    if [row["witnesses"] for row in layouts[::2]] != [8, 8, 9, 9, 10, 10]:
        raise ValueError("The enrolled incidences changed")
    common = finite_service(8, 8)
    reversed_chart = finite_service(8, 8, reader="common_reversed_chart")
    reversed_words = finite_service(8, 8, reader="common_inverted")
    blind = finite_service(8, 8, reader="common_sign_blind")
    biased = finite_service(8, 8, word=(1, 1, 2))
    capped = finite_service(8, 8, cap=134)
    if (
        common["modes"] != reversed_chart["modes"]
        or reversed_words["beta"] >= 0
        or blind["beta"] >= 0
        or biased["held_ratios"][0] <= 2
        or capped["modes"] != common["modes"]
        or not 0 < capped["overflow"] < 0.001
    ):
        raise ValueError("A reader, preparation or buffer control failed")
    return layouts


if __name__ == "__main__":
    for row in replay():
        print(row)
    for reader in ("common", "fair"):
        row = finite_service(17, 8, reader=reader)
        print(
            reader,
            "kappa",
            row["inferred_kappa"],
            "held errors",
            [float(v - 1) for v in row["held_ratios"]],
            "mean pair attempts per held odd outcome",
            row["mean_pairs_per_held_odd"],
        )
