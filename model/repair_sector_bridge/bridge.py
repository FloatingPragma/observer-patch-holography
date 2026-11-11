"""Check finite-quotient premises on actual kernels, using rational arithmetic.

The Gibbs parameter, sector energies and proposal measure are inputs. A
stationary sampler does not identify a physical Hamiltonian or clock.
"""
from __future__ import annotations

from fractions import Fraction as F


def rational(value):
    if type(value) not in (int, F):
        raise ValueError("exact integer or Fraction required")
    return F(value)


def stochastic(kernel):
    rows = tuple(tuple(rational(x) for x in row) for row in kernel)
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("nonempty square kernel required")
    if any(x < 0 for row in rows for x in row) or any(sum(row) != 1 for row in rows):
        raise ValueError("nonnegative rows of total one required")
    return rows


def normalized(weights, n):
    weights = tuple(rational(x) for x in weights)
    if len(weights) != n or any(x <= 0 for x in weights):
        raise ValueError("one positive exact weight per state required")
    return tuple(x / sum(weights) for x in weights)


def certificate(kernel, weights):
    """Exact stationarity, reversibility and a sufficient mixing condition.

Irreducibility is tested in both directions. Positive self-loops are a
sufficient aperiodicity condition; a false value does not prove periodicity.
"""
    kernel = stochastic(kernel)
    n = len(kernel)
    pi = normalized(weights, n)

    def reaches_all(transpose):
        reached, pending = {0}, [0]
        while pending:
            i = pending.pop()
            for j in range(n):
                value = kernel[j][i] if transpose else kernel[i][j]
                if value > 0 and j not in reached:
                    reached.add(j)
                    pending.append(j)
        return len(reached) == n

    stationary = all(sum(pi[i] * kernel[i][j] for i in range(n)) == pi[j] for j in range(n))
    reversible = all(pi[i] * kernel[i][j] == pi[j] * kernel[j][i] for i in range(n) for j in range(n))
    irreducible = reaches_all(False) and reaches_all(True)
    loops = all(kernel[i][i] > 0 for i in range(n))
    return {
        "stationary": stationary,
        "detailed_balance": reversible,
        "irreducible": irreducible,
        "positive_self_loops": loops,
        "unique_stationary_law": stationary and irreducible,
        "convergence_certified": stationary and irreducible and loops,
        "probabilities": [str(x) for x in pi],
    }


def lump_kernel(kernel, labels):
    """Push forward only if every representative gives the same quotient row."""
    kernel = stochastic(kernel)
    labels = tuple(labels)
    if len(labels) != len(kernel):
        raise ValueError("one quotient label per state required")
    values = tuple(dict.fromkeys(labels))
    rows = {}
    for i, label in enumerate(labels):
        row = tuple(sum(kernel[i][j] for j, target in enumerate(labels) if target == value) for value in values)
        if label in rows and rows[label] != row:
            raise ValueError(f"kernel is not lumpable at label {label!r}")
        rows[label] = row
    return values, stochastic(tuple(rows[label] for label in values))


def minorization(kernel, weights, steps):
    """Exact epsilon for K**steps(i,j) >= epsilon*pi(j).

For stationary pi, this gives TV(mu K**t,pi) <= (1-epsilon)**floor(t/steps).
The clock counts applications of this same admitted kernel.
"""
    kernel = stochastic(kernel)
    if type(steps) is not int or steps < 1:
        raise ValueError("positive integer step count required")
    n = len(kernel)
    pi = normalized(weights, n)
    if not certificate(kernel, weights)["stationary"]:
        raise ValueError("minorization target must be stationary")
    power = tuple(tuple(F(i == j) for j in range(n)) for i in range(n))
    for _ in range(steps):
        power = tuple(tuple(sum(power[i][k] * kernel[k][j] for k in range(n))
                            for j in range(n)) for i in range(n))
    return min(power[i][j] / pi[j] for i in range(n) for j in range(n))


def casimir_kernel(proposal, energies, y):
    """Casimir acceptance with complete proposed-target and rejection records.

Here y=exp(-beta/12), energies=12 C2, with rational 0<y<=1. The
stationary measure is m_i*y**energies[i] only if the proposal is reversible
for m. A different target requires the Hastings ratio in its acceptance.
"""
    proposal = stochastic(proposal)
    energies = tuple(energies)
    n = len(proposal)
    if len(energies) != n or any(type(x) is not int or x < 0 for x in energies):
        raise ValueError("one nonnegative integer energy per state required")
    y = rational(y)
    if not 0 < y <= 1:
        raise ValueError("0 < exp(-beta/12) <= 1 required")
    kernel = [[F(0) for _ in range(n)] for _ in range(n)]
    records = []
    for i in range(n):
        for j in range(n):
            acceptance = min(F(1), y ** (energies[j] - energies[i]))
            accepted = proposal[i][j] * acceptance
            rejected = proposal[i][j] - accepted
            kernel[i][j] += accepted
            kernel[i][i] += rejected
            records.append((i, j, j, "no_op" if i == j else "accepted", accepted))
            records.append((i, j, i, "rejected", rejected))
    return stochastic(kernel), tuple(records)


def su2_controls():
    """Three declared proposals on a four-label SU(2) cutoff."""
    dimensions = (1, 2, 3, 4)
    energies = tuple(3 * n * (n + 2) for n in range(4))
    y = F(9, 10)
    rows = {}
    for name, power in (("label_walk", 0), ("one_sided", 1), ("two_sided_fusion", 2)):
        proposal = [[F(0) for _ in range(4)] for _ in range(4)]
        for i in range(4):
            for j in range(4):
                if abs(i - j) != 1:
                    continue
                proposal[i][j] = (F(1, 4) if power == 0 else F(dimensions[j], 14) if power == 1
                                  else F(dimensions[j], 2 * dimensions[i]))
            proposal[i][i] = 1 - sum(proposal[i])
        kernel, records = casimir_kernel(proposal, energies, y)
        weights = [F(d) ** power * y ** e for d, e in zip(dimensions, energies)]
        row = certificate(kernel, weights)
        row["proposal_measure_power"] = power
        row["one_sided_gibbs_is_stationary"] = certificate(
            kernel, [F(d) * y ** e for d, e in zip(dimensions, energies)]
        )["stationary"]
        row["record_rows_complete"] = all(sum(r[4] for r in records if r[0] == i) == 1 for i in range(4))
        row["record_pushforward_matches_kernel"] = all(
            sum(r[4] for r in records if r[0] == i and r[2] == j) == kernel[i][j]
            for i in range(4) for j in range(4)
        )
        rows[name] = row
    identity = tuple(tuple(F(i == j) for j in range(4)) for i in range(4))
    rows["sector_preserving_control"] = certificate(identity, [1, 2, 3, 4])
    return rows


def native_quotient():
    """Extract the two-load native branch kernel and its pair-distance quotient.

    These three classical rotation orbits are not gauge irreducible sectors.
    The source's actual seam and coin outcomes, including no-ops, are counted.
    Distance is a declared full-source-state diagnostic; this function does
    not compile a distance instrument from the native current-only records.
"""
    from model.observer import coherent_continuation as source

    states, branches = source.STATES, source.branches()
    n = len(states)
    if len(branches) != n * source.EVENTS:
        raise ArithmeticError("incomplete native event register")
    kernel = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for output, *_ in branches[i * source.EVENTS:(i + 1) * source.EVENTS]:
            kernel[i][output] += F(1, source.EVENTS)
    adjacency = {i: set() for i in range(12)}
    for a, b in source.quantum_mixing.seams():
        adjacency[a].add(b)
        adjacency[b].add(a)

    def distance(a, b):
        reached, frontier = {a}, {a}
        for depth in range(12):
            if b in frontier:
                return depth
            frontier = {v for u in frontier for v in adjacency[u]} - reached
            reached |= frontier
        raise ArithmeticError("disconnected native carrier")

    pairs = [tuple(i for i, value in enumerate(state) if value) for state in states]
    if any(len(pair) != 2 for pair in pairs):
        raise ArithmeticError("two-load binary source required")
    labels = tuple(distance(*pair) for pair in pairs)
    values, quotient = lump_kernel(kernel, labels)
    sizes = tuple(labels.count(value) for value in values)
    epsilon = minorization(quotient, sizes, 2)
    return {
        "native_states": n,
        "records_per_state": source.EVENTS,
        "complete_native_records": len(branches),
        "distance_labels": list(values),
        "orbit_sizes": list(sizes),
        "kernel": [[str(x) for x in row] for row in quotient],
        "native_certificate": certificate(kernel, [1] * n),
        "quotient_certificate": certificate(quotient, sizes),
        "two_attempt_minorization": str(epsilon),
        "diagnostic_law_TV_bound_after_t_attempts": f"({1 - epsilon})**floor(t/2)",
        "interpretation": "declared full-source-state distance diagnostic; no native current-only measurement compilation or electric-sector identification",
    }


def report():
    return {
        "schema": "oph.repair-sector-bridge.v1",
        "native_quotient": native_quotient(),
        "declared_sector_controls": su2_controls(),
        "scope": "finite sampler and full-state diagnostic quotient premises only; no current-only observer instrument, physical Gibbs temperature, Hamiltonian or coupling is selected",
    }
