"""Compile a material pulse from signed-word changes and unit parity rails."""
import json
from itertools import product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import cell, federation


def signed_table():
    """Decode the actual two coin branches of each blank-collar write."""
    weights = (abs(law.write().reshape(law.N_LAB, law.N_SRC, law.N_SRC))**2).sum(axis=1)
    words = np.empty((3, 3), dtype=int)
    for column, (a, b) in enumerate(law.SRC):
        active = np.flatnonzero(weights[:, column])
        if len(active) != 2:
            raise ValueError("A native write must produce exactly two coin outcomes")
        first, second = (law.LABELS[i] for i in active)
        if first != (0, a-b) or second != (1, a-b):
            raise ValueError("Coin branches disagree with the signed endpoint record")
        np.testing.assert_allclose(weights[active, column], .5, rtol=0, atol=3e-16)
        words[a, b] = first[1]
    return words


def addresses(board, edge):
    return sorted({tuple(sorted((v, u))) for v in edge for u in board.neighbours[v]})


def changes(readers, edge, loads, words, transport=True):
    """Compute signed changes by re-reading the actual endpoint values."""
    a, b = edge
    exchange = {a: b, b: a}
    out = []
    for u, v in readers:
        before = words[loads[u], loads[v]]
        after = words[loads[exchange.get(u, u)], loads[exchange.get(v, v)]]
        if transport and (u, v) == edge:
            after = -after
        out.append(int(after - before))
    return out


def flags(deltas, metric):
    """Binary questions per address, for the equality, unary and squared readers."""
    for delta in deltas:
        if metric == "word":
            yield delta != 0
        elif metric in ("unary", "square"):
            crossings = [step < abs(delta) for step in range(2)]
            if abs(delta) > len(crossings):
                raise ValueError("Unary comparison requires a transported central word")
            if metric == "unary":
                yield from crossings
            else:
                yield from (i and j for i in crossings for j in crossings)
        else:
            raise ValueError("Unknown change reader")


def local_flags(board, edge, words, metric, transport=True):
    readers = addresses(board, edge)
    a, b = edge
    result = {}
    for d in {d for _, d in law.LABELS}:
        loads = np.zeros(board.n, dtype=int)
        loads[a], loads[b] = max(d, 0), max(-d, 0)
        result[d] = list(flags(changes(readers, edge, loads, words, transport), metric))
    return result


def incidence_check(board, words):
    """Every ternary input to each overlapping read, including shared triangle tips."""
    checks = 0
    for edge in board.seams:
        a, b = edge
        for f in addresses(board, edge):
            ports = sorted(set(edge) | set(f))
            for values in product(range(3), repeat=len(ports)):
                state = dict(zip(ports, values, strict=True))
                delta, = changes([f], edge, state, words)
                expected = 0 if f == edge else abs(state[a]-state[b])
                if abs(delta) != expected:
                    raise ArithmeticError("Witness depends on the other endpoint load")
                checks += 1
    return checks


def collar_reflection():
    eye = np.eye(1 + law.N_LAB)
    permutation = [0, *(1 + law.LABELS.index((c, -d)) for c, d in law.LABELS)]
    return eye[:, permutation]


def compiled_gate(by_word, t):
    """Run the two-state scratch circuit and uncompute every flag on every word."""
    identity = np.eye(2)
    flip = identity[::-1]
    hadamard = np.array([[1, 1], [1, -1]])/np.sqrt(2.)
    rail = (identity+flip)/2 + np.exp(-1j*t)*(identity-flip)/2
    phases = []
    for bit in (False, True):
        prepare = hadamard @ (flip if bit else identity)
        returned = prepare.T @ rail @ prepare[:, 0]
        if abs(returned[1]) > 1e-14:
            raise ArithmeticError("The phase rail retained a material record")
        phases.append(returned[0])
    diagonal = [1, *(np.prod([phases[int(bit)] for bit in by_word[d]]) for _, d in law.LABELS)]
    refl = collar_reflection()
    eye = np.eye(len(refl))
    signed_pulse = (eye+refl)/2 + np.exp(-1j*t)*(eye-refl)/2
    echo = law.echo()
    return echo @ np.kron(np.diag(diagonal) @ signed_pulse, np.eye(law.N_SRC)) @ echo.T


def material_generator(by_word):
    """Generator of the commuting rail product, in the physical echo basis."""
    refl = collar_reflection()
    energies = [0, *(sum(by_word[d]) for _, d in law.LABELS)]
    echo = law.echo()
    return echo @ np.kron((np.eye(len(refl))-refl)/2 + np.diag(energies), np.eye(law.N_SRC)) @ echo.T


def ramsey(gate, d):
    """Empty-port written seam; the balanced even sign combination beats against blank."""
    blank = np.eye(law.DIM)[:, law.SRC.index((0, 0))]
    mode = np.zeros(law.DIM)
    for sign in (d, -d):
        mode[(1+law.LABELS.index((0, sign)))*law.N_SRC] = 1/np.sqrt(2.)
    # Work in the stored-label basis and carry both prepared columns through the echo.
    u = law.echo()
    reference, signal = u @ blank, u @ mode
    evolved = gate @ ((reference+signal)/np.sqrt(2.))
    quadratures = []
    for phase in (1, 1j):
        yes = abs(np.vdot((reference+phase*signal)/np.sqrt(2.), evolved))**2
        quadratures.append([float(yes), float(1-yes)])
    return quadratures


def run():
    board, words = cell(), signed_table()
    edge = board.seams[0]
    result = {"checked_word_cases": incidence_check(board, words), "readers": {}}
    reference_exchange, reference_contrast = law.local_term_parts()  # audited comparison only
    from model.candidate_laws.interval_echo.candidate import KAPPA
    pin = reference_exchange + KAPPA*reference_contrast
    binary = [law.SRC.index(x) for x in ((0, 0), (0, 1), (1, 0), (1, 1))]
    for metric in ("word", "unary", "square"):
        by_word = local_flags(board, edge, words, metric)
        counts = [sum(by_word[d]) for d in (1, 2)]
        for other in board.seams:
            other_flags = local_flags(board, other, words, metric)
            if counts != [sum(other_flags[d]) for d in (1, 2)]:
                raise ArithmeticError("The primitive seam chart selected a strength")
        native = material_generator(by_word)
        gate = compiled_gate(by_word, 1/3)
        np.testing.assert_allclose(gate, expm(-1j*native/3), atol=3e-13, rtol=0)
        np.testing.assert_allclose(native[:, binary], pin[:, binary], atol=1e-13, rtol=0)
        error = float(np.linalg.norm(native-pin, ord=2))
        if (error < 1e-12) != (metric == "square"):
            raise ArithmeticError("Unexpected complete-collar transfer")
        result["readers"][metric] = dict(counts=counts, full_operator_error=error,
                                         phase_read=[ramsey(gate, d) for d in (1, 2)])
    square = local_flags(board, edge, words, "square")
    zero = {d: [False]*len(bits) for d, bits in square.items()}
    exchange = material_generator(zero)
    action = material_generator(square)-exchange
    # Compare complete generators, avoiding wrapped single-exposure phases.
    result["relative_rate_controls"] = {}
    for rate in (0.5, 1., 2.):
        residual = float(np.linalg.norm(exchange + rate*action - pin, ord=2))
        result["relative_rate_controls"][str(rate)] = residual
        if (residual < 1e-12) != (rate == 1.):
            raise ArithmeticError("An independently rescaled record clock escaped the full-collar test")
    counts = {}
    for level in range(3):
        actual = federation(level)
        weights = [len(addresses(actual, e))-1 for e in actual.seams]
        values, multiplicity = np.unique(weights, return_counts=True)
        counts[actual.cells] = {int(k): int(n) for k, n in zip(values, multiplicity, strict=True)}
    result["assembled_word_counts"] = counts
    # For one changed unit, a signed central word read in the old chart changes too.
    result["untransported_word"] = sum(local_flags(board, edge, words, "word", False)[1])
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
