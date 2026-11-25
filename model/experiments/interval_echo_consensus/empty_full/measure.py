"""Read the empty/full seam multiplier with two conserved loads on four native ports."""

import argparse
import json
from collections import Counter
from itertools import combinations_with_replacement, product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier

LAW = candidate.native_law
HERE = Path(__file__).resolve().parent
PORTS = (0, 1, 5, 7)
TARGET = (0, 1)
RAILS = ((0, 2), (1, 3))
WORDS = [w for w in product(range(3), repeat=4) if sum(w) == 2]
INITIAL = WORDS.index((0, 0, 1, 1))
FULL = WORDS.index((1, 1, 0, 0))
PHASES = (0, np.pi / 2, np.pi, np.pi / 3, -np.pi / 3, 3 * np.pi / 4)
SHOTS = 2**20
MODES = ("bypass", "native", "total", "which", "only_empty", "only_full")
GRID_BITS = 40
GRID = 1 << GRID_BITS


def rate_radius(outcome_bins):
    """Cover every terminal port and rejection in the enrolled histograms."""
    family_size = outcome_bins * len(MODES) * len(PHASES)
    return float(np.sqrt(np.log(2 * family_size / 0.01) / (2 * SHOTS)))


ERROR = rate_radius(len(WORDS) + 1)


def draw(shots, probabilities, rng, budget):
    """Count starts drawn from integer weights on the 2**-40 grid.

    Every outcome except the most probable has weight round(p * 2**40); the
    most probable takes the complement. Each start compares one PCG64 integer
    below 2**40 with the cumulative weights, so no floating-point operation
    enters a draw and a seed gives the same counts on every platform. Rounded
    outcomes move by at most 2**-41, the complement by the other roundings
    plus the caller's norm budget.
    """
    p = np.asarray(probabilities, dtype=float).ravel()
    weights = np.rint(p * GRID).astype(np.int64)
    top = int(np.argmax(weights))
    weights[top] = GRID - (weights.sum() - weights[top])
    bound = budget + (p.size - 1) / (2 * GRID)
    if weights.min() < 0 or np.max(abs(weights / GRID - p)) > bound:
        raise ValueError("Sampling weights leave the propagated distribution")
    starts = rng.integers(0, GRID, size=shots, dtype=np.int64)
    edges = np.cumsum(weights)
    return np.bincount(np.searchsorted(edges, starts, side="right"), minlength=p.size)


def embed(block, edge):
    """Apply a number-conserving endpoint operation with spectator ports unchanged."""
    full = np.zeros((len(WORDS), len(WORDS)), dtype=complex)
    for i, word in enumerate(WORDS):
        key = tuple(word[j] for j in edge)
        for j, output in enumerate(LAW.SRC):
            if sum(output) != sum(key):
                continue
            new = list(word)
            new[edge[0]], new[edge[1]] = output
            full[WORDS.index(tuple(new)), i] = block[j, LAW.SRC.index(key)]
    return full


def phase(port, angle):
    return np.diag([np.exp(1j * angle * w[port]) for w in WORDS])


def preparation():
    for a, b in (TARGET, *RAILS):
        if tuple(sorted((PORTS[a], PORTS[b]))) not in carrier.CELL_SEAMS:
            raise ValueError("Pulse on a non-native edge")
    # Blank memories reduce the registered Hamiltonian exactly to this load block.
    exchange, contrast = LAW.local_term_parts()
    h = (exchange + candidate.KAPPA * contrast)[: LAW.N_SRC, : LAW.N_SRC]

    def pulse(edge, duration):
        return embed(expm(-1j * duration * h), edge)

    steps = []
    for rail in RAILS:
        steps.extend([pulse(rail, np.pi / 2), phase(rail[0], -np.pi / 2)])
    half = np.pi / (4 * (candidate.KAPPA + 0.5))
    for _ in range(2):
        steps.extend([pulse(TARGET, half), phase(0, np.pi)])
    steps.append(pulse(RAILS[1], np.pi / 2))
    result = np.eye(len(WORDS))
    for u in steps:
        result = u @ result
    return result


def total_sensitive():
    """Replace the full-seam record by an orthogonal exchange-even record."""
    write = LAW.write()
    out = LAW.SRC.index((1, 1))
    write[:, out] = 0
    # Even under simultaneous coin and contrast reversal; orthogonal to both
    # symmetric coin vectors written by (2,0) and (0,2).
    for coin in (0, 1):
        for d in (-2, 2):
            i = LAW.LABELS.index((coin, d))
            write[i * LAW.N_SRC + out, out] = (1 if coin == 0 else -1) * d / 4
    np.testing.assert_allclose(write.T @ write, np.eye(LAW.N_SRC), atol=3e-15)
    return np.block(
        [
            [np.zeros((LAW.N_SRC, LAW.N_SRC)), write.T],
            [write, np.eye(len(write)) - write @ write.T],
        ]
    )


def kraus(unitary):
    blocks = unitary[:, : LAW.N_SRC].reshape(1 + LAW.N_LAB, LAW.N_SRC, LAW.N_SRC)
    return [embed(b, TARGET) for b in blocks if np.any(b)]


def record_branches(state, mode, instruments):
    if mode == "bypass":
        return [state], 0.0
    occupancy = np.array([w[0] + w[1] for w in WORDS])
    if mode in ("native", "total"):
        return [k @ state for k in instruments[mode]], 0.0
    admissions = {"which": (0, 1, 2), "only_empty": (0,), "only_full": (2,)}[mode]
    passed = [state * (occupancy == n) for n in admissions]
    rejected = float(
        sum(abs(state[i]) ** 2 for i, n in enumerate(occupancy) if n not in admissions)
    )
    return [k @ s for s in passed for k in instruments["native"]], rejected


def experiment():
    prep = preparation()
    state = prep[:, INITIAL]
    instruments = {
        name: kraus(u)
        for name, u in (("native", LAW.echo()), ("total", total_sensitive()))
    }
    rng = np.random.default_rng(914730)
    output = {}
    for mode in MODES:
        accepted, rejected = record_branches(state, mode, instruments)
        banks = []
        for angle in PHASES:
            reader = prep.conj().T @ phase(0, angle)
            p = np.r_[sum(abs(reader @ s) ** 2 for s in accepted), rejected]
            np.testing.assert_allclose(p.sum(), 1.0, atol=4e-14)
            counts = draw(SHOTS, p, rng, 4e-14)
            banks.append(dict(phase=float(angle), counts=counts.tolist()))
        output[mode] = banks
    cal, held = decode(output)
    return dict(
        ports=PORTS,
        words=WORDS,
        shots=SHOTS,
        sampling_grid_bits=GRID_BITS,
        simultaneous_rate_error=ERROR,
        prepared_probabilities=np.round(abs(state) ** 2, 13).tolist(),
        banks=output,
        calibration=cal,
        held=held,
    )


def decode(banks):
    calibration, held = {}, {}
    for mode, rows in banks.items():
        rates = np.array([r["counts"][INITIAL] / SHOTS for r in rows])
        baseline = 0.25 if mode.startswith("only_") else 0.5
        x, y = rates[:2] - baseline
        calibration[mode] = {"cosine": float(x), "sine": float(y)}
        held[mode] = []
        for angle, observed in zip(PHASES[2:], rates[2:], strict=True):
            c, s = np.cos(angle), np.sin(angle)
            expected = baseline + c * x + s * y
            budget = ERROR * (1 + abs(c) + abs(s))
            held[mode].append(
                dict(
                    forecast=float(expected),
                    observed=float(observed),
                    joint_error=float(budget),
                )
            )
    return calibration, held


def exact_channel():
    """Count equal record labels on every matrix unit using integer amplitudes."""
    cell_words = list(combinations_with_replacement(range(12), 2))
    keys = [(w.count(0), w.count(1)) for w in cell_words]

    def channel(a, b, total):
        out = Counter()
        for c in (0, 1):
            ya, yb = LAW.native_event(*a, c), LAW.native_event(*b, c)
            ja, jb = a[0] - ya[0], b[0] - yb[0]
            la = (ja, sum(a) if total else sum(a) % 2)
            lb = (jb, sum(b) if total else sum(b) % 2)
            if la == lb:
                out[ya, yb] += 1
        # Multiplicities in half units, with the unchanged spectators implicit.
        return out

    def blocks(u):
        # Scale the amplitudes to 2 or sqrt(2) times integers. Deterministic
        # equal-coin and odd-coin columns use sqrt(2); the special foil uses 2.
        arr = u[:, : LAW.N_SRC].reshape(-1, LAW.N_SRC, LAW.N_SRC)
        ans = {}
        for a in set(keys):
            A = arr[:, :, LAW.SRC.index(a)]
            factor = (
                2 if np.isclose(np.linalg.norm(A, ord=np.inf), 0.5) else np.sqrt(2.0)
            )
            v = np.rint(factor * A).astype(int)
            np.testing.assert_allclose(factor * A, v, atol=3e-15)
            ans[a] = (v, factor == 2)
        return ans

    native, total = blocks(LAW.echo()), blocks(total_sensitive())

    def match(data, a, b, expected):
        va, sa = data[a]
        vb, sb = data[b]
        lhs = va.T @ vb
        if sa != sb:
            return not expected and not lhs.any()
        rhs = np.zeros_like(lhs)
        for (ya, yb), count in expected.items():
            rhs[LAW.SRC.index(ya), LAW.SRC.index(yb)] = count * (2 if sa else 1)
        return np.array_equal(lhs, rhs)

    e_ok = t_ok = different = 0
    for a, b in product(keys, repeat=2):
        e, t = channel(a, b, False), channel(a, b, True)
        e_ok += match(native, a, b, e)
        t_ok += match(total, a, b, t)
        different += e != t
    # For every local key of any total the parity-record and contrast-record
    # coincidences are equivalent, with c retained and endpoint current j.
    extended = all(
        (a[0] - a[1] == b[0] - b[1])
        == (
            sum(a) % 2 == sum(b) % 2
            and a[0] - LAW.native_event(*a, c)[0] == b[0] - LAW.native_event(*b, c)[0]
        )
        for a, b, c in product(LAW.SRC, LAW.SRC, (0, 1))
    )
    return dict(
        cell_words=len(keys),
        matrix_units=len(keys) ** 2,
        native_matches_E=int(e_ok),
        changed_matches_T=int(t_ok),
        distinguishing_units=different,
        native_label_identity_all_qutrit_keys=extended,
    )


def run():
    obs = experiment()
    obs["exact_channels"] = exact_channel()
    single = LAW.SRC.index((1, 0))
    twice = LAW.echo() @ LAW.echo()[:, single]
    once = LAW.echo()[:, single].reshape(-1, LAW.N_SRC)
    fresh = LAW.echo()[:, : LAW.N_SRC].reshape(-1, LAW.N_SRC, LAW.N_SRC)
    fresh_return = sum(abs((k @ arm)[single]) ** 2 for k in fresh for arm in once)
    obs["two_writes_single_load"] = dict(
        same_memory_return=float(
            np.round(sum(abs(twice.reshape(-1, LAW.N_SRC)[:, single]) ** 2), 13)
        ),
        fresh_memory_return=float(np.round(fresh_return, 13)),
    )
    return json.loads(json.dumps(obs))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = run()
    if args.output:
        args.output.write_text(json.dumps(results, indent=2) + "\n")
    else:
        if results != json.loads((HERE / "results.json").read_text()):
            raise ValueError("Acquisition replay differs from stored outcomes")
    print(results["exact_channels"])
    print(
        {
            m: [r["counts"][INITIAL] for r in rows]
            for m, rows in results["banks"].items()
        }
    )
