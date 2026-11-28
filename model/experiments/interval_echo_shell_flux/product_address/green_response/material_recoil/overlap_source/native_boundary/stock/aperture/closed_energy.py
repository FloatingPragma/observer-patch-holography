"""Native local and whole-board energies of a closed repair at fixed head marginals."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np
from scipy.optimize import linprog

from . import closed_pairing as pair
from . import exchange_geometry
from .. import parked

HERE = pair.HERE


def returned(w):
    stock = pair.stock
    owners, seams, _, _, *_ = stock.carrier()
    board = stock.observations.federation(0)
    home = stock.observations.MATTER
    flat, total_blank = parked.ordinary(board, owners, home)
    ix = w[:, :4] @ 3 ** np.arange(4)
    flat = flat[ix]
    h = stock.kernels()[2]
    diag = h.diagonal().reshape(11, 9)
    terms = []
    for u, v in board.seams:
        code = 3 * flat[:, u] + flat[:, v]
        mark = w[:, 4 + seams.index((u, v))] if (u, v) in seams else 0
        terms.append(diag[mark, code])
    private = 32 * (flat**2) * (6 - board.degree)
    base = total_blank[ix].copy()
    for j, (u, v) in enumerate(seams):
        code = 3 * flat[:, u] + flat[:, v]
        base += diag[w[:, 4 + j], code] - diag[0, code]
    assert np.array_equal(np.sum(terms, axis=0) + private.sum(1), base)
    return (
        board,
        (12 * np.array(owners) + home).tolist(),
        seams,
        np.array(terms),
        private,
        base,
    )


def compressed_hamiltonian(w, board, home, seams, diagonal):
    """Four times the returned Hamiltonian restricted to the supplied sharp records."""
    h = pair.stock.kernels()[2]
    lookup = {tuple(word): i for i, word in enumerate(w)}
    locations = {x: j for j, x in enumerate(home)}
    marks = {e: j for j, e in enumerate(seams)}
    matrix = np.diag(diagonal.copy())
    for u, v in board.seams:
        j, k = locations.get(u), locations.get(v)
        s = marks.get((u, v))
        if j is None and k is None and s is None:
            assert not h[:, 0].any()
            continue
        for i, word in enumerate(w):
            left = word[j] if j is not None else 0
            right = word[k] if k is not None else 0
            mark = word[4 + s] if s is not None else 0
            old = 9 * mark + 3 * left + right
            for new in np.flatnonzero(h[:, old]):
                if new == old:
                    continue
                label, ports = divmod(new, 9)
                v1, v2 = divmod(ports, 3)
                if (j is None and v1) or (k is None and v2) or (s is None and label):
                    continue
                out = word.copy()
                if j is not None:
                    out[j] = v1
                if k is not None:
                    out[k] = v2
                if s is not None:
                    out[4 + s] = label
                dest = lookup.get(tuple(out))
                if dest is not None:
                    matrix[dest, i] += h[new, old]
    assert np.array_equal(matrix, matrix.T)
    return matrix


def all_supported_prices(w, cells, board, home, seams, diagonal):
    """Bound total energy for every preparation on the accepted returned support."""
    matrix = compressed_hamiltonian(w, board, home, seams, diagonal)
    x, y = np.nonzero(matrix)
    assert np.array_equal(cells[x], cells[y])
    radius = np.abs(matrix).sum(1) - np.abs(diagonal)
    t = np.bincount(cells, minlength=4).reshape(2, 2)
    goals = np.outer(t.sum(1), t.sum(0)).ravel()
    n = len(w)
    avg = Q(int(diagonal.sum()), 4 * n)

    def shift(costs, pick):
        return (
            sum(
                Q(int(g) * int(pick(costs[cells == k])), 4 * n * n)
                for k, g in enumerate(goals)
            )
            - avg
        )

    return dict(
        sharp_min=str(shift(diagonal, min)),
        sharp_max=str(shift(diagonal, max)),
        coherent_lower=str(shift(diagonal - radius, min)),
        coherent_upper=str(shift(diagonal + radius, max)),
        maximum_offdiagonal_row_sum_over_four=str(Q(int(radius.max()), 4)),
    )


def least_change(cells):
    t = np.bincount(cells, minlength=4).reshape(2, 2)
    p = np.outer(t.sum(1), t.sum(0))
    return np.array(
        [Q(int(p.ravel()[k]), int(t.sum() * t.ravel()[k])) for k in cells], dtype=object
    )


def extreme_total(w, cells, codes, h, ids, groups, total):
    base = pair.constraints(w, cells, codes, h, "energy")
    A = pair.group_matrix(np.vstack((base, w[:, :4].T, cells == 3)), ids, groups)
    b = A.sum(1)
    b[-1] = 144
    e = pair.group_matrix(total[None, :], ids, groups)[0]
    result = {}
    for sign, tag in ((-1, "maximum"), (1, "minimum")):
        opt = linprog(sign * e, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
        assert opt.success
        z = np.array(pair.rational(opt.x), object)
        dual = np.array(pair.rational(opt.eqlin.marginals), object)
        assert all(z >= 0) and np.all(A @ z == b)
        assert np.all(dual @ A <= sign * e) and dual @ b == sign * (e @ z)
        result[tag] = dict(
            energy_increment=str((e @ z - int(total.sum())) / (4 * len(w))),
            weights=[str(x) for x in z],
            dual=[str(x) for x in dual],
        )
    return result


def spectrum(draws, dim, rng):
    # A sharp n against a blank zero neighbor gives exact levels 0,8,9,32,33.
    # The two nontrivial sharp inputs split equally between exchange parities.
    dist = np.bincount(dim, weights=np.asarray(draws, float), minlength=3) / sum(draws)
    p = np.array([dist[0], dist[1] / 2, dist[1] / 2, dist[2] / 2, dist[2] / 2], float)
    n = 2**32
    hist = rng.multinomial(n, p / p.sum())
    return hist.tolist(), float(hist @ np.array([0, 8, 9, 32, 33]) / n)


def run():
    _, w, cells, _, codes, h, ids, groups = pair.setup()
    board, home, seams, terms, private, total = returned(w)
    n = len(w)
    z = least_change(cells)
    shift = (z - 1) / n
    assert all(shift @ w[:, :4] == 0)
    # For n=0,1,2 the extra sharp energy on the unmarked, empty-neighbor
    # contact equals (31/4)(n-1)^2+(65/4)n-31/4, in native clock units.
    middle = w[:, 1]
    assert np.array_equal(middle != 1, cells == 0)
    ordinary_edges = [
        i for i, e in enumerate(board.seams) if home[1] in e and e not in seams
    ]
    for i in ordinary_edges:
        np.testing.assert_array_equal(
            terms[i], 31 * (middle - 1) ** 2 + 65 * middle - 31
        )
    inc = shift @ total / 4
    edge = shift @ terms[ordinary_edges[0]] / 4
    assert edge == Q(31, 56)
    # Native control on the same read setting with contrast stiffness removed.
    exch = pair.stock.kernels()[3][:9, :9]
    assert np.array_equal(exch.diagonal()[3 * middle], 2 * (middle != 0))
    exchange_price = shift @ exch.diagonal()[3 * middle] / 4
    rng = np.random.default_rng(5129367)
    shots = 2**32
    orig, lo = spectrum(np.ones(n), middle, rng)
    repaired, hi = spectrum(z, middle, rng)
    err = 2 * 33 * np.sqrt(np.log(400) / (2 * shots))
    measured = [hi - lo - err, hi - lo + err]
    assert measured[0] < float(edge) < measured[1]
    head_counts = np.bincount(cells, minlength=4).reshape(2, 2)
    p = head_counts / n
    q = np.outer(p.sum(1), p.sum(0))
    # KL chain rule: keeping the native conditional law given the two heads
    # attains the least possible loss of entropy from the uniform reference.
    entropy_cost = float(np.sum(q * np.log(q / p)))
    # On the blank link park any pre-existing target word in its blank side rail,
    # send the double, read its arrival, return it, and restore the target.
    swap = exchange_geometry.contact_swap()
    slots = []
    for load in (0, 1, 2):
        state = [2, load, 0]
        pulses = 0

        def pulse(i, j):
            nonlocal pulses
            word = 3 * state[i] + state[j]
            dest = int(np.argmax(swap[:, word]))
            assert abs(swap[dest, word] - 1) < 1e-12
            state[i], state[j] = divmod(dest, 3)
            pulses += 1

        if load:
            pulse(1, 2)
        pulse(0, 1)
        assert state[1] == 2
        pulse(0, 1)
        if load:
            pulse(1, 2)
        assert state == [2, load, 0]
        slots.append(pulses // 2)
    assert slots == [1, 2, 2]
    timing = []
    for weights in (np.ones(n), z):
        p = float(sum(weights[middle != 0]) / sum(weights))
        timing.append(int(rng.binomial(shots, p)))
    dt = (timing[1] - timing[0]) / shots
    tw = 2 * np.sqrt(np.log(400) / (2 * shots))
    assert dt - tw < -1 / 28 < dt + tw
    time = {}
    matched = json.loads((HERE / "closed_pairing.json").read_bytes())["menus"]["energy"]
    matched_weights = np.array([Q(x) for x in matched["weights"]])[ids]
    for tags, zz in (("simple", z), ("energy", matched_weights)):
        d = (zz - 1) / n
        time[tags] = dict(
            total_increment=str(d @ total / 4),
            home_means=[str(x) for x in d @ w[:, :4]],
            per_seam_nonzero={
                f"{u},{v}": str(x)
                for (u, v), x in zip(board.seams, (terms @ d) / 4, strict=True)
                if x
            },
            private_nonzero={
                str(i): str(x) for i, x in enumerate(d @ private / 4) if x
            },
        )
    source_bank = json.loads((HERE / "pair_contact.json").read_bytes())
    from .pair_probe import log_factor

    source_band = np.asarray(log_factor(source_bank["source"]), float)
    energy_band = (31 / 4) * source_band
    echo_band = -source_band[::-1] / 2
    assert max(measured[0], energy_band[0]) <= min(measured[1], energy_band[1])
    assert max(dt - tw, echo_band[0]) <= min(dt + tw, echo_band[1])
    return dict(
        forecast_from_source=dict(
            displacement=source_band.tolist(),
            bond_energy=energy_band.tolist(),
            half_echo=echo_band.tolist(),
        ),
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "CLOSED_PAIRING.md",
                HERE / "closed_pairing.json",
                HERE / "pair_contact.json",
                HERE.parent / "parked.json",
            )
        },
        ordinary_home_ports=home,
        center=home[1],
        ordinary_edges=[list(board.seams[i]) for i in ordinary_edges],
        protected_echo=dict(
            sender=min(set(board.seams[ordinary_edges[0]]) - {home[1]}),
            receiver=home[1],
            half_roundtrip_by_load=slots,
            slow_counts=timing,
            trials_per_arm=shots,
            measured_half_roundtrip_change=[dt - tw, dt + tw],
            predicted_half_roundtrip_change="-1/28",
        ),
        trials_per_spectrum=shots,
        native_counts=orig,
        repaired_counts=repaired,
        levels=[0, 8, 9, 32, 33],
        measured_bond_difference=measured,
        mandatory_bond_increment=str(edge),
        exchange_only_readout_difference=str(exchange_price),
        canonical_total_increment=str(inc),
        all_supported_prices=all_supported_prices(w, cells, board, home, seams, total),
        native_total_energy=str(Q(int(total.sum()), 4 * n)),
        minimum_whole_word_entropy_loss=entropy_cost,
        preparations=time,
        matched_energy_extremes=extreme_total(w, cells, codes, h, ids, groups, total),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    r = run()
    path = HERE / "closed_energy.json"
    if opt.write:
        path.write_text(json.dumps(r, indent=2, allow_nan=False) + "\n")
    else:
        assert r == json.loads(path.read_bytes())
    print(
        r["mandatory_bond_increment"],
        r["measured_bond_difference"],
        r["canonical_total_increment"],
        {k: v["energy_increment"] for k, v in r["matched_energy_extremes"].items()},
    )
