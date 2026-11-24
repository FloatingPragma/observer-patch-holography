"""Sum over every physical belt orientation when the walls are released on the intact cell."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier, rotations
from . import belt, multiplets


def colour_certificate(sec, h, actions):
    """Bound repeated copies of a spatial irrep at an energy by exact modular rank."""
    # Sum of the five irreducible A5 characters at orders 1, 2, 3, 5.
    character_sum = {1: 16, 2: 0, 3: 1, 5: 1}
    total = 0
    identity = np.arange(len(sec.states))
    for perm in actions:
        g = identity
        for order in range(1, 6):
            g = g[perm]
            if np.array_equal(g, identity):
                break
        total += character_sum[order] * np.count_nonzero(perm == identity)
    copies, rem = divmod(int(total), 60)
    assert rem == 0

    integer_h = np.rint(2 * h).astype(np.int64)
    np.testing.assert_array_equal(2 * h, integer_h)
    operator = sparse.csr_matrix(integer_h)
    prime = 1009
    vector = np.random.default_rng(12013).integers(
        1, prime, len(sec.states), dtype=np.int64
    )
    pivots = []
    for _ in range(copies):
        column = vector.copy()
        for i, row in pivots:
            column = (column - column[i] * row) % prime
        nonzero = np.flatnonzero(column)
        if not len(nonzero):
            break
        i = int(nonzero[0])
        row = (column * pow(int(column[i]), -1, prime)) % prime
        pivots.append((i, row))
        vector = (operator @ vector) % prime
    distinct_lower = len(pivots)
    return dict(
        spatial_irrep_copies=copies,
        integer_operator_scale=2,
        prime=prime,
        source_seed=12013,
        krylov_rank=distinct_lower,
        maximum_same_energy_spatial_multiplicity=1 + copies - distinct_lower,
    )


def run():
    axes = belt.belts()
    permutation = np.array(multiplets.run()["fifth_turn_permutation"])
    result = {}
    for n in (3, 4):
        guide = BinarySector(belt.ring(10), n)
        sec = BinarySector(carrier.cell(), n)
        sharp = [guide.states[i] for i in belt.selected(guide)]
        charts = defaultdict(list)
        for axis, order in axes.items():
            for s in sharp:
                a, g, _ = belt.wall(s, 10)
                idx = sec.index[tuple(sorted(order[j] for j in s))]
                charts[idx].append(
                    dict(axis=axis, core_ends=sorted((order[a], order[(a + g) % 10])))
                )
        chosen = sorted(charts)
        h = sec.hamiltonian(candidate.KAPPA)
        e, v = np.linalg.eigh(h)
        size = len(chosen)
        overlaps = np.linalg.svd(v[chosen, -size:], compute_uv=False)
        rotation = np.array(
            [sec.index[tuple(sorted(permutation[j] for j in s))] for s in sec.states]
        )
        np.testing.assert_allclose(h[np.ix_(rotation, rotation)], h)
        order = next(iter(axes.values()))
        (source,) = [s for s in sharp if belt.wall(s, 10)[0] == 0]
        (dest,) = [
            s
            for s in sharp
            if belt.signed_distance(belt.wall(s, 10)[2], belt.wall(source, 10)[2], 10)
            == 3
        ]
        indices = [
            sec.index[tuple(sorted(order[j] for j in s))] for s in (source, dest)
        ]
        psi = np.zeros(len(e), complex)
        psi[indices] = np.array([1, 1j]) / np.sqrt(2)
        modes = v.T @ psi
        reads = []
        for t in belt.TIMES:
            p = abs(v @ (np.exp(-1j * t * e) * modes)) ** 2
            reads.append(dict(time=t, any_chart=p[chosen].sum()))
        actions = np.asarray(
            [
                [sec.index[tuple(sorted(g[j] for j in s))] for s in sec.states]
                for g in rotations()
            ]
        )
        assert len(actions) == 60
        levels = []
        commutants = []
        for ix in multiplets.groups(e):
            sub = v[:, ix]
            characters = np.einsum("si,gsi->g", sub, sub[actions])
            inner = np.mean(np.abs(characters) ** 2)
            assert abs(inner - round(inner)) < 1e-9
            commutants.append(round(inner))
            if ix[-1] < len(e) - size:
                continue
            action = sub.T @ sub[rotation]
            levels.append(
                dict(
                    energy=e[ix].mean(),
                    dimension=len(ix),
                    spatial_character_norm=inner,
                    source_weight=(abs(modes[ix]) ** 2).sum(),
                    chart_weight=(abs(sub[chosen]) ** 2).sum() / len(ix),
                    fifth_turn_trace=np.trace(action),
                    fifth_turn_identity_distance=np.linalg.norm(
                        action - np.eye(len(ix)), 2
                    ),
                )
            )
        result[str(n)] = dict(
            sharp_words=size,
            exact_colour_certificate=colour_certificate(sec, h, actions),
            all_energy_commutant_dimensions=dict(Counter(str(m) for m in commutants)),
            minimum_occupied_seams=int(sec.agreed.min()),
            exhaust_intact_minimum=chosen == list(map(int, belt.selected(sec))),
            chart_multiplicity=dict(Counter(str(len(r)) for r in charts.values())),
            core_pair_multiplicity=dict(
                Counter(
                    str(len({tuple(c["core_ends"]) for c in r}))
                    for r in charts.values()
                )
            ),
            source_charts=charts[indices[0]],
            gap_below=e[-size] - e[-size - 1],
            every_time_any_chart_lower=max(0.0, 1 - 4 * (1 - overlaps[-1] ** 2)),
            levels=levels,
            reads=reads,
        )
    return currents._portable(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(encoded)
