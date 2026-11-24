"""Energy selection and a polarization clock on a degree-preserving native ring."""

import argparse
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier, rotations


def geometry():
    board = carrier.federation(0)
    cell_edges = set(carrier.cell().seams)
    degrees = board.degree[:12]
    g = min(
        tuple(t)
        for t in rotations()
        if all(
            t[i] != i and t[t[t[i]]] == i and degrees[i] == degrees[t[i]]
            for i in range(12)
        )
    )
    orbits = sorted({tuple(sorted((i, g[i], g[g[i]]))) for i in range(12)})
    rings = {}
    for order in permutations(orbits):
        for start in product(*order):
            if start[0] != 0:
                continue
            seq = (
                tuple(start)
                + tuple(g[i] for i in start)
                + tuple(g[g[i]] for i in start)
            )

            def edges(d):
                return {
                    tuple(sorted((p, seq[(i + d) % 12]))) for i, p in enumerate(seq)
                }

            if edges(1) <= cell_edges:
                chords = edges(2) & cell_edges
                rings[seq] = (edges(1), chords)
    size = max(len(e[1]) for e in rings.values())
    best = min(seq for seq, e in rings.items() if len(e[1]) == size)
    nearest, chords = rings[best]
    full_square = nearest | {
        tuple(sorted((p, best[(i + 2) % 12]))) for i, p in enumerate(best)
    }
    return (
        board,
        best,
        dict(
            rotation=g,
            directed_cycles=len(rings),
            maximal_cycles=sum(len(e[1]) == size for e in rings.values()),
            chords=size,
        ),
        {
            "native_cycle": nearest,
            "native_with_chords": nearest | chords,
            "full_cell": cell_edges,
            "completed_square": full_square,
        },
    )


def state_reads(columns, twist, crystal, orders):
    probabilities = np.sum(abs(columns) ** 2, axis=1)
    assert abs(probabilities.sum() - 1) < 1e-9
    z = probabilities @ twist
    return dict(
        visibility=abs(z),
        phase_turns=(
            None
            if abs(z) < 1e-10
            else 0.5
            if abs(z.imag) < 1e-12 and z.real < 0
            else np.angle(z) / (2 * np.pi)
        ),
        ramsey_bright_xy=[(1 + z.real) / 2, (1 + z.imag) / 2],
        crystalline=probabilities[crystal].sum(),
        order_squared=probabilities @ abs(orders) ** 2,
    )


def panel(board, ring, name, edges, k):
    test = carrier.Carrier(name, 12, edges)
    sec = BinarySector(test, 4)
    occupation = np.array([[int(p in w) for p in ring] for w in sec.states])
    coord = np.arange(12)
    twist = np.exp(2j * np.pi * (occupation @ coord) / 12)
    order = occupation @ np.exp(2j * np.pi * coord / 3) / 4
    crystal = abs(order) > 1 - 1e-10
    assert crystal.sum() == 3
    ref = (
        np.full(12, 6)
        if name in ("completed_square", "native_chords_uniform")
        else test.degree + 6 - board.degree[:12]
    )
    h = sec.hamiltonian(k) + np.diag(k * (occupation @ ref[list(ring)]))
    e, v = np.linalg.eigh(h)
    band = np.flatnonzero(e >= e[-3] - 1e-8)
    groups = np.split(np.arange(len(e)), np.flatnonzero(np.diff(e) > 1e-8) + 1)

    def max_crystal(group):
        loc = v[crystal][:, group]
        return float(np.linalg.eigvalsh(loc @ loc.T)[-1])

    best = max(groups, key=max_crystal)
    # An unbiased four-load ensemble accepted by a top-energy filter.
    initial = v[:, band] / np.sqrt(len(band))
    raw = state_reads(initial, twist, crystal, order)
    readings = []
    for site in (0, 4, 8):
        w = occupation[:, site, None] * initial
        yes = np.linalg.norm(w) ** 2
        modes = v.T @ (w / np.sqrt(yes))
        free = []
        for t in (0, 1, 4, 16, 64):
            cols = v @ (np.exp(-1j * t * e)[:, None] * modes)
            free.append(dict(time=t, **state_reads(cols, twist, crystal, order)))
        readings.append(dict(site=ring[site], acceptance=yes, free=free))
    # Sharp translated preparations give the same phase ratio at any density
    # profile with this total. They have different period-three detector records.
    controls = {}
    for tag, start in [("crystal", (0, 3, 6, 9)), ("compact", (0, 1, 2, 3))]:
        ids = [
            sec.index[tuple(sorted(ring[(i + 4 * j) % 12] for i in start))]
            for j in range(3)
        ]
        modes = v[ids].T
        release = []
        for t in (0, 1, 4, 16, 64):
            cols = v @ (np.exp(-1j * t * e)[:, None] * modes)
            release.append(
                dict(
                    time=t,
                    rotations=[
                        state_reads(col[:, None], twist, crystal, order)
                        for col in cols.T
                    ],
                )
            )
        controls[tag] = release
    return dict(
        edges=len(edges),
        maximum_stationary_crystalline=max_crystal(best),
        maximizing_energy=e[best[0]],
        maximizing_multiplicity=len(best),
        energy_band_dimension=len(band),
        selection_rate=len(band) / len(e),
        gap=e[band[0]] - e[band[0] - 1],
        bandwidth=e[-1] - e[band[0]],
        unconditioned=raw,
        site_read=readings,
        translated_controls=controls,
    )


def run():
    board, ring, meta, graphs = geometry()
    graphs["native_chords_uniform"] = graphs["native_with_chords"]
    graphs["completed_square_native_reference"] = graphs["completed_square"]
    return currents._portable(
        dict(
            ring=ring,
            **meta,
            unit_four_site_twist_turns=4 / 12,
            panels={
                name: {
                    str(k): panel(board, ring, name, edges, k)
                    for k in (int(candidate.KAPPA), 0)
                }
                for name, edges in graphs.items()
            },
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = run()
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
