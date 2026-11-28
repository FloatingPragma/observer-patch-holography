"""Unrestricted round-trip ranging and landmark bearings on the native release."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np
from scipy.sparse.csgraph import laplacian
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

from . import closed_passage as passage
from . import collective_geometry as geom
from . import landmark_collective
from . import collective_source as source

HERE = source.tidal.HERE


def grounded_inverse(adj):
    # Ground the last vertex. Check the entire integer inverse on its
    # actual sparse Laplacian before using any of its distance entries.
    L = laplacian(adj)[:-1, :-1]
    matrix = DomainMatrix.from_list(L.tolist(), ZZ).to_sparse()
    inv, den = matrix.inv_den(method="rref")
    assert matrix * inv == matrix.eye(matrix.shape, ZZ) * den
    grounded = np.zeros(adj.shape, dtype=object)
    grounded[:-1, :-1] = inv.to_list()
    return grounded, int(den)


def ranges(base, contacts, chances):
    H, den = base
    for (u, v), p in zip(contacts, chances, strict=True):
        z = H[:, u] - H[:, v]
        scale = p.denominator * den + p.numerator * (z[u] - z[v])
        H = scale * H - p.numerator * np.outer(z, z)
        den *= scale
        common = math.gcd(den, *H.ravel())
        H, den = H // common, den // common
    d = H.diagonal()
    return d[:, None] + d[None, :] - 2 * H, den


def release(R, slots):
    # Twice the blank-star generator acts on the landing position (row).
    raw, den = R
    field = np.zeros_like(raw)
    for neighbor in slots.T:
        field += raw[neighbor] - raw
    return 2 * field, den


def row_vector(tab, i):
    return [Q(x, tab[1]) for x in tab[0][i]]


def difference(A, B):
    a, d = A
    b, e = B
    return a * e - b * d, d * e


def endpoint_counts(slots, starts, ends, seed):
    """Count first hits using blocks of twelve exact native slot selections."""
    n = len(slots)
    a, b = ends
    closed = slots.copy()
    closed[[a, b]] = np.array([[a], [b]])
    ticks = 12
    block = np.eye(n, dtype=np.int64)
    for _ in range(ticks):
        block = block[closed].sum(1)
    tickets = 6**ticks
    assert np.all(block.sum(1) == tickets)
    # Offsetting consecutive rows permits one vectorized integer search.
    cdf = (block.cumsum(1) + np.arange(n)[:, None] * tickets).ravel()
    rng = np.random.default_rng(seed)
    trials = 2**17
    rows = {}
    for start in [*starts, "excursion"]:
        live = (
            slots[a, rng.integers(6, size=trials)]
            if start == "excursion"
            else np.full(trials, start)
        )
        hit = int((live == b).sum())
        live = live[(live != a) & (live != b)]
        for _ in range(1000):
            if not len(live):
                break
            outcomes = live * tickets + rng.integers(tickets, size=len(live))
            dest = np.searchsorted(cdf, outcomes, side="right") % n
            hit += int((dest == b).sum())
            live = dest[(dest != a) & (dest != b)]
        pending = len(live)
        err = np.sqrt(np.log(2 * (len(starts) + 1) / 0.003) / (2 * trials))
        rows[str(start)] = dict(
            b=hit,
            a=trials - hit - pending,
            pending=pending,
            bounds=[
                max(0.0, hit / trials - err),
                min(1.0, (hit + pending) / trials + err),
            ],
        )
    return dict(
        trials=trials,
        seed=seed,
        ticks_per_block=ticks,
        maximum_blocks=1000,
        family_error=0.003,
        launches=rows,
    )


def precursor(slots, base, contact, measurements, contrast):
    u, v = contact
    H, den = base
    z = H[:, u] - H[:, v]
    R = Q(z[u] - z[v], den)
    h = [Q(z[u] - x, z[u] - z[v]) for x in z]
    stars = {x: sorted(set(slots[x]) - {x}) for x in (25, 104)}
    starts = sorted(set(stars).union(*stars.values()))
    data = endpoint_counts(slots, starts, contact, 144732)
    rows = data["launches"]
    excursion = rows["excursion"]["bounds"]
    rlo, rhi = 1 / (6 * excursion[1]), 1 / (6 * excursion[0])
    assert rlo <= float(R) <= rhi
    for x in starts:
        lo, hi = rows[str(x)]["bounds"]
        assert lo <= float(h[x]) <= hi

    def prices(key):
        return measurements[key]["measured_without_coincidence"]

    plo, phi = prices("original")
    qlo, qhi = prices("source")
    assert phi < qlo

    # R²(q-p)/[(1+pR)(1+qR)] increases in q and R and decreases in p.
    def scale(p, q, r):
        return (q - p) * r * r / ((1 + p * r) * (1 + q * r))

    factors = scale(phi, qlo, rlo), scale(plo, qhi, rhi)
    fields = {}
    for x, neighbors in stars.items():
        lo, hi = rows[str(x)]["bounds"]
        gradients = [rows[str(y)]["bounds"] for y in neighbors]
        lower = sum(max(0.0, lo - b, a - hi) ** 2 for a, b in gradients)
        upper = sum(max(abs(lo - b), abs(hi - a)) ** 2 for a, b in gradients)
        assert np.all(contrast[0][x] == contrast[0][x, x])
        expected = Q(contrast[0][x, x], contrast[1])
        change = [-2 * factors[1] * upper, -2 * factors[0] * lower]
        assert change[0] <= float(expected) <= change[1]
        fields[str(x)] = dict(
            neighbors=list(map(int, neighbors)),
            exact=str(expected),
            measured_interval=change,
        )
    return dict(
        records=data, bypass_resistance=str(R), bypass_bounds=[rlo, rhi], fields=fields
    )


def shell_profile(distance, eligible, repair, single):
    # Use ordinary, unmodified round-trip lengths as the radial enrollment.
    radii = np.array([float(x) for x in distance])
    upper = float(radii[eligible].max())
    width = upper / 10
    labels = np.minimum((radii / width).astype(int), 9)
    shells = []
    for k in sorted(set(labels[eligible])):
        ids = np.flatnonzero(eligible & (labels == k))
        shells.append(
            dict(
                bin=int(k),
                count=len(ids),
                radius=float(np.mean(radii[ids])),
                inward=float(sum(-repair[i] for i in ids) / len(ids)),
                **(
                    {"ordinary": float(sum(-single[i] for i in ids) / len(ids))}
                    if single is not None
                    else {}
                ),
            )
        )
    return dict(bin_width=width, shells=shells)


def panel(board, layout, counts):
    marks = [s for _, seams in layout for s in seams]
    contacts = [seams[-1] for _, seams in layout]
    adj, slots = passage.contacts(board, marks)
    survey = geom.survey(board, layout, False)
    eligible = survey.clearance >= 2
    ids = np.flatnonzero(eligible)
    base = grounded_inverse(adj)
    m = len(contacts)
    old = source.chances(board.n, ["original"] * m)
    new = source.chances(board.n, ["source"] * m)
    assert all(x < y for x, y in zip(old, new, strict=True))
    a, b, c = source.ordinary_coefficients()
    equivalent = a * board.n + b * m + c
    assert source.chances(board.n, ["original"] * m, equivalent) == new
    R0 = ranges(base, contacts, old)
    A0 = release(R0, slots)
    arms = dict(source=new, single=source.chances(board.n, ["original"] * m, m))
    if m == 1:
        arms["matched_charge"] = [Q(counts["low"]["success"])]
    differences = {}
    result = {}
    for name, ps in arms.items():
        R1 = ranges(base, contacts, ps)
        dA = difference(release(R1, slots), A0)
        field, den = dA
        # Every column is a different observer's fixed ranging landmark.
        assert np.all(field[ids] == field[ids, :1])
        # Independently apply the grounded resistance embedding to the two
        # landmark bearing. A free probe has zero initial bearing curvature
        # away from either beacon, measured separately in both arms.
        for raw, divisor in (R0, R1):
            for i, j in ((104, 116), (25, 36)):
                h2 = raw[i] - raw[j]
                f = (h2[slots] - h2[:, None]).sum(1)
                assert np.all(f[ids[(ids != i) & (ids != j)]] == 0)
                for k, sign in ((i, 1), (j, -1)):
                    if eligible[k]:
                        assert f[k] == sign * 2 * divisor
        values = [Q(int(x), den) for x in field[:, 0]]
        inward = [int(i) for i in ids if values[i] < 0]
        outward = [int(i) for i in ids if values[i] > 0]
        result[name] = dict(
            probabilities=list(map(str, ps)),
            inward=len(inward),
            outward=len(outward),
            silent=len(ids) - len(inward) - len(outward),
            maximum_inward_clearance=max(
                (int(survey.clearance[i]) for i in inward), default=0
            ),
            mean_change=str(sum(values[i] for i in ids) / len(ids)),
            independent_landmark_difference=0,
        )
        differences[name] = values
        if m == 1 and name == "source":
            result["held_precursor"] = precursor(slots, base, contacts[0], counts, dA)
            H, denom = base
            u, v = contacts[0]
            z = H[:, u] - H[:, v]
            rb = Q(z[u] - z[v], denom)
            factor = (ps[0] - old[0]) / ((1 + ps[0] * rb) * (1 + old[0] * rb))
            for i in ids:
                energy = Q(sum((z[j] - z[i]) ** 2 for j in slots[i]), denom**2)
                assert values[i] == -2 * factor * energy
    return dict(
        paths=[p for p, _ in layout],
        contacts=contacts,
        stations=board.n,
        eligible=len(ids),
        landmark_launch_comparisons=board.n * len(ids),
        original_probabilities=list(map(str, old)),
        matching_ordinary_singles=str(equivalent),
        arms=result,
        center=survey.center,
        profiles={
            str(anchor): shell_profile(
                row_vector(R0, anchor),
                eligible,
                differences["source"],
                differences["single"],
            )
            for anchor in (survey.center, *contacts[0])
        }
        if m == 1
        else {},
    )


def run():
    source.check_gate()
    # A comparison law repeats the native blank signed-load swap in every
    # head sector. Both orientations transmit either load code on every word.
    swap = source.tidal.lab.gates()[0][:9, :9]
    blind = np.kron(np.eye(11, dtype=int), swap)
    for m in range(11):
        for x in range(3):
            for y in range(3):
                assert blind[9 * m + 3 * y + x, 9 * m + 3 * x + y] == 4
    board = source.tidal.lab.energy.pair.stock.observations.federation(0)
    priors = json.loads((HERE / "exterior_tide.json").read_bytes())
    counts = priors["local_source"]["preparations"]
    results = {
        str(m): panel(board, source.packing(board, "dispersed", m), counts)
        for m in (1, 4)
    }
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "LANDMARK_RELEASE.md",
                HERE / "exterior_tide.json",
                HERE / "collective_relay_lab.json",
            )
        },
        panels=results,
        memory_blind_gate=dict(
            checked_input_columns=99, local_success=1, range_change=0
        ),
        unrestricted_circuit_launches=landmark_collective.eight_source_witness(
            json.loads((HERE / "collective_relay_lab.json").read_bytes())
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    data = json.dumps(result, indent=2, allow_nan=False) + "\n"
    path = HERE / "landmark_release.json"
    if args.write:
        path.write_text(data)
    else:
        assert path.read_text() == data
    for name, row in result["panels"].items():
        print(name, row["eligible"], row["arms"]["source"]["inward"])
