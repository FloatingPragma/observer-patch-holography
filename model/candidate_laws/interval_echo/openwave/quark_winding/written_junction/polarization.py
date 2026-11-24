"""Phase of a calibrated two-rail probe at two ends of a translated density pattern."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier
from ...electron_field.flux_string.search import Patch

X = np.arange(-96, 217)
WIDTHS = (0, 1, 2, 4, 8, 16, 32)


def layout():
    board = carrier.federation(2)
    graph = nx.Graph(board.seams)
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    gateways = {
        (u // 12, v // 12): u % 12
        for a, b in board.external
        for u, v in ((a, b), (b, a))
    }
    remaining = owners.copy()
    chain = [0]

    def extend():
        if len(chain) * 12 >= 2 * len(X) + 4:
            return True
        a = chain[-1]
        for b in sorted(owners[a]):
            if b in chain:
                continue
            remaining.remove_edge(a, b)
            if nx.is_connected(remaining):
                chain.append(b)
                if extend():
                    return True
                chain.pop()
            remaining.add_edge(a, b)
        return False

    assert extend()
    inside = nx.Graph(carrier.cell().seams)
    path = []
    for j, owner in enumerate(chain):
        a = gateways[owner, chain[j - 1]] if j else None
        b = gateways[owner, chain[j + 1]] if j + 1 < len(chain) else None
        if a is None:
            a = next(i for i in range(12) if i != b)
        if b is None:
            b = next(i for i in range(12) if i != a)
        visit = next(p for p in nx.all_simple_paths(inside, a, b) if len(p) == 12)
        path.extend(12 * owner + i for i in visit)
    sample = path[: len(X)]
    reference = path[len(X) : 2 * len(X)]
    station = path[2 * len(X) : 2 * len(X) + 4]
    for sites in (sample, reference):
        assert all(b in graph[a] for a, b in zip(sites, sites[1:]))
        graph.remove_edges_from(zip(sites, sites[1:]))
    # A blank path can route the material away while every measured head stays in place.
    assert nx.is_connected(graph)
    return sample, reference, station


def common_heads(path):
    # A positive unit moves on every third edge to the right for shift +1.
    # The -2 translation uses the complementary edges in the opposite direction.
    allowed = X[:-1] % 3 == 0
    return np.array(
        [int(x if a < b else not x) for a, b, x in zip(path, path[1:], allowed)]
    )


def translated(path, shift):
    host = (X % 3 == 0).astype(int)
    n = host.copy()
    bit = {
        tuple(sorted((a, b))): 1 << i for i, (a, b) in enumerate(zip(path, path[1:]))
    }
    patch = SimpleNamespace(bit=bit)
    steps = []
    for index in np.flatnonzero((host == 1) & (X >= 0) & (X < 120)):
        for j in range(abs(shift)):
            a = index + j * int(np.sign(shift))
            b = a + int(np.sign(shift))
            steps.append((a, b))
    seed = sum(1 << i for i, c in enumerate(common_heads(path)) if c)
    head = seed
    for a, b in steps:
        edge = tuple(sorted((path[a], path[b])))
        at, head = Patch.move(patch, path[a], head, edge, 1)
        assert n[a] == 1 and n[b] == 0 and at == path[b]
        n[a], n[b] = 0, 1
    flux = []
    for a, b in zip(path, path[1:]):
        flag = bit[tuple(sorted((a, b)))]
        flux.append(
            (int(bool(seed & flag)) - int(bool(head & flag))) * (1 if a < b else -1)
        )
    np.testing.assert_array_equal(np.diff(np.r_[0, flux, 0]) * -1, n - host)
    return (
        host,
        n,
        np.array(flux),
        dict(
            steps=len(steps),
            changed_heads=(seed ^ head).bit_count(),
            aligned_heads_hex=hex(seed),
            final_heads_hex=hex(head),
        ),
    )


def potential(centre, width, plateau=8):
    distance = np.maximum(0, np.abs(X - centre) - plateau)
    if width == 0:
        return (distance == 0).astype(float)
    return (1 + np.cos(np.pi * np.minimum(distance / width, 1))) / 2


def probe(sample, reference, profile, theta, kappa, repeats=4):
    """Two independent sharp source bits touch the two coherent rails at each address."""
    dt = theta * profile / (2 * (candidate.KAPPA + 0.5) * repeats)
    # Rows condition on the actual initial probe word 10 or 01. Columns are its two rails.
    means = np.array([[1.0, 0.0], [0.0, 1.0]])
    for _ in range(repeats):
        for i in np.flatnonzero(dt):
            transfer = np.sin(dt[i] / 2) ** 2
            means = (1 - transfer) * means + transfer * np.array(
                [sample[i], reference[i]]
            )
    zero = np.mean(np.prod(1 - means, axis=1))
    double = np.mean(np.prod(means, axis=1))
    visibility = np.prod(np.cos(dt / 2) ** 2) ** repeats
    phase = 2 * (kappa + 0.5) * (dt @ (sample - reference)) * repeats
    coherence = visibility * np.exp(1j * phase)
    single = 1 - zero - double
    out = {}
    for axis, v in [("x", coherence.real), ("y", coherence.imag)]:
        bins = [(single + v) / 2, (single - v) / 2, zero, double]
        assert min(bins) > -2e-12 and abs(sum(bins) - 1) < 1e-12
        out[axis] = bins
    return dict(outcomes=out, phase=float(np.angle(coherence)), visibility=visibility)


def run():
    source, reference, station = layout()
    cal = potential(0, 32)
    empty = np.zeros(len(X), int)
    whole = (X == 0).astype(int)
    cases = {}
    for shift in (1, -2):
        host, n, flux, heads = translated(source, shift)
        rows = []
        for name, centre, plateau in [
            ("left", 0, 8),
            ("right", 120, 8),
            ("outside", -54, 8),
            ("enclose", 60, 70),
        ]:
            for offset in (0, 1):
                for width in WIDTHS:
                    f = potential(centre + offset, width, plateau)
                    q = float(f @ (n - host))
                    field = float(np.diff(f) @ flux)
                    rows.append(
                        dict(
                            region=name,
                            offset=offset,
                            width=width,
                            material=q,
                            heads=field,
                            quarter_turn=probe(n, host, f, np.pi / 2, 8),
                            full_turn=probe(n, host, f, 2 * np.pi, 8),
                        )
                    )
        f = potential(0, 32)

        def clock(a, word):
            idx = (X >= a) & (X < a + 3)
            return [int(word[idx & (X % 3 == k)].sum()) for k in range(3)]

        heads["density_clocks"] = {str(a): clock(a, n) for a in (-6, 6, 102, 126)}
        heads["bulk_clocks"] = {str(a): clock(a, host) for a in (-6, 6, 102, 126)}
        cases[str(shift)] = dict(
            source=heads,
            reads=rows,
            zero_contrast=probe(n, host, f, 2 * np.pi, 0),
            no_write_heads=0,
            unit_share_full_turn=[
                probe(whole, empty, cal, 2 * np.pi, 8),
                probe(empty, empty, cal, 2 * np.pi, 8),
            ],
        )
    # Check both physical charge backgrounds against the full 99-state native contact.
    ex, co = law.local_term_parts()
    errors = []
    for theta in (np.pi / 2, 2 * np.pi):
        dt = theta / (2 * (candidate.KAPPA + 0.5) * 4)
        u = expm(-1j * dt * (ex + candidate.KAPPA * co))
        for base in (0, 1):
            for r in (0, 1):
                for n in (0, 1):
                    col = law.SRC.index((n + base, r + base))
                    swapped = law.SRC.index((r + base, n + base))
                    target = np.zeros(law.DIM, complex)
                    if r == n:
                        target[col] = 1
                    else:
                        phase = np.exp(-1j * dt * (candidate.KAPPA + 0.5))
                        target[col] = phase * np.cos(dt / 2)
                        target[swapped] = 1j * phase * np.sin(dt / 2)
                    errors.append(np.linalg.norm(target - u[:, col]))
    return currents._portable(
        dict(
            board_cells=320,
            sample_path=source,
            bulk_path=reference,
            station=station,
            x_range=[int(X[0]), int(X[-1])],
            copies_per_address=4,
            gate_column_error=max(errors),
            ports=["bright", "dark", "empty", "double"],
            whole_unit_quarter_turn=probe(whole, empty, cal, np.pi / 2, 8),
            whole_unit_full_turn=probe(whole, empty, cal, 2 * np.pi, 8),
            cases=cases,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(data)
