"""Find a static weak-dimer material for the unchanged heavy core and follow its emission."""

from collections import deque
from itertools import product
from hashlib import sha256
from pathlib import Path
import argparse
import json
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply, eigsh
from . import preparation
from .measure import HERE
from ..joined_decay.measure import core_pulse
from ..passive_decay.measure import compare_records, hamiltonians, physical
from ..passive_decay.readers import circuit
from ..passive_decay.material import columns
from .material import cross_owner, edge, route
from ..passive_decay.material import commission


def ring_layout(n, multicell, closed=True, board=None):
    if board is None:
        b, d, old, _, _ = route((3,), crossing=True)
    else:
        b = board
        source, _ = cross_owner(b)
        if source is None:
            return None
        d, old = source
    first = (d[-2], d[-1][::-1])
    neighbors = [set(ns) for ns in b.neighbours]
    reserved = {p for bond in (*d[:3], *old) for p in bond} - set(sum(first, ()))
    ring_ports = set(range(b.n)) - reserved
    pairs = sorted(
        (a, z)
        for a in ring_ports
        for z in neighbors[a] & ring_ports
        if b.degree[a] == b.degree[z]
    )
    links = {}
    pair_set = set(pairs)
    for a, z in pairs:
        candidates = {(u, v) for u in neighbors[a] for v in neighbors[a] & neighbors[u]}
        candidates.update(
            (u, v) for v in neighbors[a] & neighbors[z] for u in neighbors[v]
        )
        for u, v in sorted(candidates & pair_set):
            if len({a, z, u, v}) != 4:
                continue
            options = []
            if {u, v} <= neighbors[a]:
                options.append((edge(a, v), edge(a, u), z))
            if {a, z} <= neighbors[v]:
                options.append((edge(a, v), edge(z, v), u))
            links[(a, z), (u, v)] = options
    successors = {p: [] for p in pairs}
    predecessors = {p: [] for p in pairs}
    for a, z in links:
        successors[a].append(z)
        predecessors[z].append(a)
    dist = {first[0]: 0}
    queue = deque([first[0]])
    while queue:
        q = queue.popleft()
        for prev in predecessors[q]:
            if prev not in dist:
                dist[prev] = dist[q] + 1
                queue.append(prev)

    def wiring(seq):
        pairs_to_wire = list(zip(seq[1:], seq[2:]))
        if closed:
            pairs_to_wire.append((seq[-1], seq[0]))
        choices = [links[a, z] for a, z in pairs_to_wire]
        sites = set(sum((*d[:3], *seq), ()))
        for option in product(*choices):
            direct = {b for hi, lo, _ in option for b in (hi, lo)}
            if len(direct) != 2 * len(option) or direct & set(old):
                continue
            heads = dict(old)
            for hi, lo, rail in option:
                heads[hi], heads[lo] = 33, 32
                pins = sorted(
                    z for z in neighbors[rail] - sites if edge(rail, z) not in heads
                )
                if not pins:
                    break
                heads[edge(rail, pins[0])] = 32
            else:
                return heads

    def visit(seq, stock):
        if len(seq) == n:
            if closed and (seq[-1], first[0]) not in links:
                return
            if multicell and all(
                a // 12 == first[0][0] // 12 or z // 12 == first[0][0] // 12
                for a, z in (seq if closed else [seq[-1]])
            ):
                return
            heads = wiring(seq)
            if heads:
                return seq, heads
            return
        for q in successors[seq[-1]]:
            if set(q) & stock or (closed and dist.get(q, n + 1) > n - len(seq)):
                continue
            found = visit((*seq, q), stock | set(q))
            if found:
                return found

    result = visit(first, set(sum(first, ())))
    if result is None:
        return
    seq, heads = result
    dimers = (*d[:3], *seq)
    stock = commission(b, dimers, heads)
    return b, dimers, heads, stock


WAITS = (0, 4, 8, 16, 32)
AFTER = (0, 4, 8)


def fields(n):
    words = np.arange(2**n)
    bits = (words[:, None] >> np.arange(n - 1, -1, -1)) & 1
    flip = sparse.coo_matrix(
        (
            np.full(n * len(words), -0.5),
            (np.tile(words, n), np.concatenate([words ^ (1 << i) for i in range(n)])),
        ),
        shape=(len(words), len(words)),
    ).tocsr()
    return bits, flip


def strand(n):
    bits, flip = fields(n)
    return (
        flip
        + sparse.diags(-(bits[:, :-1] * (1 - bits[:, 1:])).sum(axis=1), dtype=float)
    ).tocsr()


def constructed(dimers, board, heads, name):
    n = len(dimers) - 3
    bits, flip = fields(n + 3)
    a, b, c, x = bits[:, :4].T
    mat = dict(heads)
    coupling = -(1 - c) * x
    if name == "source_cut":
        mat[edge(dimers[2][1], dimers[3][0])] = 32
        coupling = np.zeros_like(coupling)
    if name == "no_reversal":
        mat = dict.fromkeys(heads, 32)
        coupling = np.zeros_like(coupling)
    potential = physical(board, dimers, mat)
    core = 32 * a + 33 * b - 32 * a * b - a * c
    wave = -(bits[:, 3:-1] * (1 - bits[:, 4:])).sum(axis=1)
    prediction = core + wave + coupling
    if name == "no_reversal":
        prediction = 32 * a + 32 * b - 32 * a * b
    np.testing.assert_array_equal(potential - potential[0], prediction - prediction[0])
    return (
        flip + sparse.diags(potential - potential[0], dtype=float)
    ).tocsr(), coupling


def propagate(h, psi, t):
    if t == 0:
        return psi.copy()
    return expm_multiply(-1j * t * h, psi)


def terminal_local(state, port, n):
    # Trace over the complete complement before reading this actual dimer.
    first = np.moveaxis(state.reshape(-1, *((2,) * n)), port + 1, 0).reshape(2, -1)
    rho = first @ first.conj().T
    return [
        float(rho[1, 1].real),
        float(2 * rho[0, 1].real),
        float(-2 * rho[0, 1].imag),
    ]


def probes(b, dimers, heads):
    n = len(dimers) - 3
    ev, frame = eigh(strand(n).toarray())
    gaps = ev - ev[0]
    core, core_basis = eigh(hamiltonians(4)[0])
    budget = core[-1] - core[2]
    target = (gaps > 1e-10) & (gaps <= budget)
    assert target.any() and not np.all(target)
    ground = frame[:, 0]
    band = frame[:, target]
    bits, _ = fields(n)
    launch = (1 - 2 * bits.T) * ground
    # Vacuum projection after a local phase. The release analyzer can transmit
    # all relative phases between the accepted levels.
    absorbed = launch @ band
    delays = np.array(AFTER)
    terminal = np.stack(
        [(absorbed * np.exp(-1j * t * gaps[target])) @ band.T for t in delays]
    )
    known = []
    h = strand(n)
    midpoint = n // 2 - 1
    # Sever a single weak contact: replacing its 33 head by the 32 standard
    # erases the whole orientation term. Use exactly the same input stocks.
    divided = h + sparse.diags(
        bits[:, midpoint] * (1 - bits[:, midpoint + 1]), dtype=float
    )
    for name, generator in (("native", h), ("middle_cut", divided)):
        for t in WAITS:
            p, idle = propagate(generator, np.column_stack((launch[0], ground)), t).T
            tag = band.T @ p
            free = np.abs(launch @ p) ** 2
            conditional = np.abs(absorbed @ tag) ** 2
            raw = terminal_local(p, n - 1, n)
            sham = terminal_local(idle, n - 1, n)
            known.append(
                dict(
                    setting=name,
                    time=t,
                    target=float(np.vdot(tag, tag).real),
                    absorption=free.tolist(),
                    band_absorption=conditional.tolist(),
                    end_read=raw,
                    idle_end_read=sham,
                )
            )
    # Native real state preparation and its inverse implement each spectral
    # reflection. Restrict the existing physical addressing graph to this strand.
    connections, _ = preparation.mount(b, dimers, heads)
    ring_links = {
        (i - 3, j - 3): v for (i, j), v in connections.items() if min(i, j) >= 3
    }
    gates = []
    for i in (0, *np.flatnonzero(target)):
        error, tally = preparation.synthesize(frame[:, i], ring_links)
        gates.append([int(i), float(error), sum(tally.values())])
    return (
        dict(
            parent_levels=(core - core[0]).tolist(),
            max_affordable=float(budget),
            gaps=gaps.tolist(),
            target=np.flatnonzero(target).tolist(),
            sharp_overlap_by_dimer=(absorbed**2).tolist(),
            gates=gates,
            independent_launch=known,
        ),
        core_basis,
        band,
        terminal,
        gaps,
    )


def panel(state, core_basis, band, terminal, coupling, n):
    psi = np.einsum("ab,acx->bcx", core_basis, state.reshape(4, 2, -1))
    population = (abs(psi) ** 2).sum(axis=(1, 2))
    tag = psi @ band
    tagged = (abs(tag) ** 2).sum(axis=(1, 2))
    # Near and remote reads operate on the same initially enrolled component.
    # These are joint probabilities of the selected energy sector, a local
    # phase and the terminal vacuum, on separate copies for each dimer.
    amplitudes = np.einsum("bcx,tdx->tbcd", psi, terminal)
    success = (abs(amplitudes) ** 2).sum(axis=2)
    off_work = -np.vdot(state, coupling * state).real
    return dict(
        core=population.tolist(),
        core_and_band=tagged.tolist(),
        after_core_site_absorption=success.tolist(),
        end_read=terminal_local(state, n - 1, n),
        coupling_removal=float(off_work),
    )


def read_source(layout, calibration):
    b, d, heads, _ = layout
    n = len(d) - 3
    _, core_basis, band, terminal, gaps = calibration
    pulse = circuit(core_pulse(core_basis), 32)
    h, _ = constructed(d, b, heads, "native")
    ev, g = eigsh(h, k=1, which="SA", tol=5e-14, v0=np.ones(len(h.diagonal())))
    ground = g[:, 0]
    residual = np.linalg.norm(h @ ground - ev[0] * ground)
    assert residual < 2e-11
    # Perron positivity fixes the ground sign for a reproducible preparation.
    ground *= np.sign(ground.sum())
    links, _ = preparation.mount(b, d, heads)
    error, tally = preparation.synthesize(ground, links)
    launch = (pulse @ ground.reshape(4, -1)).ravel()
    trials = {}
    for name in ("native", "source_cut", "no_reversal"):
        h, coupling = constructed(d, b, heads, name)
        rows = {arm: [] for arm in ("core_pulse", "idle")}
        for t in WAITS:
            phases = propagate(h, np.column_stack((launch, ground)), t)
            for arm, vec in zip(rows, phases.T, strict=True):
                rows[arm].append(panel(vec, core_basis, band, terminal, coupling, n))
        trials[name] = rows
    return dict(
        ground_residual=float(residual),
        prepared_ground_error=float(error),
        cold_gate_count=sum(tally.values()),
        trials=trials,
    )


def acquire(calibration_only=False):
    results = {}
    for n in (6, 8):
        closed = [ring_layout(n, m) is not None for m in (True, False)]
        assert closed == [False, False]
        layout = ring_layout(n, True, False)
        assert layout is not None
        b, d, heads, stock = layout
        calib = probes(b, d, heads)
        results[str(n)] = dict(
            dimers=d,
            heads=[[a, z, w] for (a, z), w in sorted(heads.items())],
            stock=stock,
            periodic_multicell_and_any=closed,
            calibration=calib[0],
        )
        if not calibration_only:
            results[str(n)]["source"] = read_source(layout, calib)
    return dict(
        protocol_sha256={
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("WEAK_RING.md", "WEAK_STRAND.md")
        },
        column_residual=[columns(), columns(True)],
        settings=results,
    )


def decode(data):
    for n, v in data["settings"].items():
        ref = v["calibration"]
        i = ref["target"]
        print(
            n,
            "ground gaps",
            np.array(ref["gaps"])[i],
            "limit",
            ref["max_affordable"],
            "parent",
            ref["parent_levels"],
        )
        for of in ref["independent_launch"]:
            if of["setting"] == "native":
                print(
                    "known",
                    of["time"],
                    np.round(np.array(of["end_read"]) - of["idle_end_read"], 5),
                    of["band_absorption"][-1],
                )
        if "source" not in v:
            continue
        for name, arms in v["source"]["trials"].items():
            for ix, t in enumerate(WAITS):
                p, z = (arms[a][ix] for a in ("core_pulse", "idle"))
                diff = np.array(p["core_and_band"]) - z["core_and_band"]
                joint = (
                    np.array(p["after_core_site_absorption"])
                    - z["after_core_site_absorption"]
                )
                print(
                    name,
                    t,
                    "band/lower, lower terminal after [0,4,8]",
                    np.round(
                        [diff.sum(), diff[1:3].sum(), *joint[:, 1:3, -1].sum(1)], 7
                    ),
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration-only", action="store_true")
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    result = acquire(args.calibration_only)
    decode(result)
    if args.acquire:
        with args.acquire.open("x") as out:
            json.dump(result, out, allow_nan=False, separators=(",", ":"))
            out.write("\n")
    elif not args.calibration_only:
        compare_records(json.loads((HERE / "weak_strand.json").read_text()), result)
