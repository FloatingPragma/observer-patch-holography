"""Install the heavy source and its one-resident exterior on the same physical board."""

import argparse
from collections import deque
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, splu
from . import strand_exterior as exterior, weak_strand as strand
from .measure import HERE
from ..passive_decay.material import commission
from ..passive_decay.measure import compare_records
from .stationary_edge import classify

DESIGN = "JOINT_EXTERIOR.md"


def assembly(board, dimers, heads, path):
    named = {p for bond in (*dimers, *heads) for p in bond}
    owners = {p // 12 for p in named}
    blocked = {path[-1] // 12}
    lane = set(path[1:-1])
    bridges = {i: set() for i in range(board.cells)}
    for p, q in board.seams:
        if p // 12 != q // 12:
            bridges[p // 12].add(q // 12)
            bridges[q // 12].add(p // 12)
    queue = deque(sorted(owners))
    seen = owners | blocked

    def working():
        return sorted(
            p for i in owners for p in range(12 * i, 12 * i + 12) if p not in lane
        )

    while len(set(working()) - named) < 2 * len(heads) + len(dimers):
        i = queue.popleft()
        for j in sorted(bridges[i] - seen):
            owners.add(j)
            seen.add(j)
            queue.append(j)
            if len(set(working()) - named) >= 2 * len(heads) + len(dimers):
                break
    ports = working()
    idx = {p: j for j, p in enumerate(ports)}
    edges = [(idx[a], idx[b]) for a, b in board.seams if a in idx and b in idx]
    local = SimpleNamespace(n=len(ports), seams=edges, degree=board.degree[ports])
    material = commission(
        local,
        tuple((idx[a], idx[b]) for a, b in dimers),
        {(idx[a], idx[b]): cost for (a, b), cost in heads.items()},
    )
    # Initial homogeneous outside storage is empty except for its separately
    # declared resident. The compact assembly contains every donor for its heads.
    assert material["total_load"] == len(ports) - len(dimers)
    return owners, ports, material


def single(board, sites):
    sites = np.asarray(sites, int)
    ids = np.full(board.n, -1, dtype=int)
    ids[sites] = np.arange(len(sites))
    a, b = np.asarray(board.seams).T
    keep = (ids[a] >= 0) & (ids[b] >= 0)
    x, y = ids[a[keep]], ids[b[keep]]
    degree = np.bincount(np.r_[x, y], minlength=len(sites))
    potential = 8 * (degree + 6 - board.degree[sites])
    h = (
        sparse.diags(potential + degree / 2, format="csr")
        + sparse.coo_matrix(
            (np.full(2 * len(x), -0.5), (np.r_[x, y], np.r_[y, x])),
            shape=(len(sites), len(sites)),
        ).tocsr()
    )
    return h, potential


def inertia(h, shift):
    # Symmetric pivot ordering gives a scalar LDL^T factorization. Verify the
    # symmetry and congruence numerically before counting negative pivots.
    a = (h - sparse.eye(h.shape[0], format="csr") * shift).tocsc()
    lu = splu(
        a,
        permc_spec="MMD_AT_PLUS_A",
        diag_pivot_thresh=0,
        options={"SymmetricMode": True, "Equil": False},
    )
    np.testing.assert_array_equal(lu.perm_r, lu.perm_c)
    diag = lu.U.diagonal()
    assert abs(diag).min() > 1e-9
    order = np.argsort(lu.perm_r)
    defect = a[order, :][:, order] - lu.L @ lu.U
    err = float(np.max(np.asarray(abs(defect).sum(axis=1))))
    symmetry = lu.U - sparse.diags(diag) @ lu.L.T
    assert np.max(abs(symmetry.data), initial=0.0) < 1e-8
    assert err < 2e-9
    return int(np.sum(diag < 0)), err


def approach_certificate(h, potential, corridor):
    rest = np.array([i for i in range(h.shape[0]) if i not in corridor])
    a = h[corridor, :][:, corridor].toarray()
    cross = h[corridor, :][:, rest]
    assert cross.nnz == 1 and cross.data[0] == -0.5
    floor = int(min(potential[rest]))
    # B is its nonnegative principal Laplacian plus its diagonal potential.
    # Interlacing handles a one-port approach. For two ports, below floor
    # the Schur correction is bounded above by ||C||^2/(floor - energy).
    levels = np.linalg.eigvalsh(a)
    norm = abs(float(cross.data[0]))
    assert levels[-1] < floor
    second_lower = (
        floor if len(levels) == 1 else levels[1] - norm**2 / (floor - levels[1])
    )
    return dict(
        approach_block=a.tolist(),
        rest_potential_floor=floor,
        joining_norm=norm,
        ground_upper=float(levels[0]),
        second_lower=float(second_lower),
        first_gap_lower=float(second_lower - levels[0]),
        ground_outside_probability_upper=float((norm / (floor - levels[0])) ** 2),
    )


def eigenpanel(board, sites, path, interval, source, complete):
    h, potential = single(board, sites)
    # Shift-and-invert below the rigorous potential floor locates the cold
    # cluster without a guessed eigenvector or simulation time cutoff.
    low, v = eigsh(
        h,
        k=8,
        sigma=float(min(potential)) - 0.1,
        which="LM",
        v0=np.ones(len(sites)),
        tol=1e-12,
    )
    order = np.argsort(low)
    low, v = low[order], v[:, order]
    ground = float(low[0])
    g = v[:, 0]
    residual = max(
        np.linalg.norm(h @ v[:, j] - low[j] * v[:, j]) for j in range(len(low))
    )
    assert residual < 1e-8
    assert inertia(h, ground - 1e-7)[0] == 0
    assert inertia(h, ground + 1e-7)[0] >= 1
    centre = (interval[0] + interval[1]) / 2
    near, vec = eigsh(
        h,
        k=6,
        sigma=ground + centre,
        which="LM",
        v0=np.cos(np.arange(len(sites))),
        tol=1e-12,
    )
    positive = np.flatnonzero(near - ground > 1e-7)
    j = positive[np.argmin(abs(near[positive] - ground - centre))]
    same = np.abs(near - near[j]) < 1e-7
    # The 1e-7 ground bracket safely expands the two measured line intervals.
    lower = inertia(h, ground + interval[0] - 1e-7)
    upper = inertia(h, ground + interval[1] + 1e-7)
    affordable_count = inertia(h, ground + source["budget"] + 1e-7)[0] - 1
    corridor = [sites.index(p) for p in path[1:-1] if p in sites]
    terminal = sites.index(path[1]) if path[1] in sites else None
    joint = None
    if terminal is not None:
        p = float(abs(g[terminal]) ** 2)
        a = source["incoming_port_probability"]
        # In the separated one-load-per-dimer and one-exterior-load sector,
        # exchange sends an unequal port pair outside the charge coincidence.
        # Its compressed energy adds half a unit to the contrast of eight.
        amplitude = float(
            abs(17 * source["number_transition"] * g[terminal])
            * np.linalg.norm(vec[terminal, same])
        )
        # Exact projective response of an independent site-phase on this
        # separately supplied cold exterior. Sum the first resolved cluster.
        first = abs(low - low[1]) < 1e-7
        assert not first[-1]
        first_phase = float(4 * p * np.sum(v[terminal, first] ** 2))
        joint = dict(
            ground_at_entrance=p,
            sharp_input_opening_work=8.5 * (a * (1 - p) + (1 - a) * p),
            density_coupling_to_nearest_positive_group=amplitude,
            phase_first_excited_probability=first_phase,
            phase_cold_probability=(1 - 2 * p) ** 2,
            idle_first_excited_probability=0.0,
            target_state_count=upper[0] - lower[0],
        )
    _, full_mode = complete
    common = np.array(sites)
    uniform = np.ones(len(sites)) / np.sqrt(len(sites))
    projected = full_mode[common]
    survival = float(projected @ projected)
    projected /= np.sqrt(survival)

    def observed_energy(state):
        energy = float(state @ (h @ state))
        return [energy - ground, float(np.linalg.norm(h @ state - energy * state))]

    kinetic = h - sparse.diags(potential, format="csr", dtype=float)
    connected = sparse.csgraph.connected_components(kinetic, directed=False)[0]
    assert connected == 1
    np.testing.assert_allclose(kinetic @ uniform, 0, atol=1e-14)
    changed, _ = eigsh(
        kinetic,
        k=2,
        sigma=-0.0001,
        which="LM",
        v0=np.cos(np.arange(len(sites))),
        tol=1e-12,
    )
    changed.sort()
    # A uniform mixture of any ground states unresolved at the archive's
    # energy precision makes the stock control invariant under basis changes.
    cold = v[:, abs(low - ground) < 1e-7]
    acted = kinetic @ cold
    mean = float(np.sum(cold * acted) / cold.shape[1])
    spread = float(np.sqrt(np.sum((acted - mean * cold) ** 2) / cold.shape[1]))
    return dict(
        ports=len(sites),
        potential_floor=int(min(potential)),
        ground=ground,
        first_gaps=(low - low[0]).tolist(),
        eigen_residual=residual,
        recoil_interval=interval,
        target_inertia=[lower[0], upper[0]],
        max_factorization_row_error=max(lower[1], upper[1]),
        affordable_count_above_one_ground=affordable_count,
        unresolved_ground_count=inertia(h, ground + 1e-7)[0],
        contrast_removed=dict(
            ground=float(changed[0]),
            first_gap=float(np.diff(changed)[0]),
            true_native_cold_energy=mean,
            true_native_cold_spread=spread,
        ),
        approach_certificate=approach_certificate(h, potential, corridor)
        if corridor
        else None,
        nearest_computed_positive_gap=float(near[j] - ground),
        candidate_residual=float(np.linalg.norm(h @ vec[:, j] - near[j] * vec[:, j])),
        ground_corridor_probability=float(np.dot(g[corridor], g[corridor])),
        opening=joint,
        bare_mode_retained_norm2=survival,
        projected_uniform_energy_spread=observed_energy(uniform),
        projected_bare_mode_energy_spread=observed_energy(projected),
    )


def acquire():
    from ..passive_decay.measure import hamiltonians

    measured = json.loads((HERE / "strand_lines.json").read_text())["compared"]["8"]
    mid = measured[1]["gap"] - measured[0]["gap"]
    uncertainty = measured[0]["radius"] + measured[1]["radius"]
    interval = [mid - uncertainty, mid + uncertainty]
    se, sv = eigh(strand.strand(8).toarray())
    ce = eigh(hamiltonians(4)[0], eigvals_only=True)
    prior_board, prior_dimers, prior_heads, _ = strand.ring_layout(8, True, False)
    reference = {
        label: strand.constructed(prior_dimers, prior_board, prior_heads, label)[0]
        for label in ("native", "source_cut", "no_reversal")
    }
    data = {}
    for level in range(3):
        board = exterior.federation(level)
        layout = strand.ring_layout(8, True, False, board=board)
        if layout is None:
            data[str(board.cells)] = dict(geometry_refused=True)
            continue
        _, dimers, heads, _ = layout
        path, _, _ = exterior.geometry(board, dimers, heads)
        storage, local_ports, prepared = assembly(board, dimers, heads, path)
        for label, ref in reference.items():
            h, _ = strand.constructed(dimers, board, heads, label)
            diff = h - ref
            assert not diff.nnz or max(abs(diff.data)) < 1e-13
        affordable = (se - se[0] > 1e-9) & (se - se[0] < ce[-1] - ce[2])
        pulses = exterior.source_stocks(board, dimers, heads, sv, affordable)
        portbit = strand.fields(8)[0][:, -1]
        if path[0] == dimers[-1][1]:
            portbit = 1 - portbit
        source = dict(
            budget=float(ce[-1] - ce[2]),
            native_lower_core_and_second=float(pulses["native_pulse"][1, 1].real),
            source_cut_lower_core_and_second=float(
                pulses["source_cut_pulse"][1, 1].real
            ),
            native_idle_lower_core_and_second=float(pulses["native_idle"][1, 1].real),
            cut_idle_lower_core_and_second=float(pulses["source_cut_idle"][1, 1].real),
            incoming_port_probability=float((sv[:, 2] ** 2) @ portbit),
            number_transition=float(sv[:, 1] @ (portbit * sv[:, 2])),
        )
        full, _ = single(board, list(range(board.n)))
        ev, w = eigsh(
            full,
            k=5,
            sigma=48 + mid,
            which="LM",
            v0=np.sin(np.arange(board.n) + 1),
            tol=1e-12,
        )
        idx = int(np.argmin(abs(ev - 48 - mid)))
        complete = (float(ev[idx] - 48), w[:, idx])
        bare_counts = [inertia(full, 48 + edge)[0] for edge in interval]
        bare_matches = sorted(
            float(x - 48) for x in ev if interval[0] <= x - 48 <= interval[1]
        )
        assert bare_counts[1] - bare_counts[0] == len(bare_matches)
        outside = [p for p in range(board.n) if p // 12 not in storage]
        phases = {
            "outside": outside,
            "with_approach": sorted(set(outside) | set(path[1:-1])),
        }
        data[str(board.cells)] = dict(
            dimers=dimers,
            heads=[[p, q, x] for (p, q), x in sorted(heads.items())],
            path=path,
            reserved_owners=sorted(storage),
            manufacture_ports=len(local_ports),
            manufacture=prepared,
            source=source,
            prior_bare_gap=complete[0],
            bare_interval_inertia=bare_counts,
            bare_gaps_in_interval=bare_matches,
            exterior={
                label: eigenpanel(board, ports, path, interval, source, complete)
                for label, ports in phases.items()
            },
        )
    return dict(
        protocol_sha256=sha256((HERE / DESIGN).read_bytes()).hexdigest(),
        stationary_protocol_sha256=sha256(
            (HERE / "STATIC_EDGE.md").read_bytes()
        ).hexdigest(),
        stationary_boundary=classify(),
        boards=data,
    )


def decode(result):
    for n, rec in result["boards"].items():
        print(
            n,
            "reserved",
            rec["reserved_owners"],
            "path",
            rec["path"],
            "source",
            rec["source"],
        )
        for label, b in rec["exterior"].items():
            print(
                label,
                "ground, first, affordable count, interval count, approach fraction, native opening",
                b["ground"],
                b["first_gaps"][1],
                b["affordable_count_above_one_ground"],
                b["target_inertia"],
                b["ground_corridor_probability"],
                b["opening"],
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    r = acquire()
    decode(r)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(r, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "joint_exterior.json").read_text()), r)
