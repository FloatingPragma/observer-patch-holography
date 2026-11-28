"""Compare the actual emitted energy packet with independent bare exterior standards."""

import argparse
from collections import Counter, defaultdict, deque
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from model.executable_model import federation
from model.candidate_laws.interval_echo import candidate
from . import outlet, weak_strand as ws
from .material import edge
from .measure import HERE
from ..joined_decay.measure import core_pulse
from ..passive_decay.measure import compare_records, hamiltonians
from ..passive_decay.material import columns
from ..passive_decay.readers import circuit


def standards():
    standards = []
    levels = eigh(hamiltonians(4)[0], eigvals_only=True)
    affordable = levels[-1] - levels[2]
    for level in range(3):
        b = federation(level)
        a, z = np.asarray(b.seams).T
        outer = a // 12 != z // 12
        ports = np.r_[a[outer], z[outer]]
        cross = np.bincount(ports // 12, minlength=b.cells)
        assert np.all(cross == 3)
        e = np.eye(12)
        # Every individual owner must have the same internal gap, independent
        # of port numbering and of whether a held refinement has new bridges.
        floor = []
        cell_lines = Counter()
        for i in range(b.cells):
            internal = (a // 12 == i) & (z // 12 == i)
            inc = e[a[internal] % 12] - e[z[internal] % 12]
            assert inc.shape == (30, 12)
            h = inc.T @ inc / 2
            floor.append(eigh(h, eigvals_only=True)[1])
            # Independently inventory every isolated receiver mount, including
            # distant owner types unavailable to the original two short paths.
            degrees = b.degree[12 * i : 12 * i + 12]
            one = eigh(h + np.diag(8 * (5 + 6 - degrees)), eigvals_only=True)
            gaps = one - one[0]
            soft = gaps[(gaps > 1e-10) & (gaps <= affordable)]
            cell_lines[tuple(np.round(soft, 10))] += 1
        single_floor = (5 - np.sqrt(5)) / 2
        np.testing.assert_allclose(floor, single_floor, atol=5e-14)
        # On owner constants, every bridge adds |c_u-c_v|^2/(2*12).
        # This Rayleigh quotient is at most 3/12. All bridges are positive.
        max_low = float(max(cross) / 12)
        direct = None
        if level == 0:
            h = np.diag(0.5 * b.degree.astype(float))
            h[a, z] = h[z, a] = -0.5
            spectrum = eigh(h, eigvals_only=True)
            direct = [
                float(spectrum[0]),
                float(spectrum[b.cells - 1]),
                float(spectrum[b.cells]),
            ]
            assert direct[1] <= max_low and direct[2] >= single_floor - 1e-12
        # For a single bare load all contrast + private energies add to 8*6.
        np.testing.assert_equal(8 * b.degree + 8 * (6 - b.degree), np.full(b.n, 48))
        # The exact one-vacancy/excess potential is >= 80 and exchange is PSD.
        assert min(b.degree) == 5 and max(b.degree) == 6
        v = next(p for p in range(b.n) if b.degree[p] == 5)
        x = next(p for p in range(b.n) if b.degree[p] == 6 and p not in b.neighbours[v])
        numer = np.ones(b.n, dtype=int)
        numer[v], numer[x] = 0, 2
        observed = 8 * (
            (numer * numer - 1) @ (6 - b.degree) + ((numer[a] - numer[z]) ** 2).sum()
        )
        assert observed == 80
        standards.append(
            dict(
                cells=b.cells,
                ports=b.n,
                one_resident_ground=48,
                low_count=b.cells,
                upper_low=max_low,
                next_lower=single_floor,
                smallest_board_direct=direct,
                isolated_affordable_census=[
                    dict(gaps=list(key), cells=count)
                    for key, count in sorted(cell_lines.items())
                ],
                one_vacancy_potential_floor=int(np.min(8 * (2 * b.degree - 6))),
                one_pair_potential_floor=int(observed),
                pair_witness=[v, x],
            )
        )
    return standards


def geometry(board, dimers, heads):
    reserved = {p for bond in (*dimers, *heads) for p in bond}
    busy_owners = {p // 12 for p in reserved}
    queue = deque([(p,) for p in sorted(dimers[-1])])
    seen = set(dimers[-1])
    while queue:
        path = queue.popleft()
        if path[-1] // 12 not in busy_owners:
            break
        for q in sorted(board.neighbours[path[-1]]):
            if q not in reserved | seen:
                seen.add(q)
                queue.append((*path, q))
    else:
        raise ValueError("No separate exterior cell")
    cell = tuple(range(12 * (path[-1] // 12), 12 * (path[-1] // 12 + 1)))
    internal = [
        edge(a, z) for a in cell for z in board.neighbours[a] if z in cell and a < z
    ]
    assert not set(cell) & reserved and len(internal) == 30
    return path, cell, internal


def mobile_generator(board, dimers, heads, path, cell, internal, setting, number=2):
    # Keep all static strand heads. The previous core and C,D contact are
    # parked. An exterior unit pin contributes its actual diagonal energy.
    strand = dimers[3:]
    sites = tuple(sorted(set((*strand[-1], *path, *cell))))
    contact_path = {"open": path, "link_cut": path[1:], "parked": ()}[setting]
    contacts = [edge(*strand[-1]), *internal, *zip(contact_path, contact_path[1:])]
    words, tail = outlet.binary_generator(board, sites, contacts, number)
    m = len(words)
    bits, flip = ws.fields(len(strand) - 1)
    h = sparse.kron(flip, sparse.eye(m), format="csr") + sparse.kron(
        sparse.eye(len(bits)), sparse.csr_matrix(tail), format="csr"
    )
    variable = set(sum(strand, ()))
    earlier = set(sum(dimers[:3], ()))
    # Remove C,D's two triangle heads and the C-side pin (the latter is
    # constant on this subsystem). All purely source heads likewise drop.
    written = {
        pair: cost
        for pair, cost in heads.items()
        if set(pair) & variable and not set(pair) & earlier
    }
    loads = {}
    for i, (a, z) in enumerate(strand[:-1]):
        loads[a] = bits[:, i, None]
        loads[z] = 1 - loads[a]
    for p in strand[-1]:
        loads[p] = np.array([[int(p in word) for word in words]])
    potential = np.zeros((len(bits), m))
    for (a, z), cost in written.items():
        potential += cost * (1 - loads.get(a, 1) * loads.get(z, 1))
    # Matching-degree internal dimers contribute a scalar from their private
    # standards; the mobile ports already have their actual private energies.
    h += sparse.diags((potential - potential[0, 0]).ravel(), format="csr", dtype=float)
    return words, h


def injected(words, dimers, cell, cold, columns):
    lookup = {frozenset(w): i for i, w in enumerate(words)}
    nr, channels = columns.shape[0] // 2, columns.shape[1]
    output = np.zeros((nr, len(words), channels))
    initial = columns.reshape(nr, 2, channels)
    for e in (0, 1):
        for p, amp in zip(cell, cold, strict=True):
            j = lookup[frozenset((dimers[-1][1 - e], p))]
            output[:, j] = amp * initial[:, e]
    np.testing.assert_allclose(
        output.reshape(-1, channels).T @ output.reshape(-1, channels),
        np.eye(channels),
        atol=5e-13,
    )
    return output.reshape(-1, channels)


def measure(state, words, cell, sectors):
    m = len(words)
    wave = state.reshape(-1, m, state.shape[1]).transpose(1, 0, 2)
    groups = defaultdict(list)
    cell = set(cell)
    for j, w in enumerate(words):
        inside = tuple(q for q in w if q in cell)
        outside = tuple(q for q in w if q not in cell)
        groups[len(inside), outside].append((inside, j))
    ncols = state.shape[1]
    scores = {
        k: np.zeros((len(e), ncols, ncols), complex) for k, (_, e, _) in sectors.items()
    }
    for (count, _), seen in groups.items():
        labels, _, frame = sectors[count]
        by_word = dict(seen)
        region = wave[[by_word[w] for w in labels]]
        amps = np.einsum("wk,wsc->ksc", frame, region)
        scores[count] += np.einsum("ksc,ksd->kcd", amps.conj(), amps)
    return scores


def source_stocks(b, d, heads, spatial, affordable, core_levels=slice(1, 3)):
    n = len(d) - 3
    native, _ = ws.constructed(d, b, heads, "native")
    _, v = eigsh(native, k=1, which="SA", tol=5e-14, v0=np.ones(native.shape[0]))
    idle = v[:, 0]
    _, core = eigh(hamiltonians(4)[0])
    active = (circuit(core_pulse(core), 32) @ idle.reshape(4, -1)).ravel()
    band = spatial[:, affordable]
    cohorts = {}
    for name in ("native", "source_cut"):
        h, _ = ws.constructed(d, b, heads, name)
        after = ws.propagate(h, np.column_stack((active, idle)), 32)
        for arm, state in zip(("pulse", "idle"), after.T, strict=True):
            # Retain the original probability of the selected lower core
            # levels and trace both unobserved C orientations.
            actual = np.einsum(
                "ab,acx->bcx", core[:, core_levels], state.reshape(4, 2, -1)
            )
            coeff = actual.reshape(-1, band.shape[0]) @ band
            cohorts[name + "_" + arm] = coeff.T @ coeff.conj()
    bits, _ = ws.fields(n)
    known = band.T @ ((1 - 2 * bits[:, 0]) * spatial[:, 0])
    cohorts["known_D_phase"] = np.outer(known, known)
    return cohorts


def outputs(scores, rho):
    retain, reset = {}, {}
    for k, bank in scores.items():
        retain[k] = np.einsum("kab,ba->k", bank[:, 1:, 1:], rho).real
        reset[k] = bank[:, 0, 0].real * np.trace(rho).real
    norm = np.trace(rho).real
    for outcome in (retain, reset):
        np.testing.assert_allclose(
            sum(x.sum() for x in outcome.values()), norm, atol=3e-10
        )
    return dict(
        retain={str(k): x.tolist() for k, x in retain.items()},
        reset={str(k): x.tolist() for k, x in reset.items()},
    )


def one_board(n):
    b, d, heads, _ = ws.ring_layout(n, True, False)
    path, cell, internal = geometry(b, d, heads)
    sectors = outlet.receiving_cell(b, cell, internal)
    e, g = eigh(ws.strand(n).toarray())
    core = eigh(hamiltonians(4)[0], eigvals_only=True)
    budget = core[-1] - core[2]
    affordable = (e - e[0] > 1e-10) & (e - e[0] <= budget)
    modes = np.column_stack((g[:, 0], g[:, affordable]))
    receiving = sectors[1][1] - sectors[1][1][0]
    limited = np.flatnonzero((receiving > 1e-10) & (receiving <= budget))
    stocks = source_stocks(b, d, heads, g, affordable)
    _, _, cold = sectors[1]
    apparatus = {}
    stages = {}
    analyzer_bounds = {}
    for key in ("open", "link_cut"):
        words, h = mobile_generator(b, d, heads, path, cell, internal, key)
        loaded = injected(words, d, cell, cold[:, 0], modes)
        apparatus[key] = h
        analyzer_bounds[key] = []
        pages = []
        current = loaded
        previous = 0
        for t in ws.WAITS:
            current = ws.propagate(h, current, t - previous)
            previous = t
            score = measure(current, words, cell, sectors)
            soft = score[1][limited].sum(axis=0)
            contrast = soft[1:, 1:] - np.eye(soft.shape[0] - 1) * soft[0, 0]
            extremes = eigh(contrast, eigvals_only=True)
            analyzer_bounds[key].append([float(extremes[0]), float(extremes[-1])])
            pages.append({name: outputs(score, rho) for name, rho in stocks.items()})
        stages[key] = pages
    _, parked = mobile_generator(b, d, heads, path, cell, internal, "parked")
    # Before opening the path, every enrolled strand eigenvector tensored with
    # a cold cell is an eigenvector in the full physical transfer Hilbert space.
    isolated = loaded.T @ (parked @ loaded)
    np.testing.assert_allclose(
        isolated - isolated[0, 0] * np.eye(len(isolated)),
        np.diag(np.r_[0, e[affordable] - e[0]]),
        atol=5e-12,
    )
    np.testing.assert_allclose(parked @ loaded, loaded @ isolated, atol=1e-11)
    toggle = {}
    for setting, h in apparatus.items():
        delta = loaded.T @ ((h - parked) @ loaded)
        toggle[setting] = {
            k: dict(
                keep=float(np.trace(delta[1:, 1:] @ r).real),
                reset=float(delta[0, 0] * np.trace(r).real),
            )
            for k, r in stocks.items()
        }
    return dict(
        path=path,
        cell=cell,
        dimension=len(loaded),
        strand_lines=(e[affordable] - e[0]).tolist(),
        receiving_gaps=receiving.tolist(),
        receiving_affordable=limited.tolist(),
        retained_probabilities={k: float(np.trace(v).real) for k, v in stocks.items()},
        shutter_work=toggle,
        all_enrolled_state_contrast_bounds=analyzer_bounds,
        settings=stages,
    )


def acquire():
    lines = json.loads((HERE / "strand_lines.json").read_text())
    certified = {}
    for n, rows in lines["compared"].items():
        intervals = [
            [r["gap"] - r["radius"], r["gap"] + r["radius"]] for r in rows[:-1]
        ]
        assert all(lo > 0.25 and hi < (5 - np.sqrt(5)) / 2 for lo, hi in intervals)
        certified[n] = intervals
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    local = exchange + candidate.KAPPA * contrast
    error = 0.0
    for i, (a, b) in enumerate(law.SRC):
        expected = np.zeros(len(local))
        expected[i] = candidate.KAPPA * (a - b) ** 2 + 0.5
        expected[law.SRC.index((b, a))] -= 0.5
        error = max(error, float(np.linalg.norm(local[:, i] - expected)))
    assert error < 1e-12
    return dict(
        blank_column_residual=error,
        protocol_sha256=sha256((HERE / "STRAND_EXTERIOR.md").read_bytes()).hexdigest(),
        source_intervals=certified,
        blank_standards=standards(),
        checked_written_columns=[columns(), columns(True)],
        waits=list(ws.WAITS),
        panels={str(n): one_board(n) for n in (6, 8)},
    )


def finite_soft_contrast(panel, observation):
    gaps = np.asarray(panel["receiving_gaps"])
    levels = eigh(hamiltonians(4)[0], eigvals_only=True)
    limit = levels[-1] - levels[2]
    weight = outlet.clock(gaps, 0, limit, labels=32768, low_edge=outlet.PERIOD / 32768)[
        :, 1
    ]
    return float(
        weight @ (np.array(observation["retain"]["1"]) - observation["reset"]["1"])
    )


def compare_exterior(saved, observed):
    # Ordinary-energy readings resolve degenerate eigenspaces. LAPACK may
    # choose different bases inside the same space on different systems.
    for name in saved["panels"]:
        board, dimers, heads, _ = ws.ring_layout(int(name), True, False)
        _, cell, internal = geometry(board, dimers, heads)
        groups = {}
        for count, (_, values, _) in outlet.receiving_cell(
            board, cell, internal
        ).items():
            groups[str(count)] = np.r_[0, np.flatnonzero(np.diff(values) > 1e-8) + 1]
        for record in (saved, observed):
            for pages in record["panels"][name]["settings"].values():
                for page in pages:
                    for trials in page.values():
                        for arm in trials.values():
                            for count, weights in arm.items():
                                arm[count] = np.add.reduceat(
                                    weights, groups[count]
                                ).tolist()
    compare_records(saved, observed)


def decode(record):
    for n, r in record["panels"].items():
        print(
            n,
            "path",
            r["path"],
            "dimension",
            r["dimension"],
            "positive cell lines",
            r["receiving_gaps"][:4],
        )
        low = r["receiving_affordable"]
        print(
            "existing finite-clock terminal retain-minus-reset",
            finite_soft_contrast(r, r["settings"]["open"][-1]["native_pulse"]),
        )
        for setting, pages in r["settings"].items():
            d = []
            for page in pages:
                d.append(
                    {
                        k: [
                            float(np.array(v[a]["1"])[low].sum())
                            for a in ("retain", "reset")
                        ]
                        for k, v in page.items()
                    }
                )
            print(
                setting,
                "native pulse keep/reset by time",
                np.round([v["native_pulse"] for v in d], 8),
                "final all cases",
                d[-1],
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    data = acquire()
    decode(data)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(data, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_exterior(json.loads((HERE / "strand_exterior.json").read_text()), data)
