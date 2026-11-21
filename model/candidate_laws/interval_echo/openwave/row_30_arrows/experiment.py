"""One freely propagating unit and 630 coherent sign-even collars with sharp coins."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation


HERE = Path(__file__).resolve().parent
CAP = 150000
SHOTS = 262144
TIMES = (0.0, 1 / 32, 1 / 16, 1 / 8, 0.25, 0.5, 1.0)
ARMS = ("ascending", "reversed_cell_labels", "descending")


def local_check():
    """All six sign-even, |d|=1, vacuum/one-unit columns and the blank-reference echo."""
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    embeddings = []
    for coin in (0, 1):
        for loads in ((0, 0), (1, 0), (0, 1)):
            u = np.zeros(law.DIM)
            for sign in (-1, 1):
                row = (1 + law.LABELS.index((coin, sign))) * law.N_SRC + law.SRC.index(
                    loads
                )
                u[row] = 1 / np.sqrt(2)
            embeddings.append(u)
    v = np.column_stack(embeddings)
    phi = np.zeros(6)
    phi[[2, 4]] = 1 / np.sqrt(2)
    local = candidate.KAPPA * (np.eye(6) - np.outer(phi, phi))
    err = np.max(abs((h0 + candidate.KAPPA * h1) @ v - v @ local))
    if err > 1e-11:
        raise ArithmeticError("path hop differs from a native retained-collar column")
    incoming = [law.SRC.index(x) for x in [(0, 0), (1, 0)]]
    marked = np.arange(law.DIM) // law.N_SRC
    marked = np.isin(
        marked, [1 + i for i, (_, d) in enumerate(law.LABELS) if abs(d) == 1]
    )
    echo = law.echo()
    read = echo[:, incoming]
    residue = np.max(
        abs(
            echo.T @ (marked[:, None] * read)
            - np.eye(law.DIM)[:, incoming] @ np.diag([0, 1])
        )
    )
    if residue > 1e-11:
        raise ArithmeticError(
            "the blank zero-reference read does not report port occupation"
        )
    return float(err), float(residue)


def sector(car, order):
    rank = np.argsort(order)
    edges = {tuple(pair): i for i, pair in enumerate(car.seams)}
    neighbours = [
        [(w, edges[min(v, w), max(v, w)]) for w in car.neighbours[v]]
        for v in range(car.n)
    ]
    forward = [
        [(w, e) for w, e in row if rank[w] > rank[v]]
        for v, row in enumerate(neighbours)
    ]
    path_counts = [0] * car.n
    path_counts[0] = 1
    longest = [0] * car.n
    for v in order:
        if path_counts[v]:
            for w, _ in forward[v]:
                path_counts[w] += path_counts[v]
                longest[w] = max(longest[w], longest[v] + 1)
    # Each state stores its port and the bitset of the edges whose coins it has flipped.
    ports, parents, masks = [0], [-1], [0]
    layer = range(1)
    depth = 0
    while True:
        size = sum(len(forward[ports[i]]) for i in layer)
        if size == 0 or len(ports) + size > CAP:
            break
        start = len(ports)
        for i in layer:
            for w, e in forward[ports[i]]:
                ports.append(w)
                parents.append(i)
                masks.append(masks[i] ^ (1 << e))
        layer = range(start, len(ports))
        depth += 1
    complement = sum(path_counts) - len(ports)
    if complement < 0 or (complement == 0 and depth != max(longest)):
        raise ArithmeticError("path expansion and independent directed count disagree")
    # An independent labelled-coin check reconstructs all allowed transitions on the
    # materialized states. A flip set fixes the endpoint by divergence of the arrows.
    by_mask = {m: i for i, m in enumerate(masks)}
    if len(by_mask) != len(masks):
        raise ArithmeticError("different paths share a memory configuration")
    degrees = []
    for i, (v, mask) in enumerate(zip(ports, masks, strict=True)):
        eligible = [
            (w, e)
            for w, e in neighbours[v]
            if (rank[w] > rank[v]) != bool(mask & (1 << e))
        ]
        if len(eligible) != len(forward[v]) + (i != 0):
            raise ArithmeticError("arrow degree differs from append/retract degree")
        for w, e in eligible:
            j = by_mask.get(mask ^ (1 << e))
            if j is None and complement and i in layer:
                continue
            if j is None or ports[j] != w or (parents[j] != i and parents[i] != j):
                raise ArithmeticError(
                    "an enabled native hop is missing from the configuration tree"
                )
        degrees.append(len(eligible))

    sites = np.asarray(ports, dtype=int)
    children = np.arange(1, len(sites))
    parent = np.asarray(parents[1:], dtype=int)
    adjacency = sparse.csr_matrix(
        (np.ones(len(children)), (children, parent)), shape=(len(sites),) * 2
    )
    # Remove the common energy 630*kappa. Degree references stay blank at load zero.
    hopping = -candidate.KAPPA / 2
    diagonal = hopping * np.array(degrees) + candidate.KAPPA * (6 - car.degree[sites])
    h = hopping * (adjacency + adjacency.T) + sparse.diags(diagonal, format="csr")
    report = dict(
        reachable_states=sum(path_counts),
        reachable_ports=sum(c > 0 for c in path_counts),
        reachable_cells=sorted({p // 12 for p, c in enumerate(path_counts) if c}),
        paths_per_port=path_counts,
        longest_path=max(longest),
        simulated_states=len(sites),
        simulated_depth=depth,
        omitted_states=complement,
        distinct_coin_words=len(by_mask),
    )
    return h.tocsr(), sites, report


def tail_norm_bound(t, first_missing):
    if t == 0:
        return 0.0
    z = 7 * candidate.KAPPA * t  # Gershgorin bound on shifted H and its compression
    if z >= first_missing + 1:
        return 2.0
    log_term = first_missing * math.log(z) - math.lgamma(first_missing + 1)
    return min(2.0, 2 * math.exp(log_term) / (1 - z / (first_missing + 1)))


def run():
    car = federation(0)
    source_order = list(range(car.n))
    readout_edge = min(e for e in car.external if e[0] < 12)
    target = readout_edge[1] // 12
    readers = np.array(
        [p for p in range(12 * target, 12 * target + 12) if car.degree[p] == 5]
    )
    assert len(readers) == 9
    orders = [
        source_order,
        sorted(source_order, key=lambda p: (p // 12, (-p) % 12)),
        source_order[::-1],
    ]
    err, reference_err = local_check()
    experiments = {}
    for arm, (label, order) in enumerate(zip(ARMS, orders, strict=True)):
        h, sites, report = sector(car, order)
        source = np.zeros(len(sites))
        source[0] = 1.0
        readings = []
        for j, t in enumerate(TIMES):
            bound = (
                tail_norm_bound(t, report["simulated_depth"] + 1)
                if report["omitted_states"]
                else 0.0
            )
            if bound > 0.001:
                readings.append(dict(t=t, norm_bound=bound, resolved=False))
                continue
            amplitudes = expm_multiply(-1j * t * h, source) if t else source
            p = np.bincount(sites, weights=abs(amplitudes) ** 2, minlength=car.n)
            if abs(p.sum() - 1) > 1e-10 or p.min() < 0:
                raise ArithmeticError("the propagated probabilities are invalid")
            quiet = np.ones(car.n, dtype=bool)
            quiet[readers] = False
            outcomes = np.r_[p[readers], p[quiet].sum()]
            rng = np.random.default_rng([232081, arm, j])
            counts = rng.multinomial(SHOTS, outcomes / outcomes.sum())
            prob = float(outcomes[:-1].sum())
            clicked = int(counts[:-1].sum())
            # Simultaneous one-sided Clopper-Pearson lower endpoints across all arms/times.
            lower = (
                float(
                    beta.ppf(
                        1e-6 / (len(orders) * len(TIMES)), clicked, SHOTS - clicked + 1
                    )
                )
                if clicked
                else 0.0
            )
            readings.append(
                dict(
                    t=t,
                    norm_bound=bound,
                    resolved=True,
                    probabilities=outcomes.tolist(),
                    counts=counts.tolist(),
                    patch_click_probability=prob,
                    patch_clicks=clicked,
                    patch_probability_lower_from_counts=max(0.0, lower - bound),
                    exterior_probability=float(p[12:].sum()),
                    resolved_transfer=bool(min(lower, prob) - bound > 0.001),
                )
            )
        report["readings"] = readings
        experiments[label] = report
    return dict(
        source_port=0,
        candidate=candidate.NAME,
        kappa=candidate.KAPPA,
        active_seams=len(car.seams),
        cells=car.cells,
        local_column_error=err,
        reference_read_error=reference_err,
        storage_cap_states=CAP,
        receiver_cell=target,
        receiver_gateway=list(readout_edge),
        reference_ports=readers.tolist(),
        outcomes=readers.tolist() + ["quiet"],
        shots_per_arm_time=SHOTS,
        orientations=experiments,
    )


def canonical(data):
    if isinstance(data, dict):
        return {k: canonical(v) for k, v in data.items()}
    if isinstance(data, list):
        return [canonical(v) for v in data]
    if isinstance(data, float):
        return round(data, 11)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="create a new evidence file")
    output = parser.parse_args().output
    data = canonical(run())
    if output:
        with output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    elif data != json.loads((HERE / "results.json").read_text()):
        raise SystemExit("Retained-arrow receipt differs from the archive")
    for label, panel in data["orientations"].items():
        print(
            label,
            panel["reachable_states"],
            panel["reachable_cells"],
            max(r.get("patch_click_probability", 0) for r in panel["readings"]),
        )
