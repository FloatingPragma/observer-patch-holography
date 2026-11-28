"""Test local round-trip clock periods against native releases on their own ruler."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from . import sector_return_clocks as clocks

HERE = clocks.HERE
HOMES = [8, 20, 32, 44]


def cycles(K, station, neighbors, rng, count=8192):
    """Run individual slot walks until both legs finish, including all idle slots."""
    tickets = np.array([np.repeat(np.arange(len(K)), row) for row in K])
    marks = rng.choice(neighbors, size=count)
    ids = np.arange(count)
    pos = np.full(count, station)
    done_out = np.zeros(count, bool)
    duration = np.zeros(count, int)
    t = 0
    while len(ids):
        t += 1
        pos = tickets[pos, rng.integers(36, size=len(pos))]
        done_out |= pos == marks[ids]
        finished = done_out & (pos == station)
        duration[ids[finished]] = t
        pos, ids, done_out = pos[~finished], ids[~finished], done_out[~finished]
    return {
        "cycles": count,
        "total_slots": int(duration.sum()),
        "mean": float(duration.mean()),
        "standard_error": float(duration.std(ddof=1) / np.sqrt(count)),
        "longest": t,
    }


def transfer(G, u, v, affected):
    a, b = u[affected], v[affected]
    columns = G[:, a] - G[:, b]
    return columns[u] - columns[v]


def predict_edges(transfer_currents, affected, periods):
    """Woodbury update from the ordinary four-terminal transfer measurements."""
    root = np.sqrt(1 - 1 / periods[affected])
    w = transfer_currents * root
    core = np.eye(len(root)) - root[:, None] * transfer_currents[affected] * root
    return np.einsum("ij,ji->i", w, np.linalg.solve(core, w.T))


def clock_experiment():
    (ports, clear, adj, uv, cfg, ids, periods), refs, mu, delta, error = (
        clocks.pairs.apparatus()
    )
    n = len(ports)
    u, v = uv.T
    deg = adj.sum(1)
    homes = np.searchsorted(ports, HOMES)
    affected = np.flatnonzero(np.isin(u, homes) | np.isin(v, homes))
    base = refs["ordinary"][0]
    empty = base * (1 - cfg[:, 1])
    empty /= empty.sum()
    baseline = {
        "cross": np.outer(mu, mu),
        "full": np.outer(base, base),
        "ordinary": np.outer(empty, base),
    }
    ctrl = json.loads((HERE / "sector_register_balls.json").read_text())[
        "matched_local_clock_control"
    ]
    delta["controlled"] = clocks.signed_arms(cfg, ids, ctrl)
    baseline["controlled"] = np.zeros_like(delta["controlled"])
    for a, b in ctrl["word_pairs_by_arm"][0]:
        baseline["controlled"][ids[a], ids[b]] += 0.5
    error["controlled"] = 0.0
    member = np.array([periods[tuple(x)] for x in cfg])
    cases, group = np.unique(
        (member[:, None] + member).reshape(-1, len(uv)), axis=0, return_inverse=True
    )
    group = group.ravel()
    _K0, G0, R0, eps0 = clocks.kernels(adj, uv, np.full(len(u), 2))
    inc = np.zeros((len(u), n), int)
    inc[np.arange(len(u)), u] = inc[np.arange(len(u)), v] = 1
    curr = transfer(G0, u, v, affected)
    blank_sum = R0[u, v] @ inc
    near = np.searchsorted(ports, [48, 82])
    board = clocks.pairs.chart.apparatus()[1]
    native_dist = shortest_path(board.adjacency(), unweighted=True)[
        np.ix_(ports, ports)
    ]
    ii, jj = np.nonzero(
        np.triu((clear[:, None] >= 2) & (clear >= 2) & (native_dist >= 3), 1)
    )
    L = csr_matrix(np.diag(deg) - adj)
    costs, excursion, direct, paid = [], [], [], []
    checks = {"transfer": 0.0, "foster": 0.0, "release": 0.0, "hitting_equation": 0.0}
    witnesses = []
    rng = np.random.default_rng(7465319)
    max_err = eps0
    for k, t in enumerate(cases):
        K, G, R, eps = clocks.kernels(adj, uv, t)
        max_err = max(max_err, eps)
        H = 6 * n * (G.diagonal() - G)
        poisson = (np.eye(n) - K / 36) @ H - 1
        np.fill_diagonal(poisson, 0)
        checks["hitting_equation"] = max(
            checks["hitting_equation"], float(abs(poisson).max())
        )
        field = R - R0
        edge = field[u, v]
        predicted = predict_edges(curr, affected, t / 2)
        checks["transfer"] = max(checks["transfer"], float(abs(predicted - edge).max()))
        price = (1 - 2 / t) @ R[u, v]
        checks["foster"] = max(checks["foster"], float(abs(edge.sum() - price)))
        assert min(edge) >= -4 * (eps + eps0) and np.all(blank_sum >= 1 - 1e-10)
        assert price <= (t / 2 - 1).sum() + 1e-9
        incident = edge @ inc
        flight = clocks.tide.project(L, field)[ii, jj]
        from_ticks = (incident[ii] + incident[jj]) / 2
        checks["release"] = max(
            checks["release"], float(abs(from_ticks - flight).max())
        )
        costs.append(incident)
        excursion.append(R[u, v])
        direct.append(flight)
        paid.append(price)
        if k in (0, len(cases) - 1):
            for station in (homes[0], near[0], near[1]):
                sample = cycles(K, station, np.flatnonzero(adj[station]), rng)
                mean = 6 * n * (blank_sum[station] + incident[station]) / deg[station]
                witnesses.append(
                    dict(
                        case=k,
                        station=int(ports[station]),
                        expected=float(mean),
                        **sample,
                    )
                )
    assert max(checks.values()) < 1e-8
    costs = np.array(costs)
    sums = blank_sum + costs
    canonical = 6 * n * sums / deg
    flight = np.array(direct)
    readings = json.loads((HERE / "sector_return_clock_lab.json").read_text())
    np.testing.assert_array_equal(cfg, np.array(readings["reported_periods"]) - 1)
    emp, _, jk = clocks.lab.estimates(readings["peer"], readings["size"])
    null, _, null_jk = clocks.lab.estimates(readings["ordinary"], readings["size"])
    gain = (6 * n * costs / deg)[group].reshape(len(cfg), len(cfg), n)
    sampled = {
        "peer": clocks.jackknife(gain[:, :, near], emp, jk).tolist(),
        "ordinary_null": clocks.jackknife(gain[:, :, near], null, null_jk).tolist(),
    }

    def paired_redshift(blocks):
        reg = sum(block["registers"] for block in blocks)
        counts = np.sum([b["total_counts"] for b in blocks], axis=0)
        products = np.sum([b["sum_of_count_products"] for b in blocks], axis=0)
        size = readings["size"]
        d, _ = clocks.lab.pair_law(reg, counts, products, size)
        ref = (np.outer(counts, counts) - products) / (reg * (reg - 1) * size**2)
        table = canonical[group][:, near].reshape(len(cfg), len(cfg), 2)
        period = np.einsum("ij,ijk->k", ref, table)
        change = np.einsum("ij,ijk->k", d, gain[:, :, near])
        return float(np.diff(-np.log1p(change / period))[0])

    for arm, blocks in (
        ("peer", readings["peer"]),
        ("ordinary_null", readings["ordinary"]),
    ):
        jack = np.array(
            [paired_redshift(blocks[:i] + blocks[i + 1 :]) for i in range(len(blocks))]
        )
        sampled[arm + "_log_rate_82_relative_to_48"] = [
            paired_redshift(blocks),
            float(
                np.sqrt((len(jack) - 1) / len(jack) * sum((jack - jack.mean()) ** 2))
            ),
        ]
    hop = shortest_path(csr_matrix(adj), unweighted=True).astype(int)
    center = int(np.argmin(hop[:, homes].max(axis=1)))
    out, shells = {}, {}
    for key, D in delta.items():
        d = np.bincount(group, weights=D.ravel(), minlength=len(cases))
        z = np.bincount(group, weights=baseline[key].ravel(), minlength=len(cases))
        joint = z + d
        assert (
            min(joint) >= -1e-15 and abs(z.sum() - 1) < 1e-10 and abs(d.sum()) < 1e-10
        )
        # Normalize the dyadic reference. Charge its adjustment to the source
        # error, retaining the measured contrast when both sums agree.
        joint = np.maximum(joint, 0.0)
        drift = abs(joint.sum() - 1) + abs(z.sum() - 1)
        assert drift < 1e-10
        d = d / joint.sum() + z * ((z.sum() - joint.sum()) / (z.sum() * joint.sum()))
        z = z / z.sum()
        error[key] += drift
        diff = d @ costs
        initial = z @ sums
        lograte = -np.log1p(diff / initial)
        quenched = np.log1p((d @ (1 / canonical)) / (z @ (1 / canonical)))
        path_accel = d @ flight
        width = (
            8 * (max_err + eps0) * abs(d).sum() * deg
            + error[key] * np.ptp(costs, axis=0) / 2
            + 1e-13
        )
        # The derivative of log is at most one on incident sums >=1.
        # Uncertainty in the baseline only enters the difference through
        # the product with diff. The native marginal is certified in full.
        zerr = 2 * error["full"] if key == "cross" else 0.0
        logerr = width + (abs(diff) + width) * (
            4 * max_err * deg + zerr * np.ptp(costs, axis=0)
        )
        bound = float(2 * abs(d) @ np.array(paid))
        assert abs(lograte).sum() <= bound + 1e-10
        clock_a = (diff[ii] + diff[jj]) / 2
        ea = (width[ii] + width[jj]) / 2
        assert abs(clock_a - path_accel).max() < 1e-10
        ext = clear >= 2
        sat = np.abs(lograte[ext]) > logerr[ext]
        edge_contrast = -np.log1p(
            (d @ (np.array(excursion) - R0[u, v])) / (z @ np.array(excursion))
        )
        out[key] = {
            "total_absolute_log_rate": float(abs(lograte).sum()),
            "foster_upper": bound,
            "certified_inward": int(sum(clock_a < -ea)),
            "certified_outward": int(sum(clock_a > ea)),
            "exterior_resolved_clocks": int(sum(sat)),
            "exterior_clock_range": [
                float(lograte[ext].min()),
                float(lograte[ext].max()),
            ],
            "ports_48_82_period_change": (6 * n * diff[near] / deg[near]).tolist(),
            "ports_48_82_log_rate": lograte[near].tolist(),
            "ports_48_82_error": logerr[near].tolist(),
            "quenched_log_rate_at_48_82": quenched[near].tolist(),
            "directed_commute_log_range": [
                float(edge_contrast.min()),
                float(edge_contrast.max()),
            ],
        }
        source_center = int(homes[1]) if key == "ordinary" else center
        rings = hop[source_center]
        divergence = L @ lograte
        derror = deg * logerr + adj @ logerr
        # Ordinary unweighted flux of the log-rate perturbation through
        # each cut. Test constancy only past the last changed home.
        cut_flux = []
        for radius in range(rings[homes].max(), rings.max()):
            lo = rings <= radius
            cut = lo[u] != lo[v]
            signed = np.where(lo[u[cut]], 1, -1)
            flux = float(signed @ (lograte[u[cut]] - lograte[v[cut]]))
            delta_flux = float(sum(logerr[u[cut]] + logerr[v[cut]]))
            assert abs(flux - divergence[lo].sum()) < 1e-12
            cut_flux.append([int(radius), flux, delta_flux])
        out[key]["blank_laplacian_flux"] = {
            "cuts": cut_flux,
            "exterior_nonharmonic": int(sum(ext & (abs(divergence) > derror))),
        }

        def average(x, rings=rings):
            return (np.bincount(rings, weights=x) / np.bincount(rings)).tolist()

        shells[key] = {
            "center": int(ports[source_center]),
            "counts": np.bincount(rings).tolist(),
            "log_rate": average(lograte),
            "error": average(logerr),
        }
    # At clean ports the two-site relative phase slopes follow the ordinary
    # contact block. The scalar energy subtracted from each sector cancels.
    raw = clocks.pairs.chart.source.thermal.pair.stock.kernels()[2]
    active = clear > 0
    clean = clear[active] >= 2
    material = {}
    for load in (1, 2):
        J = clocks.tide.release.doubled_generator(board, ports[active], raw, load)
        kin = board.degree[ports[active]]
        assert np.array_equal(J.diagonal()[clean], kin[clean])
        launch = np.searchsorted(ports[active], [48, 82])
        pair_block = J[np.ix_(launch, launch)] / 2
        y = np.array([[0, -1j], [1j, 0]])
        psi = np.array([1.0, 1.0]) / np.sqrt(2)
        slope = float((1j * psi @ (pair_block @ y - y @ pair_block) @ psi).real)
        assert abs(slope - np.diff(-pair_block.diagonal())[0]) < 1e-14
        material[str(load)] = {
            "clean_sites": int(clean.sum()),
            "native_diagonal_minus_ordinary_degree": [0, 0],
            "initial_two_station_phase_slope": slope,
            "change_in_local_relative_phase_slope": 0,
        }
    return {
        "stations": n,
        "configurations": len(cases),
        "release_pairs": len(ii),
        "input_sha256": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "SECTOR_NEIGHBOR_TICKS.md",
                "sector_register_balls.json",
                "sector_return_clock_lab.json",
            )
        },
        "checks": checks,
        "direct_round_trips": witnesses,
        "source_inferred_period_change_at_48_82": sampled,
        "arms": out,
        "shells": shells,
        "clean_material_phase": material,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    dst = HERE / "sector_neighbor_ticks.json"
    data = clock_experiment()
    if cli.parse_args().write:
        dst.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        clocks.pairs.chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_text()), data
        )
    print(data["checks"], data["arms"]["cross"])
