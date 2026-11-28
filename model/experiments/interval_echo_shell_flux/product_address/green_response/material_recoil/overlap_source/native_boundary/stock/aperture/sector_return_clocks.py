"""Predict the free-probe separation and stationary rates from registered return records."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.linalg import solve
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from . import sector_passage_tide as tide
from . import sector_register_balls as survey
from . import sector_register_pairs as pairs
from . import sector_return_clock_lab as lab
from . import sector_return_sources as sources

HERE = lab.HERE


def kernels(adj, uv, t):
    u, v = uv.T
    W = np.zeros_like(adj)
    assert np.all(12 % t == 0)
    W[u, v] = W[v, u] = 12 // t
    D = np.diag(W.sum(1)) - W
    Z, err = pairs.chart.inverse(D[:-1, :-1])
    # One residual correction on the dyadic inverse. For E=I-DZ the
    # new error is (D^-1-Z.T)E, with an integer certificate for E.
    A = D[:-1, :-1]
    grid = 2**46
    numerator = np.rint(Z * grid).astype(np.int64)
    assert np.array_equal(numerator / grid, Z)
    residual = np.eye(len(Z), dtype=np.int64) * grid - A @ numerator
    leak = abs(residual).sum(0).max() / grid
    corr = Z.T @ (residual / grid)
    Z = Z + corr
    err = err * leak + 4 * np.finfo(float).eps * (
        1 + abs(Z).max() + len(Z) * leak * abs(Z).max()
    )
    ground = np.zeros_like(adj, dtype=float)
    ground[:-1, :-1] = 6 * Z
    # Project the grounded inverse onto the space perpendicular to constants.
    G = ground - ground.mean(1)[:, None] - ground.mean(0)[None, :] + ground.mean()
    d = G.diagonal()
    R = d[:, None] + d - 2 * G
    K = 36 * np.eye(len(adj), dtype=int) - D
    assert min(K.diagonal()) >= 0 and np.all(K.sum(1) == 36) and np.array_equal(K, K.T)
    return K, G, R, 24 * err + 1e-12


def killed_moments(K, i):
    keep = np.arange(len(K)) != i
    # Solve first and second passage-time equations against an absorbing
    # station, independently of the resistance or renewal identity.
    M = 36 * np.eye(keep.sum()) - K[np.ix_(keep, keep)]
    age = solve(M, np.full(keep.sum(), 36.0), assume_a="pos")
    second = solve(M, 36 * (2 * age - 1), assume_a="pos")
    one = 1 + K[i, keep] @ age / 36
    two = 1 + K[i, keep] @ (2 * age + second) / 36
    return one, two - one**2


def signed_arms(cfg, native_ids, raw):
    delta = np.zeros((len(cfg), len(cfg)))
    for sign, arm in zip([-1, 1], raw["word_pairs_by_arm"], strict=True):
        for i, j in arm:
            delta[native_ids[i], native_ids[j]] += sign / 2
    return delta


def jackknife(values, d, jk):
    k = np.asarray(values)
    val = np.einsum("ij,ij...->...", d, k)
    obs = np.einsum("bij,ij...->b...", jk, k)
    err = np.sqrt((len(jk) - 1) / len(jk) * np.sum((obs - obs.mean(0)) ** 2, axis=0))
    return np.stack((val, err), axis=-1)


def run():
    (ports, clear, adj, uv, cfg, ids, periods), refs, _, delta, errors = (
        pairs.apparatus()
    )
    balls = json.loads((HERE / "sector_register_balls.json").read_text())
    full_chart = pairs.chart.apparatus()
    board = full_chart[1]
    physical = shortest_path(board.adjacency(), unweighted=True)[np.ix_(ports, ports)]
    n = len(ports)
    idx = np.triu((clear[:, None] >= 2) & (clear[None, :] >= 2) & (physical >= 3), 1)
    ii, jj = np.nonzero(idx)
    sel = {int(p): int(np.searchsorted(ports, p)) for p in (48, 82)}
    times = np.array([periods[tuple(c)] for c in cfg])
    unique, group = np.unique(
        (times[:, None] + times).reshape(-1, len(uv)), axis=0, return_inverse=True
    )
    delta["controlled"] = signed_arms(cfg, ids, balls["matched_local_clock_control"])
    errors["controlled"] = 0
    obs = json.loads((HERE / "sector_return_clock_lab.json").read_text())
    np.testing.assert_array_equal(np.array(obs["reported_periods"]) - 1, cfg)
    emp, _, jk = lab.estimates(obs["peer"], obs["size"])
    ordinary_emp, _, ordinary_jk = lab.estimates(obs["ordinary"], obs["size"])
    labels = ((cfg[:, None, :] + cfg[None, :, :]) @ np.array([1, 3, 9, 27])).ravel()
    g, r, rl, rh = survey.prices(labels, *refs["retained"])
    _, e, el, eh = survey.prices(labels, *refs["exchanged"])
    q = (r - e)[g].reshape(emp.shape)
    emp_scores = np.bincount(g, weights=emp.ravel())
    qlo = emp_scores.clip(min=0) @ (rl - eh) + emp_scores.clip(max=0) @ (rh - el)
    qhi = emp_scores.clip(min=0) @ (rh - el) + emp_scores.clip(max=0) @ (rl - eh)
    measured_charge = {
        "nominal_and_sampling_se": jackknife(q, emp, jk).tolist(),
        "reference_calibration_range_at_observed_counts": [qlo, qhi],
        "ordinary_null_and_sampling_se": jackknife(
            q, ordinary_emp, ordinary_jk
        ).tolist(),
    }
    L = csr_matrix(np.diag(adj.sum(1)) - adj)
    incident, variances, rates, accel = [], [], [], []
    greens = []
    max_resid = 0.0
    _blankK, blankG, blankR, err0 = kernels(adj, uv, np.full(len(uv), 2))
    bound = 0.0
    killed_error = 0.0
    for t in unique:
        K, G, R, err = kernels(adj, uv, t)
        # Work with differences on the same carriers. The baseline is a
        # common numerical origin for the source/control contrasts.
        field = R - blankR
        inc = (adj * field).sum(1) / 2
        dV = 12 * n**2 * (G - blankG).diagonal()
        a = tide.project(L, field)
        from_clock = -(L @ dV) / (24 * n**2)
        max_resid = max(
            max_resid,
            float(np.max(abs(a[ii, jj] - inc[ii] - inc[jj]))),
            float(np.max(abs(a[ii, jj] - from_clock[ii] - from_clock[jj]))),
        )
        for station in sel.values():
            m, v = killed_moments(K, station)
            vk = 12 * n**2 * G[station, station] + n - n**2
            assert abs(m - n) < 1e-8
            killed_error = max(killed_error, abs(v - vk))
        greens.append(G - blankG)
        incident.append(inc)
        variances.append(dV)
        rates.append((36 - K.diagonal()) / (36 * n))
        accel.append(a[idx])
        bound = max(bound, err)
    assert killed_error < 1e-5 and max_resid < 1e-8
    shape = (len(cfg), len(cfg))
    gain = np.asarray(incident)[group].reshape(*shape, n)
    V = np.asarray(variances)[group].reshape(*shape, n)
    f = np.asarray(rates)[group].reshape(*shape, n)
    exact = np.asarray(accel)[group].reshape(*shape, len(ii))
    controls = {}
    ends = np.column_stack((ports[ii], ports[jj]))
    old = int(np.flatnonzero((ends == [48, 82]).all(1))[0])
    exterior = np.flatnonzero(clear >= 2)
    matrices = np.array(greens)[:, exterior][:, :, exterior]
    # Grounding cancels identically from resistance, so four raw entry
    # errors bound a chord. Each exterior launch uses at most six links.
    chord = 4 * (bound + err0)
    pred = {}
    error_by_pair = {}
    for arm, d in delta.items():
        if arm == "cross":
            label = "distinct_native"
        else:
            label = arm
        q = np.einsum("ij,ijk->k", d, gain)
        pressure = q[ii] + q[jj]
        direct = np.einsum("ij,ijk->k", d, exact)
        mix = np.bincount(group, weights=d.ravel(), minlength=len(matrices))
        average = np.einsum("k,kij->ij", mix, matrices)
        centered = average - average.mean(0) - average.mean(1)[:, None] + average.mean()
        _eig, basis = np.linalg.eigh((centered + centered.T) / 2)
        modes = []
        for k in (0, -1):
            vector = np.rint(1000 * basis[:, k] / abs(basis[:, k]).max()).astype(int)
            vector[np.argmin(abs(vector))] -= vector.sum()
            assert vector.sum() == 0 and abs(vector).sum() > 0
            weights = vector / abs(vector).sum()
            energy = np.einsum("i,kij,j->k", weights, matrices, weights)
            val = float(mix @ energy)
            width = float(
                (bound + err0) * (abs(mix).sum() + errors[arm])
                + np.ptp(energy) * errors[arm] / 2
                + 2e-13
            )
            modes.append(
                {
                    "port_weights": vector.tolist(),
                    "energy_change_at_unit_l1": [val - width, val + width],
                }
            )
        local_error = (
            6 * chord * (abs(d).sum() + errors[arm])
            + errors[arm]
            * (np.ptp(gain, axis=(0, 1))[ii] + np.ptp(gain, axis=(0, 1))[jj])
            / 2
            + 2e-13
        )
        assert np.all(abs(pressure - direct) < local_error)
        pred[arm] = pressure
        error_by_pair[arm] = local_error
        lo, hi = int(pressure.argmin()), int(pressure.argmax())
        picks = [old, lo, hi]
        k = exact[:, :, picks]
        readings = jackknife(k, emp, jk) if arm == "cross" else None
        departure = np.einsum("ij,ijk->k", d, f)
        assert abs(departure[clear >= 2]).max() < 2e-14
        controls[label] = {
            "extremal_zero_sum_port_tests": modes,
            "acceleration_range": [float(pressure[lo]), float(pressure[hi])],
            "maximum_absolute_error": float(max(local_error)),
            "inward": int((pressure < -local_error).sum()),
            "outward": int((pressure > local_error).sum()),
            "unresolved": int((abs(pressure) <= local_error).sum()),
            "pairs_in_order_baseline_min_max": ends[picks].tolist(),
            "acceleration_examples": pressure[picks].tolist(),
            "local_variance_change_at_48_82": (d.ravel() @ V.reshape(-1, n))[
                [sel[48], sel[82]]
            ].tolist(),
            "departure_frequency_change_range": [
                float(departure.min()),
                float(departure.max()),
            ],
        }
        if readings is not None:
            controls[label]["fresh_clock_prediction_and_block_jackknife_se"] = (
                readings.tolist()
            )
            controls[label]["independent_ordinary_null_and_sampling_se"] = jackknife(
                k, ordinary_emp, ordinary_jk
            ).tolist()
    # Resolve which information is missing from individual means. At each
    # home the accepted-move probability on its incident links is 1, 2/3,
    # 1/2 for totals 0, 1, 2. For equal single-member marginals the
    # surviving coefficient of the two-member same-home coincidence is 1/6.
    h = np.searchsorted(ports, [8, 20, 32, 44])
    persistence = cfg[:, None, :] * cfg[None, :, :]
    expected = np.einsum("ij,ijk->k", delta["cross"], persistence)
    observed = jackknife(persistence, emp, jk)
    dep = adj[h].sum(1) / (36 * n) * expected
    direct = np.einsum("ij,ijk->k", delta["cross"], f[:, :, h])
    np.testing.assert_allclose(dep, direct, atol=1e-15, rtol=0)
    # All arms share the same one-member law in this diluted contrast,
    # including the common baseline from the equal-charge experiment.
    balance = json.loads((HERE / "sector_register_escape.json").read_text())[
        "common_baseline_equal_mean_charge"
    ]
    pick = int(np.argmin(pred["controlled"]))
    equal = {}
    for src, key in (("cross", "distinct_native"), ("controlled", "controlled")):
        scale = np.array(balance["interventions"][key]["fraction_of_its_own_half"]) / 2
        band = pred[src][pick] + np.array([-1, 1]) * error_by_pair[src][pick]
        vals = np.multiply.outer(band, scale)
        equal[key] = [float(vals.min()), float(vals.max())]
    assert equal["distinct_native"][0] > 0 > equal["controlled"][1]
    active = clear > 0
    raw = pairs.chart.source.thermal.pair.stock.kernels()[2]
    native_release = {}
    for load in (1, 2):
        J = tide.release.doubled_generator(board, ports[active], raw, load)
        integer, _ = tide.release.derivatives(J)
        tested = clear[active] >= 2
        assert np.array_equal(
            integer[tested] / 4, -L[active][:, active].toarray()[tested] / 2
        )
        native_release[str(load)] = int(tested.sum())
    return {
        "equal_charge_free_release": {
            "endpoints": ends[pick].tolist(),
            "extra_nats": balance["target_increase_nats"],
            "relative_acceleration_change": equal,
        },
        "exterior_ports": ports[exterior].tolist(),
        "fresh_pooled_thermal_charge": measured_charge,
        "source_calibration": sources.calibrate(ports, adj, cfg, exact),
        "swaps_at_fixed_single_home_clock_distributions": {
            "occupancy_classes": cfg.tolist(),
            "endpoints": ends[old].tolist(),
            **sources.compare_swaps(
                cfg,
                (rl - eh)[g].reshape(emp.shape),
                (rh - el)[g].reshape(emp.shape),
                exact[:, :, old],
                6 * chord + 2e-13,
            ),
        },
        "native_release_checked_launches_by_load": native_release,
        "input_sha256": {
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in (
                "SECTOR_RETURN_CLOCKS.md",
                "sector_register_balls.json",
                "sector_register_escape.json",
                "sector_return_clock_lab.json",
            )
        },
        "stations": n,
        "eligible_pairs": len(ii),
        "timed_configurations": len(unique),
        "mean_return_ticks": n,
        "independent_killed_solve_variance_max_discrepancy_ticks2": killed_error,
        "clock_to_release_identity_max_residual": max_resid,
        "controls": controls,
        "same_home_coincidences": {
            "homes": ports[h].tolist(),
            "prediction": expected.tolist(),
            "fresh_estimate_and_standard_error": observed.tolist(),
            "departure_rate_change_per_walk_tick": dep.tolist(),
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    dst = HERE / "sector_return_clocks.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        pairs.chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_text()), data
        )
    print(data["controls"]["distinct_native"], data["same_home_coincidences"])
