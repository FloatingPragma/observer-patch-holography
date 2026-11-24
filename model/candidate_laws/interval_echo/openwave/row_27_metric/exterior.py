"""Transport of source information between two disjoint native cells."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path
from scipy.signal import fftconvolve

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls.belt import (
    belts,
)
from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label
from . import experiment as base

HERE = Path(__file__).resolve().parent
SHOTS = 2**25
DATES = np.array([8, 16, 32, 64])
HOLDS = (4, 8, 12, 16)
EPS = float(np.sqrt(np.log(2 * 20000 / 0.01) / (2 * SHOTS)))


def chart():
    board = federation(0)
    a, b = board.external[0]

    def ring(port, position):
        pole_pair, route = next((p, r) for p, r in belts().items() if port % 12 in r)
        order = np.roll(route, position - route.index(port % 12)) + 12 * (port // 12)
        return order.tolist(), [p + 12 * (port // 12) for p in pole_pair]

    mediator, poles = ring(a, 0)
    radar, _ = ring(b, 2)
    pole = next(p for p in poles if a not in board.neighbours[p])
    assert mediator[5] in board.neighbours[pole]
    for perimeter in (mediator, radar):
        assert all(
            v in board.neighbours[u]
            for u, v in zip(perimeter, perimeter[1:] + perimeter[:1])
        )
    distances = shortest_path(board.adjacency(), indices=pole, unweighted=True)[radar]
    assert distances.min() >= 2
    return {
        "pole": pole,
        "mobile_belt": mediator,
        "radar_belt": radar,
        "external_seam": [a, b],
        "native_distances_pole_to_radar": distances.astype(int).tolist(),
    }


def local_occupation_cost(k):
    law = candidate.native_law
    x, y = law.local_term_parts()
    lift = np.kron(label(2, 1, 1)[:, None], np.eye(law.N_SRC))
    h = x + k * y
    energies = np.diag(lift.T @ h @ lift)
    np.testing.assert_allclose(h @ lift, lift * energies, atol=1e-13)
    zero = energies[law.SRC.index((0, 0))]
    return {w: float(energies[i] - zero) for i, w in enumerate(law.SRC)}


def unitary(h, time):
    e, v = eigh(h)
    return (v * np.exp(-1j * e * time)) @ v.T


def belt(block):
    h = np.zeros((10, 10))
    for i in range(10):
        j = [i, (i + 1) % 10]
        h[np.ix_(j, j)] += block
    return h - np.eye(10) * np.trace(h) / 10


def static_cdf(h):
    u = unitary(h, 1)
    out, back = (base.first_click(u, i, j, 64) for i, j in ((0, 5), (5, 0)))
    return np.clip(fftconvolve(out, back).cumsum()[DATES], 0, 1)


def tick_curve(h):
    e, v = eigh(h)
    return np.clip(
        abs(np.exp(-1j * base.WAITS[:, None] * e) @ (v[0] * v[1])) ** 2, 0, 1
    )


def cdf_counts(p, rng):
    freq = rng.multinomial(SHOTS, np.maximum(0, np.diff(np.r_[0, p, 1])).round(13))
    return freq, freq.cumsum()[:-1] / SHOTS


def live_cdf(h, med):
    """The relay keeps the mediator conditional on the full local click subspace."""
    u = unitary(h, 1)
    initial = np.kron(med, np.eye(10)[0])
    out = np.outer(initial, initial.conj())
    back = np.zeros_like(out)
    near = np.arange(0, 100, 10)
    far = near + 5
    total = 0.0
    cdf = []
    for t in range(1, 65):
        out = u @ out @ u.conj().T
        back = u @ back @ u.conj().T
        total += back.diagonal()[near].real.sum()
        back[near, :] = 0
        back[:, near] = 0
        back[np.ix_(far, far)] += out[np.ix_(far, far)]
        out[far, :] = 0
        out[:, far] = 0
        if t in DATES:
            np.testing.assert_allclose(
                out.trace() + back.trace() + total, 1, atol=3e-12
            )
            cdf.append(float(total))
    return np.clip(cdf, 0, 1)


def guards(arm, k, amount):
    """Strict exterior equality for the two original localized changes."""
    blocks, inserts = base.calibration()
    h = base.generator(k, amount, arm, blocks, inserts)
    free = base.generator(k, amount, "sham", blocks, inserts)
    removed = (2, 8) if arm == "coin_even" else (2, 3)
    ext = [v for v in range(10) if v not in removed]
    # Each routine picked its own scalar energy zero. Remove that scalar.
    delta = (h - free)[np.ix_(ext, ext)]
    delta -= np.eye(len(ext)) * np.trace(delta) / len(ext)
    return float(np.max(abs(delta)))


def run():
    rng = np.random.default_rng(271126)
    pool = {}
    blocks, _ = base.calibration()
    all_holds = []
    continuous = []
    calibration = []
    for k in (8.0, 0.0):
        cost = local_occupation_cost(k)
        hm = belt(blocks[k, 1])
        for q in (1, 2):
            hr = belt(blocks[k, q])
            joint = cost[1, q]
            assert abs(cost[(0, q)]) < 1e-12 and abs(cost[1, 0]) < 1e-12
            sensors = []
            for m in (0, 1):
                p = static_cdf(
                    hr + np.diag([joint * m if x == 2 else 0.0 for x in range(10)])
                )
                raw, measured = cdf_counts(p, rng)
                clock = tick_curve(blocks[k, q] + np.diag([0.0, m * joint]))
                c = rng.binomial(SHOTS, clock.round(13))
                key = f"k{k}_q{q}_standard{m}"
                pool[key + "_radar"] = raw
                pool[key + "_clock"] = c
                sensors.append((p, clock, measured, c / SHOTS))
            calibration.append(
                {
                    "k": k,
                    "q": q,
                    "near_radar": sensors[1][2].tolist(),
                    "away_radar": sensors[0][2].tolist(),
                }
            )
            for occupation in (0, 1):
                med_h = hm.copy()
                med_h[5, 5] += cost[occupation, 1]
                h = np.kron(med_h, np.eye(10)) + np.kron(np.eye(10), hr)
                h[2, 2] += joint  # |mediator at 0, radar at 2>
                for launch in (0, 5):
                    for hold in HOLDS if k else (12,):
                        med = unitary(med_h, hold)[:, launch]
                        prob = float(abs(med[0]) ** 2)
                        key = f"k{k}_q{q}_occ{occupation}_launch{launch}_hold{hold}"
                        occ = rng.binomial(SHOTS, round(prob, 13))
                        pool[key + "_arrival"] = np.array([occ, SHOTS - occ])
                        # In the stored arm the position eigenspaces are invariant.
                        # Coherences between positions vanish on the receiver trace.
                        pred = occ / SHOTS
                        near, away = sensors[1], sensors[0]
                        measured_clock = rng.binomial(
                            SHOTS, ((1 - prob) * away[1] + prob * near[1]).round(13)
                        )
                        exact = (1 - prob) * away[0] + prob * near[0]
                        raw, rad = cdf_counts(exact, rng)
                        pool[key + "_stored_clock"] = measured_clock
                        pool[key + "_stored_radar"] = raw
                        forecast = (1 - pred) * away[2] + pred * near[2]
                        clock_forecast = (1 - pred) * away[3] + pred * near[3]
                        radar_error = max(abs(rad - forecast))
                        clock_error = max(abs(measured_clock / SHOTS - clock_forecast))
                        assert max(radar_error, clock_error) < 3 * EPS
                        all_holds.append(
                            {
                                "k": k,
                                "q": q,
                                "source": occupation,
                                "launch": launch,
                                "hold": hold,
                                "arrived_mediator": occ,
                                "clock_forecast_max_error": float(clock_error),
                                "radar_forecast_max_error": float(radar_error),
                                "measured": rad.tolist(),
                                "forecast": forecast.tolist(),
                            }
                        )
                        if hold == 12:
                            raw, rad = cdf_counts(live_cdf(h, med), rng)
                            pool[key + "_continuous_radar"] = raw
                            continuous.append(
                                {
                                    "k": k,
                                    "q": q,
                                    "source": occupation,
                                    "launch": launch,
                                    "count_by_deadline": raw.cumsum()[:-1].tolist(),
                                    "cdf": rad.tolist(),
                                }
                            )
            # No mediator anywhere. The occupied source can change only a scalar.
            vacuum = []
            for occpole in (0, 1):
                raw, x = cdf_counts(static_cdf(hr + np.eye(10) * cost[occpole, 0]), rng)
                pool[f"k{k}_q{q}_empty_mediator_occ{occpole}"] = raw
                vacuum.append(x)
            assert np.max(abs(vacuum[1] - vacuum[0])) < 2 * EPS
    obstruction = {
        arm: [guards(arm, 8, q) for q in (1, 2)] for arm in ("coin_even", "empty_write")
    }
    return {
        "chart": chart(),
        "shots": SHOTS,
        "deadlines": DATES.tolist(),
        "per_probability_error": EPS,
        "mixture_forecast_error_bound": 3 * EPS,
        "source_contrast_error_bound": 2 * EPS,
        "source_free_exterior_generator_residual": obstruction,
        "receiver_standards": calibration,
        "stored": all_holds,
        "continuous": continuous,
    }, pool


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result, bank = run()
    path = HERE / "exterior.json"
    counts = HERE / "exterior.npz"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
        np.savez_compressed(counts, **bank)
    else:
        prior = json.loads(path.read_text())

        def compare(a, b):
            if isinstance(a, dict):
                assert set(a) == set(b)
                for key in a:
                    compare(a[key], b[key])
            elif isinstance(a, list):
                for left, right in zip(a, b, strict=True):
                    compare(left, right)
            else:
                np.testing.assert_allclose(a, b, rtol=1e-7, atol=1e-8)

        compare(prior, result)
        with np.load(counts) as saved:
            assert set(saved.files) == set(bank)
            for key in bank:
                np.testing.assert_array_equal(saved[key], bank[key])
    print(f"Verified {len(bank)} exterior count banks.")


if __name__ == "__main__":
    main()
