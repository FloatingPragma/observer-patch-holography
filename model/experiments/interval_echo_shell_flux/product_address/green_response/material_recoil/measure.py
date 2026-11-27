"""Free receiver flights at the appointments of three native guides."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from .. import measure as green
from . import receiver as body

HERE = Path(__file__).parent
PLAN = json.loads((HERE / "protocol.json").read_text())
N = PLAN["source_releases"]
H = PLAN["maximum_appointments"]
POINTS = np.indices((10, 10, 10)) + 1
WORDS = tuple(product((0, 1), repeat=3))


def configurations():
    all_axes = (0, 1, 2)
    for face in PLAN["faces"]:
        yield "matched", face, 1, 8, all_axes, body.FLIGHT, False
    face = PLAN["faces"][2]
    yield "double", face, 2, 8, all_axes, body.FLIGHT, False
    yield "double_short_flight", face, 2, 8, all_axes, body.FLIGHT / 2, False
    yield "changed_stiffness", face, 1, 4, all_axes, body.FLIGHT, False
    yield "longitudinal", face, 1, 8, (face[-1],), body.FLIGHT, False
    yield "parked", face, 1, 8, all_axes, body.FLIGHT, True


def addresses(face):
    xyz, axis = np.array(face[:3]), face[3]
    for offset in (0, 1):
        centre = xyz.copy()
        centre[axis] += offset
        yield np.tensordot((4, 2, 1), POINTS == centre[:, None, None, None], axes=1)


def visits(rng, releases=N):
    body.pad_routes()
    steps, q = green.measured_steps(None)
    assert np.all(steps[-2:] == 0)
    np.testing.assert_array_equal(q, [1 / 8] * 3)
    live = np.zeros((10,) * 3, dtype=np.int64)
    for corner, n in zip(green.CORNERS, rng.multinomial(releases, [1 / 8] * 8), strict=True):
        measured = body.rails.release_corner(corner)
        live[tuple(p - 1 for p in measured)] = n
    totals = live.copy()
    for tick in range(1, H + 1):
        nxt = np.zeros_like(live)
        # Split whole cohorts, including the two idle labels; no live trial is dropped
        # except at a measured drain. Each split is the multinomial law for iid tags.
        for j, step in enumerate(steps[:-2]):
            take = rng.binomial(live, 1 / (8 - j))
            live -= take
            axis = next(i for i, s in enumerate(step) if s)
            src, dst = [slice(None)] * 3, [slice(None)] * 3
            src[axis], dst[axis] = (
                (slice(1, None), slice(None, -1))
                if step[axis] < 0
                else (slice(None, -1), slice(1, None))
            )
            nxt[tuple(dst)] += take[tuple(src)]
        live += nxt
        if not live.any():
            return totals, tick
        totals += live
    raise ValueError("source count still live at the registered time limit")


def forecast(face, q, k, axes, time, parked, u):
    # Independent Bloch formula averaged over fair masks and orientations.
    d = q * (4 * k + 0.5)
    z = (1 + np.exp(1j * body.CLOCK * d)) / 2
    omega = np.hypot(1, d)
    activity = np.zeros_like(u)
    for sign, indices in zip((-1, 1), addresses(face), strict=True):
        ps = []
        for word in WORDS:
            v = z ** sum(1 - word[i] for i in axes)
            ps.append(
                0
                if parked
                else body.exact_coefficient(q, k, time) * v.real
                + np.sin(omega * time) * v.imag / (2 * omega)
            )
        activity += sign * np.take(ps, indices) / 2
    return float(np.sum(u * activity)), float(np.max(abs(activity)))


def observations(rng, counts, face, q, k, axes, time, parked):
    near = rng.binomial(counts, 0.5)
    rights = 0
    for sign, arm, indices in zip(
        (-1, 1), (near, counts - near), addresses(face), strict=True
    ):
        for i, word in enumerate(WORDS):
            weight = int(arm[indices == i].sum())
            schedules = rng.multinomial(weight, [1 / 8] * 8)
            probs = body.probabilities(q, k, sign, word, axes, time, parked)
            # Freeze the sampling probabilities on a dyadic grid, including exact fair
            # nulls. This prevents library-level changes of side at probability 1/2.
            units = np.asarray(probs) * 2**32
            assert np.min(abs(units - np.floor(units) - 0.5)) > 0.001
            rights += int(rng.binomial(schedules, np.rint(units) / 2**32).sum())
    return rights


def confidence(max_mean):
    # A single coordinate's absorption time bounds the joint lifetime by 120 in mean.
    lifetime = 120
    var = lifetime / 4 + 2 * lifetime**2 * (max_mean + 2**-32) ** 2
    log = np.log(2 * len(list(configurations())) / 0.01)
    bound = np.sqrt(2 * var * log / N) + 2 * H * log / (3 * N)
    eig = 1 - 0.75 * (1 - np.cos(np.pi / 11))
    # The infinite and truncated means differ by at most half the omitted visits.
    return float(
        bound + 0.5 * np.sqrt(1000) * eig**H / (1 - eig) + lifetime / 2**33 + 1e-12
    )


def free_flight(archived=None):
    """Both load values evolve on the full native blank board without clamped pads."""
    from scipy.sparse import csr_matrix, diags
    from scipy.sparse.linalg import expm_multiply

    board, _, _ = body.rails.layout()
    edges = np.array(board.seams)
    left, right = edges.T
    hop_rows = np.r_[left, right]
    hop_cols = np.r_[right, left]
    init = np.zeros(board.n)
    init[body.LEFT] = 1
    matrices = []
    records = []
    rng = np.random.default_rng(91753)
    for load in (1, 2):
        energies, term = body.terms(load, 8)
        onsite = 8 * (6 - board.degree) * load**2
        onsite = onsite.astype(float)
        np.add.at(onsite, left, energies[load, 0] - energies[0, 0])
        np.add.at(onsite, right, energies[0, load] - energies[0, 0])
        h = (
            diags(onsite)
            + csr_matrix(
                (np.full(len(hop_rows), term[0, 1]), (hop_rows, hop_cols)),
                shape=(board.n, board.n),
            )
        ).tocsr()
        matrices.append(h)
        prob = (
            abs(expm_multiply(-3j * (h - h[0, 0] * diags(np.ones(board.n))), init)) ** 2
        )
        counts = rng.multinomial(2**22, prob / prob.sum())
        assert np.max(abs(counts / 2**22 - prob)) < np.sqrt(
            np.log(4 * board.n / 0.01) / (2 * 2**22)
        )
        if archived is not None:
            record = archived["banks"][load - 1]
            assert record["load"] == load and record["starts"] == 2**22
            where = record["where"]
            assert where == sorted(set(where)) and all(0 <= j < board.n for j in where)
            assert len(where) == len(record["counts"]) and all(
                type(n) is int and n > 0 for n in record["counts"]
            )
            observed = np.zeros(board.n, dtype=np.int64)
            observed[where] = record["counts"]
            assert observed.sum() == record["starts"]
            assert np.max(abs(observed / record["starts"] - prob)) < np.sqrt(
                np.log(4 * board.n / 0.01) / (2 * record["starts"])
            )
        records.append(
            dict(
                load=load,
                starts=2**22,
                where=np.flatnonzero(counts).tolist(),
                counts=counts[counts > 0].tolist(),
            )
        )
    diff = matrices[1] - matrices[0] - 144 * diags(np.ones(board.n))
    assert diff.count_nonzero() == 0
    return dict(ports=int(board.n), dwell=3, energy_offset=144, banks=records)


def run():
    rng = np.random.default_rng(PLAN["seed"])
    counts, ticks = visits(rng)
    f = np.zeros_like(counts, dtype=float)
    for c in green.CORNERS:
        f[tuple(j - 1 for j in c)] = 1 / 8
    u = green.poisson(f, (1 / 8,) * 3)[1:-1, 1:-1, 1:-1]
    total = int(counts.sum())
    rows = []
    for name, face, q, k, axes, time, parked in configurations():
        pred, max_mean = forecast(face, q, k, axes, time, parked, u)
        rights = observations(rng, counts, face, q, k, axes, time, parked)
        shift = (rights - total / 2) / N
        radius = confidence(max_mean)
        if abs(shift - pred) > radius:
            raise ArithmeticError("unconditional load count missed its forecast")
        rows.append(
            dict(
                setting=name,
                face=face,
                right_reads=rights,
                shift_per_release=shift,
                predicted=pred,
                simultaneous99=radius,
            )
        )
    return dict(
        free_flight=free_flight(),
        releases=N,
        all_receiver_reads_per_panel=total,
        last_drain=ticks,
        expected_appointments=float(u.sum()),
        source_seed=PLAN["seed"],
        panels=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        old = json.loads(path.read_text())
        flight = old["free_flight"]
        assert len(flight["banks"]) == 2
        assert all(
            flight[k] == data["free_flight"][k]
            for k in ("ports", "dwell", "energy_offset")
        )
        free_flight(flight)
        assert len(old["panels"]) == len(data["panels"])
        for key in (
            "releases",
            "source_seed",
            "all_receiver_reads_per_panel",
            "last_drain",
        ):
            assert old[key] == data[key]
        for o, n in zip(old["panels"], data["panels"], strict=True):
            assert o.keys() == n.keys()
            for key in ("face", "setting", "right_reads"):
                assert o[key] == n[key]
            np.testing.assert_allclose(
                [o[k] for k in ("predicted", "shift_per_release", "simultaneous99")],
                [n[k] for k in ("predicted", "shift_per_release", "simultaneous99")],
                atol=1e-12,
                rtol=0,
            )
        np.testing.assert_allclose(
            data["expected_appointments"],
            old["expected_appointments"],
            atol=1e-12,
            rtol=0,
        )
        print(
            "Measured all signed recoil panels, including shutter and receiver controls."
        )
