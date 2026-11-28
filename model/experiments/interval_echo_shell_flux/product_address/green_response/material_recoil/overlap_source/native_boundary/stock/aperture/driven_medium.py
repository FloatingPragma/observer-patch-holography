"""Exact remote response of a passive written seam to a steady ordinary supply."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import expm
from model.experiments.interval_echo_shell_flux import measure as flow
from . import mediated_release as lab

HERE = lab.HERE


def pulse_calibration():
    # Use the complete native contact with its other contacts held during this
    # addressed appointment, as on the independently surveyed transport clock.
    lab.cold.check_contact(lab.K, 1, 2)
    native = lab.cold.emission
    law = native.LAW
    prop = expm(-1j * np.pi / 8 * native.four_h(lab.K))
    one = law.SRC.index((1, 0))
    other = law.SRC.index((0, 1))
    head = 1 + law.LABELS.index((0, 2))
    blank = float(sum(abs(prop.reshape(11, 9, 99)[:, other, one]) ** 2))
    marked = float(sum(abs(prop.reshape(11, 9, 99)[:, other, 9 * head + one]) ** 2))
    np.testing.assert_allclose([blank, marked], [0.5, 0], atol=1e-13, rtol=0)
    return dict(
        blank_transfer=blank,
        written_transfer=marked,
        offered_neighbor_probability=1 / 8,
    )


def preparation(board, source, shells, mark):
    # flow.stationary takes the source at ordinal zero. The outer survey and
    # each set of physical collector addresses agree across the four banks.
    vertices = [source] + [v for v, r in enumerate(shells) if r <= 7 and v != source]
    index = {v: x for x, v in enumerate(vertices)}
    neighbors = [
        [index[w] for w in board.neighbours[v] if tuple(sorted((v, w))) not in mark]
        if shells[v] < 7
        else []
        for v in vertices
    ]
    cert = flow.stationary(neighbors)
    z = [Fraction(0)] * board.n
    pi = np.zeros(board.n)
    per = Fraction(cert["central"], cert["denominator"])
    for k, v in enumerate(vertices):
        z[v] = Fraction(cert["potential"][k], cert["central"])
        pi[v] = cert["stationary_numerators"][k] / cert["denominator"]
    w = flow.after_burn(neighbors, "native")
    finite = np.zeros(board.n)
    finite[vertices] = [x / flow.GRID for x in w]
    finite[vertices[-1]] += 1 - finite.sum()
    assert min(finite) >= 0
    return z, pi, per, finite


def current_range(count, c, r):
    rad = np.sqrt(lab.RISK / (2 * sum(count)))
    p = np.array(count) / sum(count)
    x = r / sum(count)
    center = c @ p
    width = np.ptp(c) * rad
    return center + np.array([-1, 1]) * width, x + np.array([-1, 1]) * rad


def run():
    board = lab.cold.federation(0)
    record = json.loads((HERE / "mediated_release.json").read_bytes())
    marks, _ = lab.cold.read_map(board, record["spectral_tape"])
    a, b = marks[0]
    shells = lab.cold.distances(board, [a, b])
    rng = np.random.default_rng(293416)
    out = {}
    profiles = {}
    native = {}
    for name, mask in (("blank", []), ("written", marks)):
        h, i, j = lab.ordinary.two_loads(board, lab.K, mask)
        kernel, bank = lab.release_bank(h, i, j, shells)
        native[name] = kernel
        for source in (a, b):
            z, pi, through, finite = preparation(board, source, shells, mask)
            tag = f"{name}_{source}"
            position = rng.multinomial(2**40, finite)
            meter = int(rng.binomial(2**40, finite[shells == 7].sum()))
            prob = bank @ finite
            counts = [
                rng.multinomial(lab.SHOTS, row / row.sum()).tolist() for row in prob
            ]
            out[tag] = dict(
                stationary_throughput=str(through),
                position_counts=position.tolist(),
                source_meter=meter,
                release_counts=counts,
                finite_tv_to_stationarity=float(sum(abs(pi - finite)) / 2),
                cumulative_acceleration=(
                    kernel @ finite / (2 * (board.n - 1))
                ).tolist(),
                grounded_response=[str(x) for x in z],
            )
            profiles[name, source] = z
    np.testing.assert_array_equal(native["blank"][1:], native["written"][1:])
    # The signed source at fixed injected current is the incidence
    # vector of the removed edge times the written-board voltage drop.
    interior = np.flatnonzero(shells < 7)
    signed = {}
    for source in (a, b):
        delta = [
            y - x
            for x, y in zip(
                profiles["blank", source], profiles["written", source], strict=True
            )
        ]
        written = profiles["written", source]
        drop = written[a] - written[b]
        assert drop and (drop > 0) == (source == a)
        div = []
        for v in interior:
            charge = sum(delta[v] - delta[q] for q in board.neighbours[v])
            assert charge == drop * (int(v == a) - int(v == b))
            div.append(charge)
        assert sum(div) == 0
        signed[source] = delta
    za, zb = profiles["written", a], profiles["written", b]
    da, db = za[a] - za[b], zb[a] - zb[b]
    assert all(x * db == y * da for x, y in zip(signed[a], signed[b], strict=True))
    distant = np.flatnonzero((shells >= 4) & (shells < 7))
    largest = max(distant, key=lambda v: abs(signed[a][v]))
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "mediated_release.json", HERE / "DRIVEN_MEDIUM.md")
        },
        seed=293416,
        ordinary_transport=pulse_calibration(),
        edge=[a, b],
        source_trials=2**40,
        position_trials=lab.SHOTS,
        arms=out,
        fixed_flux_drop={str(a): str(da), str(b): str(db)},
        selected_distant_address=int(largest),
        selected_hops=int(shells[largest]),
        distant_affected=sum(bool(signed[a][v]) for v in distant),
        distant_enrolled=len(distant),
        stationary_profile_change_per_injection=[
            str(16 * signed[a][largest]),
            str(16 * signed[b][largest]),
        ],
    )


def decode(data):
    board = lab.cold.federation(0)
    shells = lab.cold.distances(board, data["edge"])
    n = board.n
    target = data["selected_distant_address"]
    rate = np.full(n, 0.0)
    rate[target] = 1.0
    kernel = np.asarray(
        json.loads((HERE / "mediated_release.json").read_bytes())[
            "ordinary_integer_kernel"
        ]
    ) / (2 * (n - 1))
    out = {}
    for v in data["edge"]:
        positions = []
        pred = []
        flight = []
        for name in ("blank", "written"):
            arm = data["arms"][f"{name}_{v}"]
            count = np.array(arm["position_counts"])
            assert sum(count) == data["source_trials"]
            interval, through = current_range(count, rate, arm["source_meter"])
            assert through[0] > 0
            positions.append(
                np.array([interval[0] / through[1], interval[1] / through[0]])
            )
            p = count / sum(count)
            radius = np.sqrt(lab.RISK / (2 * sum(count)))
            central = kernel @ p
            width = np.ptp(kernel, axis=1) * radius
            pred.append(np.stack((central - width, central + width), axis=1))
            assert all(
                sum(row) == data["position_trials"] for row in arm["release_counts"]
            )
            cdf = np.cumsum(arm["release_counts"], axis=1) / data["position_trials"]
            flight.append(cdf - cdf[0])
        delta = flight[1] - flight[0]
        width = 4 * np.sqrt(lab.RISK / (2 * data["position_trials"]))
        out[str(v)] = dict(
            distant_occupation_per_injection_change=lab.delta(
                positions[1], positions[0]
            ).tolist(),
            exterior_acceleration=np.stack(
                (pred[1][1:, 0] - pred[0][1:, 1], pred[1][1:, 1] - pred[0][1:, 0]),
                axis=1,
            ).tolist(),
            probability_changes=delta.tolist(),
            probability_change_error=width,
            resolved_exterior_cuts_at_first_read=[
                r for r in range(1, int(shells.max())) if abs(delta[1, r]) > width
            ],
        )
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    g = cli.add_mutually_exclusive_group()
    g.add_argument("--write", action="store_true")
    g.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "driven_medium.json"
    if args.counts_only:
        data = json.loads(path.read_bytes())
        lab.cold.apparatus.check(data["decoded"], decode(data))
    else:
        data = run()
        data["decoded"] = decode(data)
        if args.write:
            path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        else:
            lab.cold.apparatus.check(json.loads(path.read_bytes()), data)
    print(
        data["selected_hops"],
        data["distant_affected"],
        [
            (k, v["distant_occupation_per_injection_change"])
            for k, v in data["decoded"].items()
        ],
    )
