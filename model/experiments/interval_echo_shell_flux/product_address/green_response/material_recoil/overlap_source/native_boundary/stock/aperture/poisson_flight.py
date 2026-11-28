"""Predict held material detector records from a weak, independently read source."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from . import weak_source as experiment

lab = experiment.lab
HERE = lab.HERE
READS = 2**54
TIME = 1 / 16


def detector(board, keep):
    h, i, j = lab.lab.ordinary.two_loads(board, lab.lab.K)
    allowed = keep[i] & keep[j]
    index = np.cumsum(keep) - 1
    return h[allowed][:, allowed].tocsr(), index[i[allowed]], index[j[allowed]]


def propagate(h, i, j, factor, shells):
    n = max(i) + 1
    sizes = [int(r.max()) + 1 for r in shells]
    final = [np.zeros(k) for k in sizes]
    wave0 = [np.zeros(k) for k in sizes]
    acceleration = [np.zeros(k) for k in sizes]
    for k in range(0, factor.shape[1], 12):
        initial = factor[j, k : k + 12] / np.sqrt(n - 1)
        hw = h @ initial
        second = 2 * (abs(hw) ** 2 - (initial.conj() * (h @ hw)).real).sum(axis=1)
        wave = expm_multiply((-1j * TIME) * h, initial)
        for s, a, b, c in zip(shells, final, wave0, acceleration, strict=True):
            a += np.bincount(
                s[i], weights=(abs(wave) ** 2).sum(axis=1), minlength=len(a)
            )
            b += np.bincount(
                s[i], weights=(abs(initial) ** 2).sum(axis=1), minlength=len(b)
            )
            c += np.bincount(s[i], weights=second, minlength=len(c))
    return [(a, b, c) for a, b, c in zip(final, wave0, acceleration, strict=True)]


def differentiated(h, i, j, phi, tangent, shell):
    n = len(phi)
    pair = np.column_stack((phi, tangent))[j, :] / np.sqrt(n - 1)
    step = expm_multiply((-1j * TIME) * h, pair)
    change = 2 * ((step[:, 0].conj() * step[:, 1]).real - pair[:, 0] * pair[:, 1])
    hw = h @ pair
    hhw = h @ hw
    curv = 2 * (
        2 * hw[:, 0] * hw[:, 1] - pair[:, 0] * hhw[:, 1] - pair[:, 1] * hhw[:, 0]
    )
    return np.bincount(shell[i], weights=change).cumsum(), np.bincount(
        shell[i], weights=curv
    ).cumsum()


def read_counts(out):
    observed = {}
    error = 4 * np.sqrt(lab.lab.RISK / (2 * READS))
    for name, rows in out["forecast"].items():
        a = out["receivers"][name][name]
        b = out["receivers"]["blank"][name]

        def displacement(bank):
            old, new = [np.array(bank[k]) for k in ("launch", "later")]
            assert old.sum() == new.sum() == READS and min(old) >= 0 and min(new) >= 0
            return (new - old).cumsum() / READS

        change = displacement(a) - displacement(b)
        observed[name] = [
            dict(
                radius=row["radius"],
                change=float(change[row["radius"]]),
                interval=(change[row["radius"]] + np.array([-1, 1]) * error).tolist(),
            )
            for row in rows
        ]
    return observed


def run():
    path = HERE / "weak_source.json"
    reference = HERE / "mediated_release.json"
    source = json.loads(path.read_bytes())
    tape = json.loads(reference.read_bytes())
    fixture = experiment.fixture(tape)
    b, keep, marks, y, h0, ho, hw, _, _ = fixture
    shells = [
        lab.lab.cold.distances(b, [y])[keep],
        lab.lab.cold.distances(b, marks[0])[keep],
    ]
    tags = ("ordinary", "written")
    H, i, j = detector(b, keep)
    p = np.asarray(source["reference"]) / source["accepted"]
    phi = np.sqrt(p)
    forecast = {}
    e, v = eigh(h0)
    psi = experiment.filtered(e, v)[0]
    half = experiment.unitary(h0, experiment.STEP * (1 - experiment.DUTY) / 2)
    gain = lab.ordinary_gain(tape)
    for name, s in zip(tags, shells, strict=True):
        dp, bound, _ = experiment.calibrated(
            source["reference"],
            gain,
            source["filter"]["excited_amplitude_bound"],
            fixture,
            name,
        )
        change, acc = differentiated(H, i, j, phi, dp / (2 * phi), s)
        source_station = marks[0] if name == "written" else [y]
        rows = []
        adj = b.adjacency()[keep][:, keep]
        for r in (2, 3, 4, 6):
            inside = s <= r
            div = adj @ (~inside).astype(int) * inside - (adj @ inside.astype(int)) * (
                ~inside
            )
            # The clean boundary excludes every side-port neighbor and source port.
            boundary = np.where(keep)[0][div != 0]
            if any(x in b.neighbours[y] or x in source_station for x in boundary):
                continue
            weight = (adj @ div) + div
            mid = float(weight @ dp) * gain.mean() / (len(dp) - 1)
            err = (
                float(abs(weight) @ bound) * gain.max()
                + abs(float(weight @ dp)) * np.ptp(gain) / 2
            ) / (len(dp) - 1)
            rows.append(
                dict(
                    radius=r,
                    linear_displacement=float(experiment.DUTY * change[r]),
                    linear_acceleration=float(experiment.DUTY * acc[r]),
                    stationary_source_acceleration=float(experiment.DUTY * mid),
                    source_input_band=(
                        experiment.DUTY * (mid + np.array([-1, 1]) * err)
                    ).tolist(),
                )
            )
        forecast[name] = rows
    rng = np.random.default_rng(135491)
    banks = {}
    laws = {}
    for name, h in (("blank", h0), ("ordinary", ho), ("written", hw)):
        cycle = half @ experiment.unitary(h, experiment.STEP * experiment.DUTY) @ half
        rho = experiment.time_density(cycle, psi)
        vals, vectors = eigh(rho)
        # Discarded probability bounded independently of any position detector.
        retain = vals > 3e-15
        discarded = sum(abs(vals[~retain]))
        assert discarded < 2e-12
        factor = vectors[:, retain] * np.sqrt(vals[retain] / sum(vals[retain]))
        results = propagate(H, i, j, factor, shells)
        arm = {}
        complete = {}
        for tag, (end, zero, acc) in zip(tags, results, strict=True):
            np.testing.assert_allclose([sum(end), sum(zero)], 1.0, atol=2e-12)
            arm[tag] = dict(
                launch=rng.multinomial(READS, zero / zero.sum()).tolist(),
                later=rng.multinomial(READS, end / end.sum()).tolist(),
            )
            complete[tag] = dict(
                displacement=(end - zero).cumsum().tolist(),
                acceleration=acc.cumsum().tolist(),
            )
        banks[name] = arm
        laws[name] = dict(
            surviving_modes=int(retain.sum()),
            discarded_probability=float(discarded),
            cuts=complete,
        )
    out = dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (reference, path, HERE / "POISSON_FLIGHT.md")
        },
        time=TIME,
        trials=READS,
        forecast=forecast,
        receivers=banks,
        laws=laws,
    )
    out["decoded"] = read_counts(out)
    for name in tags:
        for result in out["decoded"][name]:
            idx = result["radius"]
            delta = (
                laws[name]["cuts"][name]["displacement"][idx]
                - laws["blank"]["cuts"][name]["displacement"][idx]
            )
            assert result["interval"][0] <= delta <= result["interval"][1]
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "poisson_flight.json"
    if args.counts_only:
        out = json.loads(path.read_bytes())
        lab.lab.cold.apparatus.check(out["decoded"], read_counts(out))
    else:
        out = run()
        if args.write:
            path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
        else:
            lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(out["forecast"])
    print(out["decoded"])
