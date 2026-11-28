"""Resolve weak-source release from separate position and interference records."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from . import poisson_flight as flight

lab = flight.lab
experiment = flight.experiment
HERE = flight.HERE
BATCH = 2**50
BANKS = 2048
TRIALS = BATCH * BANKS
FRINGES = lab.READS


def meters(board, keep, marks, y, h0, H, i, j):
    s = keep.sum()
    nr = np.arange(s)
    adj = board.adjacency()[keep][:, keep]
    wave = sparse.csr_matrix(
        (np.ones(len(j), dtype=np.int64), (np.arange(len(j)), j)), shape=(len(j), s)
    )
    H2 = (2 * H).astype(np.int64)
    assert (H2 != 2 * H).nnz == 0
    a = H2 @ wave
    b = H2 @ a
    out = {}
    for name, src in (("ordinary", [y]), ("written", marks[0])):
        shell = lab.lab.cold.distances(board, src)[keep]
        arm = []
        for radius in (2, 3, 4, 6):
            inside = shell <= radius
            du, dv = adj[inside][:, ~inside].nonzero()
            u, v = nr[inside][du], nr[~inside][dv]
            if set(np.where(keep)[0][np.r_[u, v]]) & (
                set(board.neighbours[y]) | set(src)
            ):
                continue
            div = np.bincount(u, minlength=s) - np.bincount(v, minlength=s)
            weight = adj @ div + div
            skew = np.zeros_like(h0)
            skew[u, v] = 0.5
            skew[v, u] = -0.5
            correction = h0 @ skew - skew @ h0
            P = inside[i]
            term = wave[P].T @ b[P]
            exact = (2 * (a[P].T @ a[P]) - term - term.T).toarray()
            np.testing.assert_array_equal(
                exact, 4 * (4 * lab.lab.K * np.diag(weight) + correction)
            )
            arm.append((radius, weight, correction))
        out[name] = arm
    return out


def decode(out, src, tape, readers):
    gain = lab.ordinary_gain(tape)
    result = {}
    for name, rows in readers.items():
        bank = []
        for idx, (radius, w, c) in enumerate(rows):
            width = 2 * np.sqrt(lab.lab.RISK / (2 * src["accepted"])) * np.ptp(w)
            value = (
                w
                @ (np.array(src["arms"][name]) - np.array(src["arms"]["blank"]))
                / src["accepted"]
            )
            field = lab.lab.signed_product(gain, value + np.array([-width, width]))
            norm = np.trace(abs(c)) + np.triu(abs(c), 1).sum()
            fringes = []
            for tag in ("blank", name):
                pos, neg, other = out["fringes"][tag][name][idx]
                assert pos + neg + other == FRINGES and min(pos, neg, other) >= 0
                fringes.append((pos - neg) / FRINGES)
            corr = norm * (
                (fringes[1] - fringes[0])
                + np.array([-1, 1]) * 2 * np.sqrt(2 * lab.lab.RISK / FRINGES)
            )
            acceleration = (field + corr) / (len(w) - 1)
            displacement = 0.0
            for tag, sign in ((name, 1), ("blank", -1)):
                for clock, phase in (("later", 1), ("launch", -1)):
                    counts = out["position"][tag][name][clock]
                    assert sum(counts) == TRIALS and min(counts) >= 0
                    displacement += sign * phase * (sum(counts[: radius + 1]) / TRIALS)
            error = 4 * np.sqrt(lab.lab.RISK / (2 * TRIALS))
            bank.append(
                dict(
                    radius=radius,
                    from_source=acceleration.tolist(),
                    from_position=float(displacement),
                    terminal_interval=(
                        displacement + np.array([-1, 1]) * error
                    ).tolist(),
                )
            )
        result[name] = bank
    return result


def run(counts=None):
    srcpath = HERE / "weak_source.json"
    oldpath = HERE / "poisson_flight.json"
    refpath = HERE / "mediated_release.json"
    tape = json.loads(refpath.read_bytes())
    src = json.loads(srcpath.read_bytes())
    old = json.loads(oldpath.read_bytes())
    b, keep, marks, y, h0, ho, hw, _, _ = experiment.fixture(tape)
    H, i, j = flight.detector(b, keep)
    obs = meters(b, keep, marks, y, h0, H, i, j)
    if counts is not None:
        return decode(counts, src, tape, obs)
    e, v = eigh(h0)
    psi = experiment.filtered(e, v)[0]
    half = experiment.unitary(h0, experiment.STEP * (1 - experiment.DUTY) / 2)
    rng = np.random.default_rng(691847)
    position, fringe = {}, {}
    for name, h in (("blank", h0), ("ordinary", ho), ("written", hw)):
        cycle = half @ experiment.unitary(h, experiment.STEP * experiment.DUTY) @ half
        rho = experiment.time_density(cycle, psi)
        rows, fr = {}, {}
        for tag, source in (("ordinary", [y]), ("written", marks[0])):
            shell = lab.lab.cold.distances(b, source)[keep]
            p = np.bincount(shell, weights=(1 - np.diag(rho).real) / (keep.sum() - 1))
            # Independently generated full native terminal probabilities from the first run.
            change = old["laws"][name]["cuts"][tag]["displacement"]
            final = p + np.diff(np.r_[0, change])
            rows[tag] = {
                clock: rng.multinomial(BATCH, prob / prob.sum(), size=BANKS)
                .sum(axis=0)
                .tolist()
                for clock, prob in (("launch", p), ("later", final))
            }
            fr[tag] = [lab.fringe_counts(c, rho, rng) for _, _, c in obs[tag]]
        position[name] = rows
        fringe[name] = fr
    out = dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (srcpath, oldpath, refpath, HERE / "POISSON_PRECISION.md")
        },
        seed=691847,
        batch=BATCH,
        banks=BANKS,
        position=position,
        fringes=fringe,
    )
    out["decoded"] = decode(out, src, tape, obs)
    for name, rows in out["decoded"].items():
        a, b = [old["laws"][tag]["cuts"][name] for tag in (name, "blank")]
        for row in rows:
            r = row["radius"]
            for band, key in (
                ("from_source", "acceleration"),
                ("terminal_interval", "displacement"),
            ):
                val = a[key][r] - b[key][r]
                assert row[band][0] <= val <= row[band][1]
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "poisson_precision.json"
    if args.counts_only:
        out = json.loads(path.read_bytes())
        lab.lab.cold.apparatus.check(out["decoded"], run(out))
    else:
        out = run()
        if args.write:
            path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
        else:
            lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(out["decoded"])
