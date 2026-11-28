"""Joint sharp-clock and flight records relative to a mobile-inventory standard."""

import argparse
import json
from hashlib import sha256
from itertools import combinations

import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from scipy.sparse.linalg import expm_multiply

from . import measure as source
from ...matched_probe.signal_time import flight

HERE = source.HERE
SHOTS = 2**28
SEED = 221904
NAMES = ("baseline", "unit", "half_unit", "matched", "unmatched", "pi_6", "pi_3")
SCHEMA = (
    "clock_disagreement",
    "radius_A",
    "radius_B",
    "far_same_radius_A",
    "far_same_radius_B",
)


def inputs():
    states, _, dark, other = source.preparation()
    empty, loaded = (source.timing.get_arms()[0][k] for k in ("empty", "mass+"))
    paired = [np.zeros((2, 2)), states["P"]]

    def mix(a, b):
        return [(x + y) / 2 for x, y in zip(a, b, strict=True)]

    common = mix(empty, paired)
    standard = mix(loaded, paired)
    arms = {"baseline": common, "unit": standard, "half_unit": mix(common, standard)}
    angles = dict(matched=0, unmatched=np.pi / 2, pi_6=np.pi / 6, pi_3=np.pi / 3)
    agreement = {}
    for tag, angle in angles.items():
        v = np.cos(angle) * dark + 1j * np.sin(angle) * other
        rho = states["FG0"] + 0.75 * (np.outer(v, v.conj()) - np.outer(dark, dark))
        agreement[tag] = float(rho[1, 1].real + rho[2, 2].real)
        arms[tag] = mix(empty, [np.zeros((2, 2)), rho])
    return arms, agreement


def joint_laws(stiffness, arms):
    board, _, dist, ports = source.matched.fp.make_board()
    max_radius = int(2 * max(dist))
    laws = {k: np.zeros((3, 2, 2, max_radius + 1)) for k in arms}
    for n in range(3):
        terms, starts, radii = flight.parts(n, stiffness)
        # The same source-sector pinching therefore applies to joint terminal
        # positions and clock signs even with the second indistinguishable load.
        owners = [sum(p < 12 for p in word) for word in combinations(range(board.n), n)]
        charge = np.repeat(owners, 8) - np.tile(np.arange(8) % 2, len(owners))
        for term in terms:
            row, col = term.nonzero()
            np.testing.assert_array_equal(charge[row], charge[col])
        indices = sorted({i for ids in starts.values() for i in ids})
        init = np.zeros((len(radii), len(indices)))
        init[indices, np.arange(len(indices))] = 1
        # One common phase shift for all parities is essential to the clock signs.
        shift = float(terms[0].diagonal().mean())
        for source_parity in (0, 1):
            waves = {}
            for near in (0, 1):
                for far in (0, 1):
                    h = (
                        terms[0]
                        + source_parity * terms[1]
                        + near * terms[2]
                        + far * terms[3]
                    )
                    op = -0.5j * (h - shift * sparse.eye(len(radii), format="csr"))
                    trace = op.diagonal().sum()
                    two = expm_multiply(op, init, traceA=trace)
                    waves[near, far] = (two, expm_multiply(op, two, traceA=trace))
            for time in (0, 1):
                for nf, ff in np.ndindex(2, 2):
                    # In the sign basis, 0 means the initial sign and 1 the other.
                    psi = (
                        sum(
                            (-1) ** (nf * i + ff * j) * w[time]
                            for (i, j), w in waves.items()
                        )
                        / 4
                    )
                    for extra, ids in starts.items():
                        launch = [None, *ports].index(extra)
                        block = psi[:, [indices.index(i) for i in ids]]
                        m = n - int(extra is not None)
                        for r in np.unique(radii):
                            v = block[radii == r]
                            gram = (v.conj().T @ v).reshape(2 + 2 * m, 4, 2 + 2 * m, 4)
                            # Independent fair sharp coins on both clock contacts.
                            effect = np.trace(gram, axis1=1, axis2=3) / 4
                            for k, rho in arms.items():
                                laws[k][launch, nf, ff, r] += (
                                    np.trace(rho[m] @ effect).real / 4
                                )
    for a in laws.values():
        np.testing.assert_allclose(a.sum(axis=(1, 2, 3)), 1, atol=2e-11, rtol=0)
        assert a.min() > -1e-13
    return laws


def decode(histograms):
    h = np.asarray(histograms, float)
    h = h / h.sum(axis=(1, 2, 3), keepdims=True)
    r = np.linspace(0, 1, h.shape[-1])
    radius = (h @ r).sum(axis=(1, 2))
    tagged = (h[:, :, 0] @ r).sum(axis=1)
    return np.r_[
        h[0, 0, 1].sum() + h[0, 1, 0].sum(),
        radius[1:] - radius[0],
        tagged[1:] - tagged[0],
    ]


def witness(changes, errors):
    # Look for one recorded combination with opposite signs on the whole
    # source segment and the positive calibration ray, allowing simultaneous errors.
    a = np.vstack((changes["unit"], -changes["matched"], -changes["unmatched"]))
    b = np.column_stack((-a + errors, a + errors, np.ones(3)))
    result = linprog(
        [0.0] * 10 + [-1.0],
        A_ub=np.vstack((b, [1.0] * 10 + [0.0])),
        b_ub=[0, 0, 0, 1],
        bounds=[(0, None)] * 10 + [(None, None)],
        method="highs",
    )
    if not result.success:
        raise ValueError(result.message)
    weights = result.x[:5] - result.x[5:10]
    values = {k: float(weights @ v) for k, v in changes.items()}
    allowance = float(abs(weights) @ errors)
    return dict(
        weights=weights.tolist(),
        scores=values,
        half_width=allowance,
        separated=bool(
            values["unit"] > allowance
            and max(values["matched"], values["unmatched"]) < -allowance
        ),
    )


def inward_interval(f, g):
    """Range in w, with strict inequalities and endpoints checked separately."""
    low, high = 0.0, 1.0
    for x, y in zip(f[1:3], g[1:3], strict=True):
        if x == y:
            if y >= 0:
                return None
        elif x > y:
            high = min(high, -y / (x - y))
        else:
            low = max(low, -y / (x - y))
    if low >= high:
        return None
    return dict(
        u_bounds=[(1 + 3 * low) / 4, (1 + 3 * high) / 4],
        lower_included=bool(low == 0 and max(g[1:3]) < 0),
        upper_included=bool(high == 1 and max(f[1:3]) < 0),
    )


def run():
    arms, agreement = inputs()
    rng = np.random.default_rng(SEED)
    # Seven bounded raw expectations per source, law and geometry (clock flag
    # and two radius scores at three launch choices), plus four source reads.
    eps = float(np.sqrt(np.log(2 * (2 * len(arms) * 7 + 4) / 0.01) / (2 * SHOTS)))
    errors = eps * np.array([2, 4, 4, 4, 4])
    hits = {k: int(rng.binomial(SHOTS, np.clip(u, 0, 1))) for k, u in agreement.items()}
    out = {}
    for stiffness in (8, 0):
        laws = joint_laws(stiffness, arms)
        counts = {
            k: np.stack(
                [
                    rng.multinomial(
                        SHOTS, np.maximum(bank, 0).ravel() / np.maximum(bank, 0).sum()
                    ).reshape(bank.shape)
                    for bank in arm
                ]
            )
            for k, arm in laws.items()
        }
        averages = {k: decode(a) for k, a in laws.items()}
        observed = {k: decode(a) for k, a in counts.items()}
        exact = {k: v - averages["baseline"] for k, v in averages.items()}
        changes = {k: v - observed["baseline"] for k, v in observed.items()}
        np.testing.assert_allclose(
            laws["half_unit"], (laws["unit"] + laws["baseline"]) / 2, atol=1e-13
        )
        forecasts = {}
        for k in ("pi_6", "pi_3"):
            w = (4 * agreement[k] - 1) / 3
            np.testing.assert_allclose(
                laws[k],
                w * laws["matched"] + (1 - w) * laws["unmatched"],
                rtol=0,
                atol=1e-13,
            )
            weight = (4 * hits[k] / SHOTS - 1) / 3
            pred = weight * observed["matched"] + (1 - weight) * observed["unmatched"]
            widths = errors + (4 / 3) * eps * (
                abs(observed["matched"] - observed["unmatched"]) + errors
            )
            forecasts[k] = dict(
                predicted=pred.tolist(),
                observed=observed[k].tolist(),
                passed=bool(np.all(abs(pred - observed[k]) <= widths)),
            )
        for k in exact:
            assert np.all(abs(changes[k] - exact[k]) < errors)
        out[str(stiffness)] = dict(
            computed_changes={k: v.tolist() for k, v in exact.items()},
            measured_changes={k: v.tolist() for k, v in changes.items()},
            counts={k: v.tolist() for k, v in counts.items()},
            held_forecasts=forecasts,
            certificate=witness(changes, errors),
            both_inward=inward_interval(exact["matched"], exact["unmatched"]),
            certified_both_inward=inward_interval(
                changes["matched"] + errors, changes["unmatched"] + errors
            ),
            half_load_control=bool(
                np.all(
                    abs(changes["half_unit"] - 0.5 * changes["unit"]) <= 1.5 * errors
                )
            ),
        )
    return dict(
        protocol_sha256=sha256((HERE / "JOINT_RECORDS.md").read_bytes()).hexdigest(),
        seed=SEED,
        trials_per_bank=SHOTS,
        labels=SCHEMA,
        launches=[None, *source.matched.fp.make_board()[3]],
        max_total_radius=int(2 * max(source.matched.fp.make_board()[2])),
        radius_unit="fraction of maximum total hop radius",
        source_agreements=hits,
        changes_simultaneous_99_half_width=errors.tolist(),
        laws=out,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "joint_records.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        source.region.check(json.loads(path.read_text()), data)
    for k, law in data["laws"].items():
        print(k, law["certificate"], law["both_inward"], flush=True)
