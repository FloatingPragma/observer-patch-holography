"""Two contact spectrometers reconstruct the same signed field per ordinary load."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np
from scipy.optimize import linprog

from . import pair_work as lab

HERE = lab.HERE
pair = lab.energy.pair


def min_variance(w, cells, codes, h, ids, groups):
    variance = (w[:, 0] - 1) ** 2
    rows = np.vstack(
        (pair.constraints(w, cells, codes, h, "energy"), w[:, :4].T, (w[:, 1] - 1) ** 2)
    )
    A = pair.group_matrix(rows, ids, groups)
    b = A.sum(1)
    b[-1] = 80
    c = pair.group_matrix(variance[None, :], ids, groups)[0]
    opt = linprog(c, A_eq=A, b_eq=b, bounds=(0, None))
    assert opt.success
    z = np.array(pair.rational(opt.x), object)
    dual = np.array(pair.rational(opt.eqlin.marginals), object)
    assert (
        all(z >= 0)
        and np.all(A @ z == b)
        and np.all(dual @ A <= c)
        and dual @ b == c @ z
    )
    return (c @ z) / len(w), dict(weights=list(map(str, z)), dual=list(map(str, dual)))


def spectroscopy(spectra, n, t, rng):
    def contrast(left, right):
        levels = spectra[left, right, rng.integers(2, size=len(n))]
        return levels - levels % 8

    return (contrast(n, t) - contrast(n, 0) - contrast(0, t)) // t


def run():
    _, w, cells, _, codes, h, ids, groups = pair.setup()
    board, home, marks, *_ = lab.energy.returned(w)
    spec = lab.local_spectra()
    spectra = np.array([[spec[i, j] for j in range(3)] for i in range(3)])
    weights = dict(
        original=np.full(len(w), F(1, len(w)), object),
        repaired=lab.energy.least_change(cells) / len(w),
    )
    rng = np.random.default_rng(820194)
    bands = {}
    out = {}
    # |P|<=32 for the one-hop signed difference. Each squared reading lies
    # in [0,1024]. Resolve the difference of the two sampled means.
    width = 1024 * np.sqrt(np.log(8 / 0.01) / (2 * lab.SHOTS))
    for t in (1, 2):
        row = {}
        for key, p in weights.items():
            i = rng.choice(len(w), lab.SHOTS, p=np.array(p, float))
            x, y = w[i, :2].T
            field = spectroscopy(spectra, x, t, rng) - spectroscopy(spectra, y, t, rng)
            np.testing.assert_array_equal(field, -16 * (x - y))
            hist = np.bincount(field // 16 + 2, minlength=5)
            exact = p @ (256 * (w[:, 0] - w[:, 1]) ** 2)
            assert p @ (w[:, 0] - w[:, 1]) == 0
            sq = float(np.mean(field**2))
            assert abs(sq - float(exact)) <= width
            row[key] = dict(
                counts=hist.tolist(), mean_square=str(exact), measured_mean_square=sq
            )
        delta = (
            row["repaired"]["measured_mean_square"]
            - row["original"]["measured_mean_square"]
        )
        bands[t] = [delta - 2 * width, delta + 2 * width]
        out[t] = row
    minimum, witness = min_variance(w, cells, codes, h, ids, groups)
    floor = 256 * (minimum - F(11, 28) + 3 * F(1, 14))
    exact_delta = F(out[1]["repaired"]["mean_square"]) - F(
        out[1]["original"]["mean_square"]
    )
    assert 0 < floor < exact_delta and all(
        0 < lo < float(exact_delta) < hi for lo, hi in bands.values()
    )

    clean = sorted(set(range(board.n)) - set(home) - set(np.ravel(marks)))
    occupation = np.zeros((len(w), board.n), int)
    occupation[:, home] = w[:, :4]
    active = []
    nchanged = 0
    for port in clean:
        neighbors = board.neighbours[port]
        loads = occupation[:, neighbors]
        field = -16 * loads.sum(1)
        if np.any(field):
            active.append(port)
        assert (weights["repaired"] - weights["original"]) @ field == 0
        nchanged += int((weights["repaired"] - weights["original"]) @ (field**2) != 0)
        for t in (1, 2):
            # Private padding and all other energy costs shared with an empty
            # ordinary standard at this very port cancel in the difference.
            d = h.diagonal()
            full = (d[3 * t + loads] - d[loads] - d[3 * t]).sum(1)
            parity = 2 * ((t != loads).astype(int) - (loads != 0).astype(int) - 1).sum(
                1
            )
            np.testing.assert_array_equal(full - parity, 4 * t * field)

    held = json.loads((HERE / "pair_transfer.json").read_bytes())
    _, path = lab.mount(board, home, marks)
    deg = int(board.degree[path[1]])
    original = weights["original"]
    static = {}
    for t in (1, 2):
        # Adding an arbitrary static diagonal field translates every signed
        # detuning by the same real d. The nearest resonance may contribute
        # its whole weight. All other resonances are at least half a gap away.
        diag = h.diagonal()
        x, y = w[:, :2].T
        four_detuning = (diag[3 * t + x] - diag[x]) - (diag[3 * t + y] - diag[y])
        sharp = sorted(set(four_detuning))
        spacing = F(int(min(np.diff(sharp))), 4)
        peak = max(sum(original[four_detuning == v]) for v in sharp)
        cap = (
            F(1, 2)
            * (1 + F(1, lab.DURATION))
            * (peak + (1 - peak) / (1 + (spacing / 2) ** 2))
        )
        count = held["ordinary_transfer"][str(t)]["rational_arm_arrivals"]
        observed = count / lab.SHOTS
        lo = observed - np.sqrt(np.log(4 / 0.01) / (2 * lab.SHOTS))
        assert lo > float(cap)
        static[t] = dict(
            peak_weight=str(peak),
            min_detuning_spacing=str(spacing),
            any_independent_static_field_arrival_upper=str(cap),
            measured_tuned_arm_arrival_lower_99=lo,
            degree=deg,
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "PAIR_FIELD.md",
                HERE / "pair_transfer.json",
                HERE / "fluctuation_balls.json",
            )
        },
        signed_levels=[-32, -16, 0, 16, 32],
        shots_per_arm=lab.SHOTS,
        spectrometers=out,
        canonical_mean_square_increase=str(exact_delta),
        observed_increase_bands_99=bands,
        stationary_minimum_partner_variance=str(minimum),
        unavoidable_mean_square_increase=str(floor),
        stationary_minimum_witness=witness,
        original_clean_ports=len(clean),
        original_sensitive_ports=active,
        original_changed_mean_ports=0,
        original_changed_mean_square_ports=nchanged,
        independent_static_field=static,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "pair_field.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        # Convert integer dict keys to their JSON spellings before comparison.
        lab.check(json.loads(dst.read_bytes()), json.loads(json.dumps(result)))
    print(
        result["canonical_mean_square_increase"],
        result["unavoidable_mean_square_increase"],
        result["observed_increase_bands_99"],
        result["independent_static_field"],
    )
