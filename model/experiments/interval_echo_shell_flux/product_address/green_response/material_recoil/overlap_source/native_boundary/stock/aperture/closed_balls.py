"""Test a single signed rest density on overlapping regional readings of the returned stock."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

from mpmath import iv
import numpy as np

from . import closed_energy as energy
from . import thermal_radar as radar

HERE = energy.HERE
RADII = (2, 4, 6, 8, 10, 12)
N = 2**50
EPS = 2**-23
DYAD = 2**36


def mask_ids(w, mask):
    cols = [i for i in range(7) if mask & (1 << i)]
    return (
        np.unique(w[:, cols], axis=0, return_inverse=True)[1]
        if cols
        else np.zeros(len(w), int)
    )


def score(ids, ref, native, other, weights):
    def aggregate(counts):
        out = np.zeros(ids.max() + 1, np.int64)
        np.add.at(out, ids, counts)
        return out

    multiplicity = aggregate(np.ones(len(ids), int))
    r = aggregate(ref) / N
    delta = (aggregate(other) - aggregate(native)) / N
    # Reference is independent of both endpoint banks. For any fixed
    # score, the two endpoint expectations have Hoeffding halfwidth
    # 2*EPS times its range simultaneously over all possible masks.
    logs = -np.log(r)
    value = float(delta @ logs)
    log_error = max(-np.log1p(-EPS / r))
    absolute_l1 = min(2, sum(abs(delta)) + 2 * len(r) * EPS)
    width = 2 * EPS * np.ptp(logs) + log_error * absolute_l1 + 2e-12
    predicted = 0.0
    for k, count in enumerate(multiplicity):
        predicted -= float(sum(weights[ids == k] - 1) / len(weights)) * np.log(
            count / len(weights)
        )
    return [value - width, value + width], predicted


def energy_at_ports(board, terms, private, weights):
    change = (weights - 1) / len(weights)
    sites = change @ private / 4
    for (u, v), x in zip(board.seams, terms @ change / 8, strict=True):
        sites[u] += x
        sites[v] += x
    return np.array(sites, float)


def qiv(x):
    return iv.mpf(x.numerator) / int(x.denominator)


def geometry_test(hops, lookup, scores, predicted, port_energy):
    # If S_R=K_R*rho/(2R), set d=rho/4; all radius-two centers then
    # fix K2*d=S2. K2 is 4I+3A and has determinant 1 modulo 3.
    table = np.array([scores[m] for m in range(128)])
    k = np.maximum(0, 4 - hops**2)
    assert np.array_equal(k % 3, np.eye(len(k), dtype=int))
    identity = np.eye(len(k), dtype=np.int64)
    pinv = np.rint(np.linalg.solve(k, identity) * DYAD).astype(np.int64)
    res = DYAD * identity - pinv @ k
    err = max(abs(res).sum(1))
    assert err < DYAD
    inverse_bound = Q(int(max(abs(pinv).sum(1))), int(DYAD - err))
    base = table[lookup[2]]
    base_norm = Q(float(abs(base).max()))
    density_norm = inverse_bound * base_norm
    results = {}
    for radius in RADII[1:]:
        kr = np.maximum(0, radius**2 - hops**2)
        weights = (2.0 / radius) * kr @ pinv / DYAD
        target = table[lookup[radius]]
        center = target.mean(1)
        differences = center - weights @ base.mean(1)
        y = np.rint(weights * DYAD).astype(np.int64)
        # Exact integer residual per candidate witness, even though the
        # inverse approximation and row selection use floating point.
        mismatch = (2 * DYAD) * kr - radius * (y @ k)
        row_norm = abs(mismatch).sum(axis=1)
        errors = (abs(y) / DYAD) @ np.diff(base, axis=1).ravel() / 2 + np.diff(
            target, axis=1
        ).ravel() / 2
        errors += row_norm / (radius * DYAD) * float(density_norm) + 2e-9
        margins = abs(differences) - errors
        pick = int(np.argmax(margins))
        yi = y[pick]
        witness = iv.mpf(target[pick].tolist())
        for coeff, b in zip(yi, base, strict=True):
            if coeff:
                witness -= iv.mpf(int(coeff)) / DYAD * iv.mpf(b.tolist())
        residual = Q(int(row_norm[pick]), radius * DYAD) * density_norm
        witness += iv.mpf([-1, 1]) * qiv(residual)
        r = radar.lab.interval_output(witness)
        null_read = np.asarray(predicted)[lookup[radius]]
        # The measured native energy distribution is one particular candidate;
        # it is offered here in native units, before any geometric calibration.
        weighted_energy = kr @ port_energy / (2 * radius)
        results[str(radius)] = dict(
            center=pick,
            strongest_signed_residual=r,
            separated_centers=int(sum(margins > 0)),
            radius_entropy=target[pick].tolist(),
            source_prediction=float(null_read[pick]),
            energy_weighted_sum=float(weighted_energy[pick]),
            witness_dyadic_numerator=yi.tolist(),
            witness_row_error_l1=str(Q(int(row_norm[pick]), radius * DYAD)),
        )
    return dict(
        radius_two_inverse_maximum_row_norm_upper=str(inverse_bound),
        density_maximum_norm_upper=str(density_norm),
        weights_denominator=DYAD,
        comparisons=results,
    )


def run():
    iv.dps = 45
    _, w, cells, _, _, _, ids, _ = energy.pair.setup()
    board, home, seams, terms, private, _ = energy.returned(w)
    # The same survey bypasses the three service gateways on every round trip.
    # Terminal material and source records refer to the parked chart.
    _, hops, _ = radar.chart(board, seams)
    distances = np.column_stack(
        (hops[:, home], *[hops[:, [u, v]].max(axis=1) for u, v in seams])
    )
    lookup = {r: (distances < r) @ 2 ** np.arange(7) for r in RADII}
    raw = json.loads((HERE / "closed_pairing.json").read_bytes())
    arms = {
        "cell_reweighted": energy.least_change(cells),
        "matched_path_energy": np.array(
            [Q(x) for x in raw["menus"]["energy"]["weights"]]
        )[ids],
    }
    rng = np.random.default_rng(10946223)
    n = len(w)
    ref = rng.multinomial(N, np.full(n, 1 / n))
    native = rng.multinomial(N, np.full(n, 1 / n))
    reports = {}
    for name, z in arms.items():
        held = rng.multinomial(N, np.asarray(z, float) / n)
        scores = []
        pred = []
        for mask in range(128):
            envelope, p = score(mask_ids(w, mask), ref, native, held, z)
            assert envelope[0] <= p <= envelope[1]
            scores.append(envelope)
            pred.append(p)
        assert abs(pred[-1]) < 1e-12
        reports[name] = dict(
            held_counts=held.tolist(),
            regional_bands=scores,
            regional_predictions=pred,
            survey=geometry_test(
                hops, lookup, scores, pred, energy_at_ports(board, terms, private, z)
            ),
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "CLOSED_PAIRING.md",
                HERE / "closed_pairing.json",
                HERE / "closed_energy.json",
            )
        },
        shots_per_bank=N,
        bit_order=["home0", "home1", "home2", "home3", "head0", "head1", "head2"],
        radius_convention="strict interior; w_R(d)=(R^2-d^2)_+/(2R)",
        count_halfwidth=EPS,
        reference_counts=ref.tolist(),
        native_counts=native.tolist(),
        regions={str(r): v.tolist() for r, v in lookup.items()},
        results=reports,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run()
    path = HERE / "closed_balls.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    for name, case in result["results"].items():
        print(
            name,
            {
                r: (v["separated_centers"], v["strongest_signed_residual"])
                for r, v in case["survey"]["comparisons"].items()
            },
        )
