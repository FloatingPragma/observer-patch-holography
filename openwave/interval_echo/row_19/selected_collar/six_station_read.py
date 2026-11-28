"""Read 30 native round-trip clocks across the counted source intervention."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from mpmath import iv
import numpy as np

from matched_einstein_limit.source_implication.full_geometry.station_balls.reciprocal_screen.range_contraction import (
    held_stations as survey,
)
from . import measure

HERE = measure.HERE
COHORTS = 16
SHOTS = COHORTS * survey.STARTS
SEED = 403714


def bounds(x):
    return [
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    ]


def shapes(square):
    geom = survey.geometry
    values = []
    for sites in survey.QUARTETS:
        i, *rest = sites
        gram = np.array([[geom.dot(square, i, j, k) for k in rest] for j in rest])
        values.append(
            (
                sum(geom.face_area(square, tri) for tri in combinations(sites, 3)),
                iv.sqrt(geom.curved.det_interval(gram)) / 6,
            )
        )
    gram = np.array(
        [[geom.dot(square, 0, j, k) for k in range(1, 6)] for j in range(1, 6)]
    )
    minors = [bounds(geom.curved.det_interval(gram[:k, :k])) for k in range(1, 6)]
    return values, minors


def decode(banks):
    # 60 means from independent launches on two endpoint preparations.
    # A timeout has score 3 and the earliest completed return has score 1/2.
    width = iv.mpf("2.5") * iv.sqrt(iv.ln(12000) / (2 * SHOTS)) + iv.mpf("1e-10")
    assert float(survey.history.timing_error()[1]) < 1e-12
    ends, ranks = {}, {}
    for label, bank in banks.items():
        means = np.full((6, 6), iv.mpf(0), dtype=object)
        for route, tape in zip(survey.ROUTES, bank, strict=True):
            assert len(tape) == 11 and min(tape) >= 0 and sum(tape) == SHOTS
            center = iv.mpf(sum(i * c for i, c in enumerate(tape, 2))) / (4 * SHOTS)
            means[route] = center + iv.mpf([-1, 1]) * width
        ends[label], ranks[label] = shapes(survey.geometry.observed_squared(means))
    changes = []
    for sites, lower, upper in zip(
        survey.QUARTETS, ends["minus"], ends["plus"], strict=True
    ):
        changes.append(
            dict(
                stations=list(sites),
                area=bounds(upper[0] - lower[0]),
                volume=bounds(upper[1] - lower[1]),
            )
        )
    return dict(changes=changes, leading_gram_minors=ranks)


def echo_probabilities(tape):
    # Invert Bernstein using var <= p and sqrt(2 a p) <= p/2 + a.
    # On the joint event p <= 2*p_hat + 4*a, where a = log(18000)/N.
    a = iv.ln(18000) / SHOTS
    estimates = {}
    for route, bins in zip(survey.ROUTES, tape, strict=True):
        assert len(bins) == 11 and sum(bins) == SHOTS and min(bins) >= 0
        center = iv.mpf(sum(bins[:-1])) / SHOTS
        width = iv.sqrt(2 * a * (2 * center + 4 * a)) + 2 * a / 3 + iv.mpf("1e-10")
        band = center + iv.mpf([-1, 1]) * width
        assert band.a > 0 and band.b < 1
        estimates[route] = band
    return estimates


def regional_weights(sites, probabilities):
    # Product law of a uniform ordered route and the reference echo flag,
    # compared with their measured joint law. All misses enter the score.
    selected = [r for r in survey.ROUTES if all(i in sites for i in r)]
    mean = sum(probabilities[r] for r in selected) / len(selected)
    return {
        r: iv.ln(mean * (1 - probabilities[r]) / (probabilities[r] * (1 - mean)))
        / len(selected)
        for r in selected
    }


def read_regional_source(banks, reference_tape):
    probabilities = echo_probabilities(reference_tape)
    source_prob = {key: echo_probabilities(bank) for key, bank in banks.items()}
    readings = {}
    for sites in survey.QUARTETS:
        weights = regional_weights(sites, probabilities)
        qd = sum(
            w * (source_prob["plus"][r] - source_prob["minus"][r])
            for r, w in weights.items()
        )
        readings[sites] = (weights, qd)
    return readings


def squared_range_jets(center, slopes):
    """Enroll all six clocks and reconstruct their squared symmetric ranges."""
    from matched_einstein_limit.source_implication.full_algebra import enclosure

    dims = len(slopes[0])
    zero = np.full((dims, dims), iv.mpf(0), dtype=object)
    logs = np.zeros((6, 6), dtype=object)
    for route, value, v in zip(survey.ROUTES, center, slopes, strict=True):
        logs[route] = enclosure.Jet(value, v, zero).log()
    gains = [sum(logs[i, j] - logs[j, i] for j in range(4)) * 0.25 for i in range(4)]
    gains += [
        sum(logs[i, j] - logs[j, i] + gains[j] for j in range(4)) * 0.25 for i in (4, 5)
    ]
    square = np.zeros((6, 6), dtype=object)
    for i, j in combinations(range(6), 2):
        square[i, j] = square[j, i] = (
            logs[i, j] + logs[j, i] - gains[i] - gains[j]
        ).exp()

    return square


def geometry_slopes(center, slopes):
    """Area and volume jets from interval-valued route means and contrasts."""
    square = squared_range_jets(center, slopes)

    def dot(i, j, k):
        return (square[i, j] + square[i, k] - square[j, k]) * 0.5

    for sites in survey.QUARTETS:
        area = sum(
            (square[i, j] * square[i, k] - dot(i, j, k) * dot(i, j, k)).sqrt() * 0.5
            for i, j, k in combinations(sites, 3)
        )
        i, *others = sites
        a, b, c = [square[i, m] for m in others]
        x, y, z = [dot(i, h, m) for h, m in combinations(others, 2)]
        det = a * b * c + x * y * z * 2 - a * z * z - b * y * y - c * x * x
        volume = det.sqrt() * (iv.mpf(1) / 6)
        yield sites, area.g, volume.g


def local_calibrations(bins, ref, displacement, banks, reference_tape):
    """Independently fit entropy/area and volume coefficients on each quartet."""
    from matched_einstein_limit.source_implication.full_algebra import enclosure
    from matched_einstein_limit.source_implication.shared_radar.measure import (
        preparations,
    )

    states = preparations()
    np.testing.assert_allclose(ref, states["equilibrium"], rtol=0, atol=1e-14)
    # A local unitary rotation preserves entropy exactly. Its generator has
    # eigenvalues ±1/2, so the half-difference of its ±pi/2 turns is its tangent.
    ordinary = states["ordinary"] - ref
    rotation = (states["plus"] - states["minus"]) / 2
    directions = [ordinary, rotation, displacement]
    regional_readings = read_regional_source(banks, reference_tape)
    response = np.einsum("rtbij,dbji->rd", bins, directions[:2]).real
    channel_slopes = dict(zip(survey.ROUTES, response, strict=True))
    deficits = np.einsum("rtbij,t->rbij", bins, np.arange(2, 12) / 4 - 3)
    center = 3 + np.einsum("rbij,bji->r", deficits, ref).real
    slopes = np.einsum("rbij,dbji->rd", deficits, directions).real
    assert max(sum(np.linalg.norm(h, ord="nuc") for h in d) for d in directions) <= 2
    # The established capped-mean effect error is below 1e-12. This wider
    # band also covers the binary64 4x4-block contractions and preparations.
    eps = iv.mpf(["-1e-9", "1e-9"])
    means = [iv.mpf(float(x)) + eps for x in center]
    tangents = [np.array([iv.mpf(float(a)) + eps for a in v]) for v in slopes]
    certpath = Path(enclosure.__file__).with_name("certificate.json")
    ordinary_band = json.loads(certpath.read_bytes())["signed_transfers"][
        "ordinary_entropy"
    ]
    entropy = iv.mpf(ordinary_band.strip("[]").split(","))
    q_interval = json.loads((HERE / "reference_read.json").read_bytes())["cases"][
        str(2**38)
    ]["readings"]["equilibrium"]["interval"]
    rows = []
    for sites, areas, volumes in geometry_slopes(means, tangents):
        ao, az, ad = areas
        vo, vz, vd = volumes
        denom = ao * vz - az * vo
        row = dict(stations=list(sites), calibration_determinant=bounds(denom))
        weights, qd = regional_readings[sites]
        qo, qz = [
            sum(
                w * (iv.mpf(float(channel_slopes[r][i])) + eps)
                for r, w in weights.items()
            )
            for i in (0, 1)
        ]
        row["regional_echo_charge_from_counts"] = bounds(qd)
        if not denom.a <= 0 <= denom.b:
            local_prediction = (
                (ad * vz - az * vd) * qo + (ao * vd - ad * vo) * qz
            ) / denom
            row["regional_echo_charge_from_geometry"] = bounds(local_prediction)
            row["regional_echo_disagreement"] = bool(
                local_prediction.b < qd.a or local_prediction.a > qd.b
            )
            predicted = entropy * (ad * vz - az * vd) / denom
            estimate = bounds(predicted)
            row.update(
                ordinary_area_gain=bounds(-entropy * vz / denom),
                predicted_charge=estimate,
                disagrees_with_source=estimate[1] < q_interval[0]
                or estimate[0] > q_interval[1],
            )
        rows.append(row)
    return rows


def run():
    ref, *_ = measure.jacobi.beam.aperture.source_and_changes()
    _, endpoints, *_ = measure.select_source(ref)
    ports, bins = survey.return_bins()
    rng = np.random.default_rng(SEED)
    banks = {}
    for label, state in zip(("minus", "plus"), endpoints, strict=True):
        counts = sum(
            np.array(survey.acquire(state, bins, rng), dtype=np.int64)
            for _ in range(COHORTS)
        )
        banks[label] = counts.tolist()
    reference_tape = sum(
        np.array(survey.acquire(ref, bins, rng), dtype=np.int64) for _ in range(COHORTS)
    ).tolist()
    return dict(
        source_recipe_sha256=sha256((HERE / "results.json").read_bytes()).hexdigest(),
        ports=list(ports),
        routes=list(map(list, survey.ROUTES)),
        starts_per_route=SHOTS,
        seed=SEED,
        counts=banks,
        regional_echo_reference=reference_tape,
        local_calibrations=local_calibrations(
            bins, ref, endpoints[1] - endpoints[0], banks, reference_tape
        ),
        **decode(banks),
    )


def main():
    opts = argparse.ArgumentParser(description=__doc__)
    opts.add_argument("--write", action="store_true")
    opts.add_argument("--counts-only", action="store_true")
    args = opts.parse_args()
    outfile = HERE / "six_station_read.json"
    iv.dps = 35
    if args.counts_only:
        record = json.loads(outfile.read_bytes())
        result = decode(record["counts"])
        assert all(result[k] == record[k] for k in result)
        readings = read_regional_source(
            record["counts"], record["regional_echo_reference"]
        )
        for row in record["local_calibrations"]:
            assert (
                bounds(readings[tuple(row["stations"])][1])
                == row["regional_echo_charge_from_counts"]
            )
    else:
        result = run()
        if args.write:
            outfile.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            assert result == json.loads(outfile.read_bytes())
    for key in ("area", "volume"):
        rows = [r[key] for r in result["changes"]]
        print(
            key,
            dict(
                negative=sum(b < 0 for a, b in rows),
                positive=sum(a > 0 for a, b in rows),
            ),
            "original quartet",
            rows[0],
        )


if __name__ == "__main__":
    main()
