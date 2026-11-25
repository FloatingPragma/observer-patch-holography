"""Transport every traceless source observable through each material instrument."""

import argparse
import json
from fractions import Fraction
from math import ceil, floor

import numpy as np
from mpmath import iv, mp
from scipy.linalg import qr

from . import mounts

cuts, full = mounts.cuts, mounts.full_source


def midpoint(x):
    ratio = cuts.survey.dyadic(x.mid._mpi_[0])
    return mp.mpf(ratio.numerator) / ratio.denominator


def center(matrix):
    return mp.matrix([[midpoint(x) for x in row] for row in matrix])


def publish(values, digits=55):
    return np.array(
        [
            [mp.nstr(values[i, j], digits) for j in range(values.cols)]
            for i in range(values.rows)
        ]
    )


def intervals(decimals):
    return np.vectorize(iv.mpf, otypes=[object])(decimals)


def endpoint(x, index):
    return cuts.survey.dyadic(x._mpi_[index])


def upper(x, numbers=3):
    """Scientific notation, rounded upwards with integer arithmetic."""
    val = endpoint(x, 1)
    if val == 0:
        return "0"
    exponent = floor(np.log10(float(val))) - numbers + 1
    place = Fraction(10) ** exponent
    return f"{ceil(val / place)}e{exponent}"


def interval_rows(rows, scale):
    b = np.array(
        [[iv.mpf([lo, hi]) / scale for lo, hi in row["row"]] for row in rows],
        dtype=object,
    )
    b[:, -1] = b[:, -1] / iv.sqrt(8)
    return b


def read_native():
    data = json.loads((mounts.HERE / "effects.json").read_text())
    rows = []
    for arm, label in zip(data["charts"], mounts.MOUNTS, strict=True):
        if tuple(arm[k] for k in ("center", "target", "orientation")) != label:
            raise ValueError(
                "chart identity disagrees with the enrolled material experiment"
            )
        keys = [(c["charge"], c["reading"]) for c in arm["rows"]]
        if len(set(keys)) != 31:
            raise ValueError("each chart needs 31 distinct enrolled physical channels")
        b = interval_rows(arm["rows"], data["scale"])
        rows.append(b)
    return rows


def archive_operator():
    # Assemble the material-and-marginal reader exactly as specified on the reference chart.
    cuts_original = json.loads((cuts.HERE / "effects.json").read_text())["cuts"]
    by = {tuple(r["channel"]): r["row"] for r in cuts_original}
    region = json.loads((cuts.survey.HERE / "translation.json").read_text())["rows"]
    joint = cuts.interval_rows(region + [by[k] for k in cuts.CHANNELS])
    local = cuts.interval_rows(
        json.loads((cuts.HERE / "local_effects.json").read_text())["rows"]
    )
    old = json.loads((cuts.HERE / "full_source.json").read_text())
    w = intervals(old["material_weights"])
    m = intervals(old["marginal_and_parity_weights"])
    return np.r_[w @ joint, w @ local + m]


POOLS = ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2))


def instruments(enrolled):
    data = []
    for pool in POOLS:
        keys = [(a, i) for a in pool for i in range(31)]
        if len(pool) > 1:
            # This choice depends only on native material response, in raw units.
            table = cuts.mid(np.vstack([enrolled[a] for a in pool]))
            pivots = qr(table.T, mode="economic", pivoting=True)[2][:31]
            keys = [keys[i] for i in sorted(pivots)]
        data.append((keys, np.array([enrolled[a][i] for a, i in keys])))
    return data


def certify():
    iv.dps = 75
    mp.dps = 80
    enrolled = read_native()
    old = archive_operator()
    # All charts take the same centre of the identically specified old observable.
    # The actual modular generator and geometry are consulted only in the checks.
    target = center(old[None, :])
    expected = np.r_[cuts.exact_entropy_center(), full.local_log_center()]
    loc = json.loads((cuts.survey.HERE.parent / "certificate.json").read_text())
    radius = cuts.ri(str(loc["state_HS_error_upper"]))
    log_radius = radius / (cuts.ri(".02") - radius)
    reports = []
    midpoints = []
    encoders = []
    defects = []
    targets = []
    for keys, b in instruments(enrolled):
        mat = center(b)
        inv = mat**-1
        quantized = publish(inv)
        R = intervals(quantized)
        defect = cuts.norm(np.eye(31, dtype=int) - R @ b)
        if defect.b >= cuts.ri("1e-19"):
            raise ValueError(f"material-only source inverse is uncertified: {defect}")
        sigma = (1 - defect) / cuts.norm(R)
        # Asymmetric SVD on floating-point data can lose singular values of the
        # raw native matrix. A column of the verified inverse gives an upper bound.
        largest = max(cuts.norm(col).a for col in R.T)
        sigma_upper = (1 + defect) / largest
        prices = publish(target * inv)[0]
        w = intervals(prices)
        represented = w @ b
        discrepancy = cuts.norm(represented - old)
        total = cuts.norm(represented - expected) + log_radius
        assert total.b < cuts.ri("3e-5")
        midpoints.append(intervals(publish(mat)))
        encoders.append(R)
        defects.append(defect)
        targets.append(represented)
        reports.append(
            dict(
                chart_indices=POOLS[len(reports)],
                physical_rows=keys,
                rank=31,
                left_inverse_defect_upper=upper(defect),
                minimum_singular_lower_in_1e30=floor(endpoint(sigma, 0) * 10**30),
                minimum_singular_upper_in_1e30=ceil(endpoint(sigma_upper, 1) * 10**30),
                sensor_error_l2_to_source_hs_upper=upper(cuts.norm(R) / (1 - defect)),
                prices=prices.tolist(),
                price_l1_upper=upper(sum(abs(x) for x in w)),
                old_reader_discrepancy_upper=upper(discrepancy),
                true_modular_error_upper=upper(total),
            )
        )
    inter = []
    for i, pool in enumerate(POOLS):
        for j in range(len(enrolled)):
            if j in pool:
                continue
            outgoing = enrolled[j]
            M = midpoints[j]
            transport = M @ encoders[i]
            # R and M are published decimals; this includes their rounding and the
            # native error on both instruments before comparing arbitrary sources.
            bias = cuts.norm(outgoing - M) + cuts.norm(M) * defects[i]
            physical = cuts.norm(targets[j] - targets[i])
            assert bias.b < cuts.ri("1e-20") and physical.b < cuts.ri("1e-19")
            inter.append(
                dict(
                    origin=i,
                    destination=j,
                    detector_vector_error_per_unit_HS_upper=upper(bias),
                    detector_error_l2_amplification_upper=upper(cuts.norm(transport)),
                    modular_disagreement_per_unit_HS_upper=upper(physical),
                )
            )
    return dict(mounts=mounts.MOUNTS, instruments=reports, transfers=inter)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    report = certify()
    if args.write:
        (mounts.HERE / "transfer.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))
