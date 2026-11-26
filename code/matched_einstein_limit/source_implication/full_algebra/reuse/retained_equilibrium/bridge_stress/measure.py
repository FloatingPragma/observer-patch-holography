"""Check independent material readings on the retained constant-entropy and -volume curve."""

import json
from pathlib import Path

from mpmath import iv
from mpmath.libmp import round_ceiling, round_floor, to_int


HERE = Path(__file__).parent
SOURCE = HERE.parents[3]
MATERIAL = SOURCE / "directional_stress"
ARMS = ("equilibrium", "ordinary", "pairing", "plus", "minus")


def read(path):
    return json.loads(path.read_text(), parse_float=str)


def band(text):
    return iv.mpf(text.strip("[]").split(","))


def outward(x):
    """Strictly outward binary64-independent grid in millionths."""
    n = 10**6
    lower, upper = (x * n)._mpi_
    return [to_int(lower, round_floor), to_int(upper, round_ceiling)]


def score_counts(current_data, work_data, q, ratio, rotation):
    # Two banks, twenty raw preparations; Bonferroni gives 99% jointly.
    # We allocate 0.005 failures to each bank and form all contrasts on
    # that same simultaneous event, including any r in the certified band.
    nc, nw = current_data["trials_per_setting"], work_data["trials_per_preparation"]
    mag = current_data["binary_readout_magnitude"]
    phases = [
        sum(k * iv.sin(2 * iv.pi * j * k / 19) for k in range(1, 10)) / 19
        for j in range(19)
    ]
    curr, work = {}, {}
    for arm in ARMS:
        counts = current_data["charges"][q]["phase_settings"][arm]["plus_counts"]
        assert len(counts) == 19 and all(0 <= h <= nc for h in counts)
        curr[arm] = mag * sum(
            w * (2 * iv.mpf(h) / nc - 1) for w, h in zip(phases, counts, strict=True)
        )
        hits = work_data["charges"][q]["outcome_counts"][arm]
        assert len(hits) == 3 and min(hits) >= 0 and sum(hits) == nw
        work[arm] = iv.mpf(8 * hits[1] + 9 * hits[2]) / nw
    phase_error = mag * iv.sqrt(2 * iv.log(2 * 190 / iv.mpf(".005")) / nc)
    work_error = 9 * iv.sqrt(iv.log(2 * 10 / iv.mpf(".005")) / (2 * nw))
    # The modulus sum of the five preparation coefficients is 2+|rotation|
    # since 0<ratio<1. Form the point estimate jointly in r to avoid counting
    # the same large phase offsets several times.
    oangle, cangle = angles()

    def slope(record, name, angle):
        return (
            record[name]
            - record["equilibrium"]
            + angle * (record["plus"] - record["minus"]) / 2
        )

    def measured(record, error):
        estimate = slope(record, "pairing", cangle) - ratio * slope(
            record, "ordinary", oangle
        )
        return estimate + iv.mpf([-1, 1]) * (2 + abs(rotation)) * error

    return measured(work, work_error), measured(
        curr, sum(map(abs, phases)) * phase_error
    )


def angles():
    data = read(SOURCE / "signed_transfers.json")["tests"]
    return tuple(
        iv.mpf(next(row["angle_per_unit"] for row in data if row["source"] == name))
        for name in ("ordinary", "collar")
    )


def report():
    iv.dps = 50
    base = read(SOURCE / "full_algebra/certificate.json")
    src = base["signed_transfers"]
    ordinary, collar = map(band, (src["ordinary_entropy"], src["pairing_entropy"]))
    ratio = collar / ordinary
    assert 0 < ratio.a and ratio.b < 1 and band(src["right_adjuster_volume"]).b < 0
    delta_energy = -ratio * band(src["ordinary_energy"])
    eta = band(base["calibrated_gain"])
    gap = iv.mpf(base["global_survey"]["feasible_objective_gap_coefficient"])
    area_gap = gap * (2 * delta_energy) ** 2 / eta
    assert area_gap.a > iv.mpf("0.00068")
    jo = read(MATERIAL / "certificate.json")
    wo = read(MATERIAL / "traction/certificate.json")
    cm = read(MATERIAL / "phase_read.json")
    wm = read(MATERIAL / "traction/strain.json")
    a, c = angles()
    rotation = c - ratio * a
    rows = {}
    for mode in ("source_certificate", "simulated_counts"):
        mat = []
        for q in ("-1", "1"):
            if mode == "source_certificate":

                def matched(record):
                    return band(record["pairing"]) - ratio * band(record["ordinary"])

                work = matched(wo["charges"][q]["observables"]["full_seam"])
                current = matched(jo[q])
            else:
                work, current = score_counts(cm, wm, q, ratio, rotation)
            assert min(work.a, current.a) > 0
            mat.append((work, current))
        (wmn, jmn), (wpl, jpl) = mat
        det = wmn * jpl - wpl * jmn
        assert det.a > 0
        rows[mode] = {
            key: {
                "work": outward(w),
                "current": outward(j),
                "work_over_current": outward(w / j),
            }
            for key, (w, j) in zip(("vacancy", "excess"), mat, strict=True)
        }
        rows[mode]["determinant"] = outward(det)
    # Hold the vacancy current and survey volume without using entropy or area.
    # The two tight current ratios share the trace-error and volume-adjuster
    # coefficients inside each original source certificate.
    labels = {q: band(jo[q]["pairing_over_ordinary"]) for q in ("-1", "1")}
    held_ratio = labels["-1"]
    current_plus = band(jo["1"]["ordinary"]) * (labels["1"] - held_ratio)
    entropy = collar - held_ratio * ordinary
    current_count = {}
    for q in ("-1", "1"):
        _, current = score_counts(cm, wm, q, held_ratio, c - held_ratio * a)
        current_count[q] = outward(current)
    assert current_plus.b < 0 and entropy.a > 0
    # Calibrate again using ONLY the vacancy count bank. On the single
    # simultaneous event the true contrast is affine in r, so endpoint
    # signs bracket its root. Extrema of the held excess contrast over this
    # bracket occur at those same endpoints.
    bracket = (iv.mpf(".2805"), iv.mpf(".283"))
    endpoint = {
        q: [score_counts(cm, wm, q, r, c - r * a)[1] for r in bracket]
        for q in ("-1", "1")
    }
    assert endpoint["-1"][0].b < 0 < endpoint["-1"][1].a
    held_counts = [
        min(outward(x)[0] for x in endpoint["1"]),
        max(outward(x)[1] for x in endpoint["1"]),
    ]
    assert held_counts[1] < 0
    rows["count_only_material_calibration"] = {
        "vacancy_pairing_over_ordinary": [280500, 283000],
        "vacancy_current_at_ratio_endpoints": [outward(x) for x in endpoint["-1"]],
        "held_excess_current": held_counts,
    }
    rows["fixed_vacancy_current_and_volume"] = {
        "vacancy_current": [0, 0],
        "excess_current": outward(current_plus),
        "entropy": outward(entropy),
        "area": outward(-entropy / eta),
        "nominal_current_counts": current_count,
    }
    rows.update(
        units="outward endpoints in millionths",
        entropy_ratio=outward(ratio),
        parked_energy_derivative=outward(delta_energy),
        liminf_area_loss_per_step_squared_lower=outward(area_gap)[0],
        count_coverage=0.99,
    )
    return rows


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
