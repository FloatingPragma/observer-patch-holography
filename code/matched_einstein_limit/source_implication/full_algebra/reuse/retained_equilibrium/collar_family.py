"""Compose the certified source, retained-feedback and signed-volume bounds."""

import json
from fractions import Fraction as F
from pathlib import Path

from .control import spectrum_checked
from .probe import density_center


HERE = Path(__file__).parent


def interval(text):
    lower, upper = map(F, text.strip("[]").split(","))
    assert lower <= upper
    return lower, upper


def verify():
    previous = json.loads((HERE.parents[1] / "certificate.json").read_text(), parse_float=F)
    feedback = json.loads((HERE / "results.json").read_text(), parse_float=F)
    center, odd_weight, _, radius = density_center()
    t = feedback["twirl"]
    error = feedback["certified_absolute_error_per_partial_sum_slack"]
    assert t == F(1, 8) and radius == previous["state_HS_error_upper"]
    # ||X* - center||_1 <= sqrt(8) radius on the complete block algebra.
    trace_error = F("0.000000566")
    assert 8 * radius**2 < trace_error**2
    eigen_error = F("0.00000004")  # certified by spectrum_checked
    assert error <= F("0.000013")
    margins = []
    for block, p, record in zip(center, (1 - odd_weight, odd_weight), feedback["blocks"]):
        spectrum = [F.from_float(float(v)) for v in spectrum_checked(block)]
        assert spectrum[-1] - eigen_error - trace_error > F("0.019")
        for k, f in enumerate(record["with_bit"], start=1):
            # s_k(X*) - k tr(X*)/4. Both the Ky Fan sum and sector trace
            # can move by at most trace_error; p is the exact center weight.
            spread_upper = (
                sum(spectrum[:k]) + k * eigen_error + trace_error
                - F(k, 4) * (p - trace_error)
            )
            margins.append((1 - t) * (f - error) - t * spread_upper)
    assert min(margins) > F("0.0020454")
    # Varying the sector trace adds at most t*k/4 to the Lipschitz coefficient.
    assert min(margins) - 2 * F("0.001") > 0

    signed = previous["signed_transfers"]
    volume = interval(signed["right_adjuster_volume"])
    assert volume[1] < F("-0.00002564")
    gain = interval(previous["calibrated_gain"])
    assert 0 < gain[0]
    area = {}
    for kind in ("ordinary", "pairing"):
        lo, hi = interval(signed[f"{kind}_entropy"])
        assert lo > 0
        area[kind] = (-hi / gain[0], -lo / gain[1])
    gap = previous["global_survey"]["feasible_objective_gap_coefficient"]
    # Near the source lambda_min > .019 > exp(-4), so entropy is locally
    # 4-Lipschitz on normalized states. The collar slope is strictly > .23.
    assert interval(signed["pairing_entropy"])[0] > F(".23")
    assert gap * (F(".23") / 4) ** 2 > F(".00046")
    return margins, signed, area


if __name__ == "__main__":
    lower, _, area_slopes = verify()
    print("Post-twirl margins, lower bounds:", [float(x) for x in lower])
    for name, bounds in area_slopes.items():
        print(name, "signed area derivative:", tuple(map(float, bounds)))
