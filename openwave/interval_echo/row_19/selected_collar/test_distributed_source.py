"""Independent checks of the spatial source assignment and the ordinary controls."""

import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np

from openwave.interval_echo.row_19.selected_collar import observation_density as density

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / (name + ".json")).read_bytes())


def test_address_enrollment_preserves_the_flight():
    _, sources, *_ = density.c.trials()
    x = sources["reference"]

    def run(ports):
        return density.history(
            x,
            (0,),
            np.random.default_rng(42),
            signs=(-1,),
            read_ports=ports,
            windows=(1,),
        )[2][1]

    separate, previous = run((0, 1, 38, 39, 5, 6)), run(None)
    for key, records in separate.items():
        np.testing.assert_allclose(records[0], previous[key][0], rtol=0, atol=3e-14)
        np.testing.assert_allclose(
            records[1:].sum(axis=0), previous[key][1], rtol=0, atol=3e-14
        )


def test_every_distributed_exclusion_and_source_prediction():
    tape = read("distributed_source")
    for arm in [
        *tape["target_counts"].values(),
        *(v for z in tape["standard_counts"].values() for v in z.values()),
    ]:
        for head in arm:
            for counts in head:
                assert len(counts) == 7 and all(len(site) == 8 for site in counts)
                assert sum(map(sum, counts)) == tape["starts"]
    d = tape["result"]
    for source in d["sources"].values():
        for direct, comparisons in [
            (False, source["comparisons"]),
            (True, source["direct_comparisons"]),
        ]:
            target = dict(source["target"])
            if direct:
                target.update(source["direct_material"])
            for group, c in comparisons.items():
                if c["status"] == "excluded":
                    lhs = dict.fromkeys(d["columns"], F(0))
                    rhs = F(0)
                    for name, side, num, den in c["coefficients"]:
                        assert name in d["groups"][group] and 0 < num and 0 < den
                        j = int(side == "upper")
                        w = F(num, den) * (-1 if j else 1)
                        rhs += w * F(target[name][1 - j])
                        for hist, z in d["columns"].items():
                            lhs[hist] += w * F(z[name][j])
                    assert min(lhs.values()) >= 0 > rhs
                elif c["status"] == "compatible_relaxation":
                    w = {k: F(n, den) for k, n, den in c["allocation"]}
                    assert all(v >= 0 for v in w.values())
                    for name in d["groups"][group]:
                        lo, hi = [
                            sum(v * F(d["columns"][k][name][j]) for k, v in w.items())
                            for j in range(2)
                        ]
                        assert lo <= F(target[name][1]) and hi >= F(target[name][0])
                else:
                    assert c["status"] == "unresolved"
        mass_lo, mass_hi = map(F, source["inferred_total_energy"])
        for n, z in source["source_only_prediction"].items():
            lo, hi = map(F, z["prediction"])
            for k, v in d["columns"].items():
                elo, ehi = map(F, d["ordinary_total_energy"][k])
                endpoints = [
                    F(a) * m / e
                    for a in v[n]
                    for m in (mass_lo, mass_hi)
                    for e in (elo, ehi)
                ]
                # Archived inputs have one extra outward binary64 ulp.
                assert float(min(endpoints)) >= float(lo) - 1e-14 * max(
                    abs(float(lo)), 1e-30
                )
                assert float(max(endpoints)) <= float(hi) + 1e-14 * max(
                    abs(float(hi)), 1e-30
                )
            target = source["direct_material"].get(n, source["target"][n])
            assert z["excludes_programmed_read"] == (lo > target[1] or hi < target[0])


def test_exact_source_match_and_independent_clock_records():
    d = read("regional_source_calibration")
    experiment = read("matched_regional_flight")
    full = d["calibrations"]["energy+connected+marginal"]
    assert len(full["basis"]) == 21 and len(full["ordinary"]) == 15
    assert full["weighted_inverse_residual_upper"] < 1
    assert all(r["excluded"] for r in full["ordinary"].values())
    assert all(
        0 <= v <= 1
        for p in experiment["ordinary_probabilities"].values()
        for v in p.values()
    )
    assert all(
        sum(c) == experiment["starts"]
        for arm in experiment["clock_counts"].values()
        for c in arm
    )
    for n, z in experiment["shape"].items():
        obs = z["center_clock_tangent"]
        for pred in (
            "center_coincidence_forecast",
            "center_forecast_with_setting_error",
        ):
            lo, hi = z[pred]
            assert max(lo, obs[0]) <= min(hi, obs[1]), (n, pred, z)
    assert any(a * b > 0 for a, b in experiment["summed_pauli_contrasts"].values())


def test_missing_source_rank_certificates_as_exact_rationals():
    calibration = read("regional_source_calibration")["calibrations"]
    groups_by_source = read("regional_rank")["results"]
    transport = read("regional_transport")
    calibration["energy transport"] = transport["calibration"]
    groups_by_source["energy transport"] = transport["residual_rank"]
    for fields, groups in groups_by_source.items():
        residuals = calibration[fields]["ordinary"]
        for proof in groups.values():
            rows, cols = proof["meters"], proof["histories"]
            k = proof["rank_lower"]
            assert k == len(set(rows)) == len(set(cols)) > 0
            v = list(map(F, proof["norm_weights"]))
            assert min(v) > 0
            K = [list(map(F, row)) for row in proof["preconditioner"]]
            for i in range(k):
                bound = F(0)
                for j, h in enumerate(cols):
                    lower, upper = F(i == j), F(i == j)
                    for a, meter in zip(K[i], rows, strict=True):
                        lo, hi = sorted(
                            a * F(x)
                            for x in residuals[h]["observed_minus_prediction"][meter]
                        )
                        lower -= hi
                        upper -= lo
                    # The component interval of I-KR includes the true residual.
                    bound += max(abs(lower), abs(upper)) * v[j]
                assert bound <= F(proof["contraction_upper"]) * v[i] < v[i]
