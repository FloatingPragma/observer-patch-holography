"""Transfer an ordinary source calibration to a collar on two independent probe species."""

import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path

import numpy as np
from matched_einstein_limit.source_implication.beam_response import measure as camera

HERE = Path(__file__).resolve().parent
CHANGES = ("ordinary", "unpairing")
TIMES = (3.75, 4.0, 4.25)


def common_gain_bound():
    cert = json.loads(camera.HERE.joinpath("certificate.json").read_text())
    def interval(q):
        return [Decimal(v.strip()) for v in cert[q]["ordinary_calibrated_ratio"][1:-1].split(",")]
    with localcontext() as context:
        context.prec = 90
        lo_v, hi_v = interval("-1")
        lo_e, hi_e = interval("1")
        assert lo_e < hi_e < 0 < lo_v < hi_v < 1
        return dict(vacancy=[str(lo_v), str(hi_v)], excess=[str(lo_e), str(hi_e)],
                    lower_bound_on_max_abs_error_any_common_gain=str((lo_v - hi_e) / 2),
                    lower_bound_on_excess_error_at_unit_gain=str(1 - hi_e))


def endpoint_step(ref, tangent):
    step = 1/256
    for _ in range(50):
        if all(np.min(np.linalg.eigvalsh(ref + sign * step * tangent[key])) > 0
               for key in CHANGES for sign in (-1, 1)):
            return step
        step /= 2
    raise ValueError("no resolvable two-sided source change")


def exposure(table, ref, tangent, charges, step):
    p = camera.probabilities(table, ref)
    def record(pixels):
        law = np.column_stack((pixels.T, 1 - pixels.sum(axis=0)))
        assert np.min(law) >= 0
        np.testing.assert_allclose(law.sum(axis=1), 1, atol=2e-14, rtol=0)
        return law.tolist()
    totals = {k: charges[k]["full"] for k in CHANGES}
    slopes = {k: camera.area_derivative(p, camera.probabilities(table, tangent[k]))
              for k in CHANGES}
    def ratio(slopes):
        return slopes["unpairing"] / (slopes["ordinary"] * totals["unpairing"] / totals["ordinary"])
    endpoints = {}
    for h in (step, step/2, step/16, step/32):
        finite = {}
        for k in CHANGES:
            laws = []
            areas = []
            for sign in (-1, 1):
                state = ref + sign * h * tangent[k]
                np.testing.assert_allclose(state, state.conj().transpose(0, 2, 1), atol=1e-12, rtol=0)
                np.testing.assert_allclose(np.trace(state, axis1=1, axis2=2).sum(), 1, atol=1e-12)
                assert np.min(np.linalg.eigvalsh(state)) > 0
                pixels = camera.probabilities(table, state)
                laws.append(record(pixels))
                areas.append(np.linalg.det(camera.image(pixels)[0]))
            finite[k] = dict(minus_plus_laws=laws, slope=(areas[1] - areas[0]) / (2*h))
        finite["ratio"] = ratio({k: finite[k]["slope"] for k in CHANGES})
        endpoints[str(h)] = finite
    # Cancel the leading centered-difference error, using the two physical step sizes.
    richardson = (4 * endpoints[str(step/2)]["ratio"] - endpoints[str(step)]["ratio"]) / 3
    fine = (4 * endpoints[str(step/32)]["ratio"] - endpoints[str(step/16)]["ratio"]) / 3
    np.testing.assert_allclose(fine, ratio(slopes), rtol=2e-5, atol=1e-7)
    return dict(reference_law=record(p), area_slopes=slopes,
                collar_area_slope_from_ordinary=slopes["ordinary"] * totals["unpairing"] / totals["ordinary"],
                ratio=ratio(slopes), finite_changes=endpoints, extrapolated_ratio=richardson, finer_extrapolated_ratio=fine)


def run():
    board = camera.front.optics.pairing.formation.federation(1)
    route = camera.front.optics.pairing.formation.region(board)
    rings = camera.screens(board, route)
    ref, tangent, charges, _ = camera.aperture.source_and_changes()
    step = endpoint_step(ref, tangent)
    observations = {}
    for q in (-1, 1):
        observations[str(q)] = {str(t): exposure(camera.effects(board, route, rings, q, t),
                                                 ref, tangent, charges, step) for t in TIMES}
    bound = common_gain_bound()
    for q, key in ((-1, "vacancy"), (1, "excess")):
        lo, hi = map(float, bound[key])
        assert lo <= observations[str(q)]["4.0"]["ratio"] <= hi
    return dict(cells=board.cells, launch_ports=rings[0], receiver_ports=rings[1],
                last_outcome="outside_screen", source_slopes=charges, physical_step=step,
                exact_source_at_four=bound, probes=observations)


def compare(actual, stored):
    if isinstance(actual, dict):
        assert actual.keys() == stored.keys()
        for key in actual:
            compare(actual[key], stored[key])
    elif isinstance(actual, (int, float, list)):
        np.testing.assert_allclose(actual, stored, rtol=2e-6, atol=2e-12)
    else:
        assert actual == stored


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(path.read_text())
        assert result["exact_source_at_four"] == saved["exact_source_at_four"]
        for q in result["probes"]:
            compare(result["probes"][q], saved["probes"][q])
        for key in ("cells", "launch_ports", "receiver_ports", "last_outcome", "source_slopes", "physical_step"):
            compare(result[key], saved[key])
    print({q: [row["ratio"] for row in by_time.values()] for q, by_time in result["probes"].items()})
