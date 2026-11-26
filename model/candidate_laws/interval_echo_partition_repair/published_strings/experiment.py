"""Three same-source small-ball inputs with no feedback into Interval Echo."""
import argparse
from fractions import Fraction
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np

from ..proper_balls.measure import AXES, FRAMES, lift
from . import finite_rods
from . import geometry as geo
from matched_einstein_limit.frame import decode
from .archive import source_densities
from .memories import Pair

HERE = Path(__file__).resolve().parent
SCALES = (32, 128, 512)


@lru_cache(maxsize=128)
def _ideal_contrasts(K, ell, stretch):
    """The same product-chart quadratures with the flat constant integrated exactly."""
    if not (0 <= K <= 2*geo.KAPPA and 0 < ell <= 1/32 and 1 <= stretch <= 9):
        raise ValueError("ideal survey outside the enrolled curvature, radius or frame range")
    if K == 0:
        return geo.Ball(ell, 4*math.pi*ell**3/3, 4*math.pi*ell**2), 0., 0.
    length, transverse, angular, nodes, weights = geo.quadrature(stretch, 16)
    square = K*length**2
    # sinc(sqrt(z))²-1. On the enrolled domain z<0.13 and the first omitted
    # term is below 3e-21, with geometrically smaller subsequent terms.
    coefficients = [(-1)**j*2**(2*j+1)/math.factorial(2*j+2) for j in range(1,9)]

    def excess(radii):
        z = np.asarray(radii)[:, None, None]**2 * square
        x = transverse*z*np.polynomial.polynomial.polyval(z,coefficients)
        return np.sum(angular*x/(1+np.sqrt(1+x)), axis=(1,2))

    def volume_excess(r):
        return float(weights @ (nodes**2 * excess(r*nodes)))

    shift = 0.
    for _ in range(24):
        v = volume_excess(ell*(1+shift))
        update = math.expm1(-math.log1p(3*v)/3)
        if abs(update-shift) <= 2*math.ulp(update):
            shift = update
            break
        shift = update
    else:
        raise ArithmeticError("excess-volume radius failed to converge")
    r = ell*(1+shift)
    t2_minus_one = shift*(2+shift)
    t4_minus_one = t2_minus_one*(2+t2_minus_one)
    boundary = float(excess([r])[0])
    area = 4*math.pi*(t2_minus_one+(1+t2_minus_one)*boundary)/ell**2
    # The archived integration by parts uses r * sum(w*u*V(r*u)).
    inner = np.array([volume_excess(r*u) for u in nodes])
    weighted = 4*math.pi*(t4_minus_one/15 + (1+t4_minus_one)*float(weights @ (nodes**4*inner)))
    ball = geo.Ball(r, 4*math.pi*r**3*(1/3+volume_excess(r)), 4*math.pi*r*r*(1+boundary))
    return ball, area, weighted


def _frames():
    out = {}
    for label, u in zip(FRAMES.values(), geo.FRAME, strict=True):
        out.setdefault(label, tuple(map(float, u)))
    if out.keys() != AXES.keys():
        raise ArithmeticError("survey and theorem frame orders disagree")
    return out


def endpoint_positions(shape, ell):
    """Two material worldlines on a fixed ray, placed without any source reading."""
    direction = np.array((1., 2., 4.))/math.sqrt(21)
    crossing = ell/np.linalg.norm(np.linalg.solve(lift(shape)[1], direction))
    return [(crossing*(1 + sign/160)*direction).tolist() for sign in (-1, 1)]


def _trials(densities):
    yield from ((name, rho, collar, 1, True) for name, (rho, collar) in densities.items())
    rho, collar = densities["held"]
    yield "half_defect", rho, collar, 0.5, True
    yield "blind_pulse", rho, collar, 1, False


def _board(densities, *, finite):
    rows = []
    for denominator in SCALES:
        ell = 1/denominator
        for name, u in _frames().items():
            trials = {}
            for key, rho, collar, strength, pulse_reads_collar in _trials(densities):
                K = geo.curvature(rho, collar, collar_strength=strength)
                if finite:
                    # The downstream lottery has exact rational frequencies.
                    # Recover them from the declared 33-address source.
                    gain = finite_rods.angle_gain(Fraction(1, denominator))
                    density = Fraction(rho).limit_denominator(33) + Fraction(strength)*Fraction(collar).limit_denominator(33)
                    record = finite_rods.read_ball(gain*density, Fraction(1, denominator), u)
                    ball = geo.Ball(float(record.radius), 4*math.pi*float(record.volume_over_4pi),
                                    4*math.pi*float(record.area_over_4pi))
                    weighted_volume = 4*math.pi*float(record.weighted_over_4pi)
                    flat = finite_rods.read_ball(0, Fraction(1, denominator), u)
                    area = 4*math.pi*float((record.area_over_4pi-flat.area_over_4pi)*denominator**4)
                    weighted = 4*math.pi*float(record.weighted_over_4pi*denominator**4-Fraction(1,15))
                else:
                    ball, area, weighted = _ideal_contrasts(K, ell, float(geo.spatial_stretch(u)))
                    weighted_volume = ell**4*(4*math.pi/15+weighted)
                stress = geo.project(geo.energy_tensor(rho, collar), u)
                pulse = geo.project(geo.energy_tensor(rho, collar if pulse_reads_collar else 0), u)
                pair = Pair.from_weighted_work(pulse*weighted_volume)
                trials[key] = dict(proper_radius=ball.radius, proper_volume=ball.volume, proper_area=ball.area,
                                   weighted_volume=weighted_volume, angle=pair.angle_increment,
                                   normalized_area_contrast=area, weighted_excess=weighted, pulse_stress=pulse,
                                   action_stress=stress, metric_charge=geo.project(geo.einstein_tensor(K), u))
            rows.append(dict(shape=name, inverse_ell=denominator, trials=trials,
                             worldlines=endpoint_positions(name, ell)))
    endpoints = np.array([r["worldlines"] for r in rows])
    for i, row in enumerate(rows):
        distance = np.linalg.norm(endpoints @ np.linalg.inv(lift(row["shape"])[1]), axis=2)
        ell = 1/row["inverse_ell"]
        reference = distance < ell
        whole = np.flatnonzero(reference.all(axis=1)).tolist()
        row["whole_pairs"] = whole
        for key, trial in row["trials"].items():
            contained = distance < trial["proper_radius"]
            if (not np.array_equal(reference, contained)
                    or np.flatnonzero(contained[:, 0] != contained[:, 1]).tolist() != [i]):
                raise ArithmeticError("a source changed membership of a registered endpoint")
            pair = Pair(trial["angle"])
            interior_kl = sum(Pair(rows[j]["trials"][key]["angle"]).whole_modular_change for j in whole)
            entropy = pair.cut_entropy_change
            modular = pair.cut_modular_change + interior_kl
            area = trial.pop("normalized_area_contrast")
            pulse = trial.pop("pulse_stress")
            # Subtract the flat weighted work before adding the nonlinear cut response.
            bulk = math.fsum((pair.entropy_beyond_linear/ell**4,
                              2*math.pi*pulse*trial.pop("weighted_excess"),
                              (8*math.pi**2/15)*(pulse-trial["action_stress"])))
            trial.update(endpoint_one_probability=3/8+pair.probability_change,
                         entropy_change=entropy, vacuum_modular_change=modular,
                         relative_entropy=interior_kl+pair.cut_relative_entropy,
                         area_over_ell4=area, entropy_over_ell4=entropy/ell**4, bulk_over_ell4=bulk,
                         geometric_over_ell4=area+(4*math.pi/15)*trial["metric_charge"],
                         balance_over_ell4=entropy/ell**4+area)
    return rows


def _tomography(rows):
    result = []
    for d in SCALES:
        shots = {r["shape"]:r["trials"] for r in rows if r["inverse_ell"] == d}
        for stock in shots["rest"]:
            # Keep the input list in the theorem's signed timelike frame order.
            arm = [shots[shape][stock] for shape in FRAMES.values()]
            area_tensor = decode([-15*r["area_over_ell4"]/(4*math.pi) for r in arm])
            entropy_tensor = decode([15*r["entropy_over_ell4"]/(8*math.pi**2) for r in arm])
            stress = decode([r["action_stress"] for r in arm])
            # QuantitativeSmallBall.small_ball_tensor_bound for eta=1, Lambda=0.
            bulk, geometric, balance = [max(abs(r[field]) for r in arm) for field in
                                       ("bulk_over_ell4", "geometric_over_ell4", "balance_over_ell4")]
            theorem_bound = (17/8)*(15/(4*math.pi))*(bulk+geometric+balance)
            metric = decode([r["metric_charge"] for r in arm])
            residual = max(abs(g-2*math.pi*t) for g,t in zip(metric, stress, strict=True))
            if residual > theorem_bound + 5e-12:
                raise ArithmeticError("matched finite error budget does not cover its tensor")
            result.append(dict(inverse_ell=d, stock=stock, from_areas=list(area_tensor),
                               from_entropies=list(entropy_tensor), from_action=list(stress),
                               largest_metric_action_error=max(abs(g-2*math.pi*t) for g,t in zip(area_tensor,stress,strict=True)),
                               largest_entropy_action_error=max(abs(s-t) for s,t in zip(entropy_tensor,stress,strict=True)),
                               largest_balance_over_ell4=balance, largest_bulk_over_ell4=bulk,
                               largest_geometry_over_ell4=geometric, finite_tensor_bound=theorem_bound,
                               metric_vs_action_tensor_residual=residual))
    return result


def counted_endpoints(rows):
    row, = [r for r in rows if r["shape"] == "xy" and r["inverse_ell"] == 32]
    trials = 2**40
    rng = np.random.Generator(np.random.PCG64(112605))
    radius = math.sqrt(math.log(600)/(2*trials))
    counts = {}
    intervals = {}

    def h(p):
        return -p*math.log(p)-(1-p)*math.log1p(-p)

    for stock in ("blank", "ordinary", "held"):
        hits = int(rng.binomial(trials, row["trials"][stock]["endpoint_one_probability"]))
        counts[stock] = hits
        p = hits/trials
        if not (0 < p-radius < p+radius < 0.5):
            raise ValueError("binary entropy monotonicity domain changed")
        intervals[stock] = (h(p-radius), h(p+radius))
    a, b = intervals["held"], intervals["ordinary"]
    delta = [a[0]-b[1], a[1]-b[0]]
    t = row["trials"]
    return dict(seed=112605, trials_per_stock=trials, one_counts=counts,
                held_minus_ordinary_entropy_99=delta,
                opposite_finite_rod_area=t["ordinary"]["proper_area"]-t["held"]["proper_area"],
                opposite_half_strength_area=t["ordinary"]["proper_area"]-t["half_defect"]["proper_area"])


def measure():
    densities = source_densities()
    ideal = _board(densities, finite=False)
    actual = _board(densities, finite=True)
    return dict(schema="oph.published_string_ball.v1", frames=FRAMES, scales=list(SCALES),
                ordinary_and_coincidence_densities={k:list(v) for k,v in densities.items()},
                rows=ideal, ideal_tomography=_tomography(ideal), finite_tomography=_tomography(actual),
                finite_rows=actual, counted_endpoints=counted_endpoints(actual))


def verify_saved(result):
    saved = json.loads((HERE/"results.json").read_text())
    for key in ("schema", "frames", "scales", "ordinary_and_coincidence_densities"):
        if saved[key] != result[key]:
            raise ValueError(f"source, frame or scale mismatch in {key}")
    counted = result["counted_endpoints"]
    old_counts = saved["counted_endpoints"]
    for key in ("seed", "trials_per_stock", "one_counts"):
        if old_counts[key] != counted[key]:
            raise ValueError(f"changed counted endpoints: {key}")
    for key in ("held_minus_ordinary_entropy_99", "opposite_finite_rod_area", "opposite_half_strength_area"):
        np.testing.assert_allclose(counted[key], old_counts[key], atol=2e-15, rtol=1e-9)
    for key in ("ideal_tomography", "finite_tomography"):
        for old, new in zip(saved[key], result[key], strict=True):
            if (old["stock"], old["inverse_ell"]) != (new["stock"], new["inverse_ell"]):
                raise ValueError("tensor source changed")
            for name in old.keys() - {"stock", "inverse_ell"}:
                np.testing.assert_allclose(new[name], old[name], atol=4e-12, rtol=1e-9)
    for name in ("rows", "finite_rows"):
        _compare_boards(saved[name], result[name])


def _compare_boards(saved_rows, actual_rows):
    for old, new in zip(saved_rows, actual_rows, strict=True):
        for key in ("shape", "inverse_ell", "whole_pairs"):
            if old[key] != new[key]:
                raise ValueError(f"proper ball selection changed in {key}")
        if not np.allclose(old["worldlines"], new["worldlines"], rtol=1e-14, atol=0):
            raise ValueError("source-independent material layout changed")
        if old["trials"].keys() != new["trials"].keys():
            raise ValueError("registered controls changed")
        for key, trial in old["trials"].items():
            for field, expected in trial.items():
                reading = new["trials"][key]
                tol = 1e-25
                if field in ("geometric_over_ell4", "balance_over_ell4"):
                    # The last sum cancels two normalized O(1) readings. Budget
                    # roundoff in those two operands, independently of radius.
                    other = ((4*math.pi/15)*reading["metric_charge"] if field == "geometric_over_ell4"
                             else reading["entropy_over_ell4"])
                    tol = 64*math.ulp(abs(reading["area_over_ell4"])+abs(other))
                if not math.isclose(reading[field], expected, abs_tol=tol, rel_tol=1e-9):
                    raise ValueError(f"changed {field} on {key}, {new['shape']}, 1/{new['inverse_ell']}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    result = measure()
    if args.write:
        (HERE/"results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        verify_saved(result)
    print("Verified 10 frames, 3 proper volumes, 3 source archives and 2 changed-law controls.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
