"""Three simultaneous records of a dust source and held exterior clock transfers."""
import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import json
from math import inf, isclose, nextafter, pi, sqrt
from pathlib import Path

from model.candidate_laws.interval_echo_partition_repair.published_strings.geometry import FRAME
from model.candidate_laws.interval_echo_partition_repair.published_strings.memories import Pair
from matched_einstein_limit.frame import apply, design, inverse

from . import survey
from .body import Body, archive_bodies, collapse_rate, gravitational_mass, release_calibration


HERE = Path(__file__).resolve().parent
SCALES = (512, 2048, 8192)
NAMES = ("rest", "+x", "-x", "+y", "-y", "+z", "-z", "xy", "xz", "yz")


def _counts(value):
    return F(value).limit_denominator(66)


def _direction(frame):
    dot = sum(n*w for n, w in zip((1, 2, 4), frame[1:], strict=True))/frame[0]
    return dot, dot*dot/21


def _fraction(f):
    return D(f.numerator)/f.denominator


def measured_frame(stamps, denominator):
    """Four-velocity from six phase differences per two ticks of the moving clock."""
    if len(stamps) != 6 or any(len(arm) != 2 for arm in stamps):
        raise ValueError("six opposing phase pairs are required")
    rates = tuple(F(end-start, 2*denominator**8) for start, end in stamps)
    return (sum(rates)/6, *((rates[i+1]-rates[i])/2 for i in (0, 2, 4)))


def frame_inverse(observers, denominator):
    frames = [measured_frame(obs["frame_phase_stamps"], denominator) for obs in observers]
    return inverse(design(frames))


def inverse_gain(matrix):
    return nextafter(float(max(sum(map(abs, row)) for row in matrix)), inf)


def regional_work(tickets, weight, gamma, include_collar):
    """Weighted work from enclosed mass differences, rod weights and measured boost."""
    masses = [F((hi-lo)*64, 7*n) for n, lo, hi in tickets]
    stock = sum(masses) if include_collar else masses[0]
    return 3*_fraction(stock*gamma**2)*weight


def _crossing(frame, radius, h):
    """Common comoving ray with direction (1,2,4)/sqrt(21), including its time tilt."""
    dot, dot2 = _direction(frame)
    stretch = 1/(1-_fraction(dot2)).sqrt()
    w = stretch*_fraction(dot)/D(21).sqrt()
    return stretch*survey.comoving_distance(radius, h, w)


def _trial(body, pairs, calibration, *, read_collar=True):
    geometric = _counts(body.gravitating_count)
    intervals, rows = {}, []
    for denominator in SCALES:
        ell = D(1)/denominator
        h = collapse_rate(gravitational_mass(geometric))
        tickets = [survey.mass_tally(_counts(x), denominator) for x in (body.ordinary, body.collar)]
        masses = [F((hi-lo)*64, n*7) for n, lo, hi in tickets]
        metered = sum(masses)
        readings = []
        for label, frame in zip(NAMES, FRAME, strict=True):
            gamma = frame[0]
            ball = survey.read_ball(geometric, gamma, denominator)
            blank = survey.read_ball(F(0), gamma, denominator)
            area = 4*pi*float((ball.area_over_solid_angle-blank.area_over_solid_angle)/ell**4)
            beam_stamps, _ = survey.frame_clocks(frame, h, ell)
            velocity = measured_frame(beam_stamps, denominator)
            work = regional_work(tickets, ball.mass_weight_over_solid_angle, velocity[0], read_collar)
            memory = Pair.from_weighted_work(float(work))
            cut = (denominator, _direction(frame)[1])
            if cut in intervals and intervals[cut] != memory:
                raise ArithmeticError("two observers assigned different operations to the same memory")
            intervals[cut] = memory
            boundary = _crossing(velocity, ball.radius, ball.measured_h)
            cut_pairs = [k for k, (left, right) in pairs.items() if left < boundary < right]
            if cut_pairs != [cut]:
                raise ArithmeticError(f"expected one common material cut at {denominator}, {label}")
            enclosed = [key for key, (_, right) in pairs.items() if right < boundary]
            stress = 3*float(metered*velocity[0]**2)/(4*pi)
            curvature = float(6*gravitational_mass(geometric)*gamma**2)
            entropy = memory.cut_entropy_change/float(ell**4)
            readings.append(dict(frame=label, area_over_r4=area, entropy_over_r4=entropy,
                                 bulk_error_over_r4=entropy-(8*pi*pi/15)*stress,
                                 area_error_over_r4=area+(4*pi/15)*curvature,
                                 stationarity_over_r4=area+entropy, independent_stress=stress,
                                 calibrated_balance_over_r4=calibration[denominator]["eta"]*area+entropy,
                                 volume_stop_error=float(4*D(str(pi))*(ball.volume_over_solid_angle-ell**3/3)),
                                 frame_phase_stamps=beam_stamps,
                                 proper_mass_weight=str(ball.mass_weight_over_solid_angle),
                                 max_velocity_component_error=float(max(abs(read-w)
                                     for read, w in zip(velocity, frame, strict=True))),
                                 membership_margin_over_r=float(min(abs(spot-boundary)
                                     for ends in pairs.values() for spot in ends)/ell),
                                 cut_key=cut, whole_keys=enclosed))
            if label == "xy" and denominator == SCALES[0]:
                readings[-1]["clock_sample"] = dict(tick=str(ball.sample_tick),
                                                     stamps=list(ball.sample_timestamps))
        matrix = frame_inverse(readings, denominator)
        tensor = apply(matrix, [-15*r["area_over_r4"]/(4*pi) for r in readings])
        rows.append(dict(inverse_radius=denominator, observers=readings, mass_tallies=tickets,
                         scale_record=ball.scale_record, pulse_uses_joint_stock=read_collar,
                         measured_inverse_gain=inverse_gain(matrix), rod_einstein_tensor=list(tensor)))
    for row in rows:
        for obs in row["observers"]:
            cut = intervals[obs.pop("cut_key")]
            obs["regional_relative_entropy"] = cut.cut_relative_entropy + sum(
                intervals[k].whole_modular_change for k in obs.pop("whole_keys"))
    return rows


def registered_pairs():
    """The same twenty-one parked pairs shared by all enrolled stocks and frames."""
    with localcontext() as ctx:
        ctx.prec = 75
        pairs = {}
        for denominator in SCALES:
            for frame in FRAME:
                group = _direction(frame)[1]
                crossing = 1/(D(denominator)*(1-_fraction(group)).sqrt())
                pairs[denominator, group] = (crossing*159/160, crossing*161/160)
        return pairs


def measure():
    bodies = archive_bodies()
    pairs = registered_pairs()
    with localcontext() as ctx:
        ctx.prec = 75
        calibration = {den: release_calibration(_counts(bodies["ordinary"].ordinary), den) for den in SCALES}
        cells = {name: _trial(body, pairs, calibration)
                 for name, body in bodies.items()}
        held = bodies["held"]
        half = Body(held.ordinary, held.collar, 0.5)
        cells["half_collar_force"] = _trial(half, pairs, calibration)
        cells["blind_collar_work"] = _trial(held, pairs, calibration, read_collar=False)
    exterior = {name: body.exterior_records() for name, body in {**bodies, "half_collar_force": half}.items()}
    ordinary, held = bodies["ordinary"], bodies["held"]
    # Freeze the ordinary gain from the finer released-probe record before using the held clocks.
    ordinary_acc = exterior["ordinary"]["signed_radius_acceleration"][-1]
    gain = -4*ordinary_acc/ordinary.ordinary
    m = held.ordinary + held.collar
    forecast = dict(ordinary_release_gain=gain, held_signed_release=-gain*m/4,
                    held_clock_ratio=sqrt((1-gain*m)/(1-gain*m/2)))
    return dict(schema="oph.published_spherical_dust.v1", masses={k:[b.ordinary,b.collar] for k,b in bodies.items()},
                registered_pairs=len(pairs), ordinary_calibration={str(k):v for k,v in calibration.items()},
                ball_surveys=cells, exterior=exterior, held_forecast=forecast)


def verify(result):
    archive = json.loads((HERE/"results.json").read_text())
    if result.keys() != archive.keys() or result["masses"] != archive["masses"] or result["registered_pairs"] != 21:
        raise ValueError("unexpected material archive or enrollment")
    def compare(new, old, path):
        if isinstance(old, dict):
            if set(new) != set(old):
                raise ValueError(f"missing record at {path}")
            for key in old:
                compare(new[key], old[key], f"{path}.{key}")
        elif isinstance(old, list):
            for i, (x,y) in enumerate(zip(new, old, strict=True)):
                compare(x, y, f"{path}[{i}]")
        elif isinstance(old, float):
            # Exterior release divides a difference of binary64 radius ticks by dt².
            if ".signed_radius_acceleration" in path or ".held_forecast." in path:
                absolute = 3e-9
            elif "_over_r4" in path or "rod_einstein_tensor" in path:
                absolute = 8e-14
            else:
                absolute = 1e-75
            if not isclose(new, old, rel_tol=3e-11, abs_tol=absolute):
                raise ValueError(f"changed measured value at {path}")
        elif isinstance(old, int) and path.startswith("dust.exterior.") and any(
                field in path for field in (".received[",".released_radius_ticks[")):
            if abs(new-old) > 4:
                raise ValueError(f"changed exterior clock or radius at {path}")
        elif new != old:
            raise ValueError(f"changed record at {path}")
    compare(result, archive, "dust")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    result = measure()
    if args.write:
        (HERE/"results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        verify(result)
    held = result["ball_surveys"]["held"][-1]
    bound = max(abs(o["stationarity_over_r4"]) for o in held["observers"])
    print(f"Spherical dust: 10 rest frames, 3 scales, 21 shared pairs, held |s+a|/ell^4 <= {bound:.3g}.")
    return result


if __name__ == "__main__":
    main()
