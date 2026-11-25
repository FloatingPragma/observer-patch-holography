"""Independent offline arithmetic check of the archived OW-30 v3 headline.

Only Python's standard library is used. The meter is reconstructed from joint
integer counts. Frequencies, fitted wavelengths and spatial residuals are
inputs from the archived validation summary; their acquisition and reduction
require the separate frozen replay. No candidate admission or row award is
inferred. Saved evidence is read without alteration.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
METER = HERE.parents[1] / "row_30_collective/coupling_meter.json"
SEEDS = [785029, 514823, 690271]
TOLERANCE = 0.0048


class Rejected(ValueError):
    """A recomputed headline gate rejects the supplied measured values."""


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    def unique(pairs):
        out = dict(pairs)
        require(len(out) == len(pairs), "duplicate JSON key")
        return out

    return json.loads(
        Path(path).read_text(), object_pairs_hook=unique,
        parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f"invalid {x}")),
    )


def unwrap(values):
    out = [values[0]]
    for old, value in zip(values, values[1:]):
        delta = (value - old + math.pi) % (2 * math.pi) - math.pi
        out.append(out[-1] + delta)
    return out


def meter_from_counts(bench):
    """Reconstruct the clock and four phase slopes without producer imports."""
    counts = bench["all_counts"]
    require(len(counts) == 4 and bench["joint_words"] == [
        [i >> 2, (i >> 1) & 1, i & 1] for i in range(8)
    ], "incorrect meter channels or word order")
    samples = len(counts[0][0])
    require(samples == 161, "incorrect meter hold count")
    rates = []
    for probe in counts:
        require(len(probe) == 4, "incomplete meter orientations")
        orientations = []
        for holds in probe:
            require(len(holds) == samples, "incomplete meter holds")
            values = []
            for word in holds:
                require(len(word) == 8 and all(type(n) is int and n >= 0 for n in word),
                        "invalid joint counts")
                require(sum(word) == bench["trials"] == 2**17, "meter shot budget differs")
                values.append([sum(n for i, n in enumerate(word) if i & bit) / sum(word)
                               for bit in (4, 2, 1)])
            orientations.append(values)
        rates.append(orientations)
    calibration = bench["hand_calibration"]
    require(bench["calibration_trials"] == 2**22 and len(calibration) == 4
            and all(type(n) is int and 0 <= n <= 2**22 for n in calibration),
            "invalid hand calibration")
    cal = [n / 2**22 for n in calibration]
    offset = sum(cal) / 4
    a, b = (cal[0] - cal[2]) / 2, (cal[1] - cal[3]) / 2
    determinant = a * a + b * b
    require(determinant > 0, "unresolved hand calibration")
    phase = []
    for t in range(samples):
        x, y = [sum(rates[p][o][t][axis] for p in range(4) for o in range(4)) / 16 - offset
                for axis in (0, 1)]
        phase.append(math.atan2((b * x - a * y) / determinant,
                                (a * x + b * y) / determinant))
    times = unwrap(phase)
    center = sum(times) / samples
    times = [t - center for t in times]
    denominator = sum(t * t for t in times)
    require(denominator > 0, "unresolved clock")
    slopes = []
    for probe in rates:
        angles = unwrap([math.atan2(probe[1][t][2] - probe[3][t][2],
                                    probe[0][t][2] - probe[2][t][2])
                         for t in range(samples)])
        require(max(abs(y - x) for x, y in zip(angles, angles[1:])) < 0.9 * math.pi,
                "unresolved probe phase")
        slopes.append(sum(x * y for x, y in zip(times, angles)) / denominator)
    h, j = slopes[0] / 2, slopes[1] / 4
    require(j > h > 0 and max(abs(x) for x in slopes[2:]) < 0.005,
            "coupling meter or switched-off controls fail")
    return dict(h=h, j=j, m=2 * abs(j - h), c_squared=4 * j * h)


def evaluate(receipt, meter):
    """Compute dispersion and spatial headline gates from supplied measurements."""
    require(len(receipt["boards"]) == 3, "three enrolled boards required")
    output = []
    for level, board in enumerate(receipt["boards"]):
        enrollment, spatial = board["enrollment"], board["spatial"]
        require(enrollment["level"] == level and enrollment["seed"] == SEEDS[level]
                and enrollment["ineligible"] is None, "validation enrollment differs")
        modes = board["dispersion"]["modes"]
        n = len(enrollment["route"])
        channels = min(8, max(4, (n + 2) // 4))
        require([r["launch"] for r in modes] == list(range(1, channels + 1)),
                "missing or reordered modes")
        k = [r["k"] for r in modes]
        omega = [r["omega"] for r in modes]
        require(k == spatial["even_trained_k"] and omega == spatial["omega_from_pair_counts"],
                "spatial and temporal inputs disagree")
        require(all(type(x) in (int, float) and math.isfinite(x) for x in k + omega)
                and all(0 < x < math.pi for x in k),
                "invalid measured modes")
        m, c2 = meter["m"], meter["c_squared"]
        q2 = [2 - 2 * math.cos(x) for x in k]
        lattice = [math.sqrt(m * m + c2 * x) for x in q2]
        continuum = [math.sqrt(m * m + c2 * x * x) for x in k]
        tangent = [m + c2 * x / (2 * m) for x in q2]
        low = [i for i, x in enumerate(k) if 0 < x <= 0.50]
        lattice_error = max(abs(x - y) for x, y in zip(lattice, omega))
        continuum_error = max((abs(continuum[i] - omega[i]) for i in low), default=math.inf)
        def rms(forecast):
            return math.sqrt(sum(((x - forecast[0]) - (y - omega[0]))**2
                                 for x, y in zip(forecast[1:], omega[1:])) / (len(omega) - 1))
        pairs = 2 * channels - 3
        degrees = (n - 1) * 2 * pairs
        tail = math.log(10**5)
        cutoff = degrees + 2 * math.sqrt(degrees * tail) + 2 * tail
        require(spatial["projected_degrees_of_freedom"] == degrees
                and math.isclose(spatial["quadratic_cutoff"], cutoff, abs_tol=1e-8),
                "spatial cutoff differs from terms")
        gates = dict(
            lattice=min(omega) > 0 and lattice_error < TOLERANCE,
            continuum=len(low) >= 2 and max(k[i] for i in low) - min(k[i] for i in low) >= 0.02
                      and continuum_error < TOLERANCE,
            shape=rms(lattice) < 0.5 * rms(tangent),
            profile=0 <= spatial["candidate_q"] <= cutoff,
            launch_foil=spatial["launch_label_q"] > cutoff,
        )
        if not all(gates.values()):
            raise Rejected(f"board {level}: {[key for key, ok in gates.items() if not ok]}")
        output.append(dict(cells=20 * 4**level, frequencies=len(modes), low_k_lines=len(low),
                           lattice_max=lattice_error, continuum_max=continuum_error,
                           spatial_q=spatial["candidate_q"], spatial_cutoff=cutoff))
    return dict(scope="Archived local dispersion and spatial-summary arithmetic only",
                source_commit=receipt["draw"]["binding"]["source_commit"], meter=meter,
                boards=output, frequencies=sum(b["frequencies"] for b in output),
                low_k_lines=sum(b["low_k_lines"] for b in output),
                filing="HOLD: official candidate admission and public confirmation required")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", type=Path, default=HERE)
    parser.add_argument("--meter", type=Path, default=METER)
    args = parser.parse_args(argv)
    try:
        package, receipt = (read(args.packet / n) for n in ("PACKAGE.json", "validation.json"))
        binding = receipt["draw"]["binding"]
        require(receipt["draw"]["seeds"] == SEEDS and binding["mode"] == "validation"
                and binding["anchor"] is None, "author-side validation draw required")
        require(all(binding[key] == package[key] for key in ("source_commit", "candidate_sha256"))
                and binding["candidate_pin"] == package["source_commit"]
                and binding["protocol_sha256"] == package["task_sha256"], "packet binding differs")
        require(hashlib.sha256(args.meter.read_bytes()).hexdigest() ==
                receipt["boards"][0]["dispersion"]["meter"]["source_sha256"],
                "meter count bytes differ from the archived binding")
        result = evaluate(receipt, meter_from_counts(read(args.meter)))
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    except Rejected as error:
        print(f"headline rejected: {error}", file=sys.stderr)
        return 2
    except (ValueError, KeyError, TypeError, IndexError, OSError) as error:
        print(f"headline invalid: {error}", file=sys.stderr)
        return 4


if __name__ == "__main__":
    raise SystemExit(main())
