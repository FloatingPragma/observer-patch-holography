"""Compare the protected regional modular response with a native energy-orientation read."""

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import log, sqrt
from pathlib import Path

import numpy as np

from model.executable_model import federation
from .... import shared_variance, volume
from .. import measure as clock
from .. import word_entropy

HERE = Path(__file__).resolve().parent


def word_panel():
    _, _, _, source = volume.radar_source_and_hand(federation(0))
    laws, _ = shared_variance.source_words(source)
    reverse, fresh, shared = np.rint(2592 * laws).astype(int)
    np.testing.assert_array_equal(np.stack([reverse, fresh, shared]) / 2592, laws)
    delta = fresh - reverse
    np.testing.assert_array_equal(delta, shared - fresh)

    # The reader uses the cyclic order 0 -> 32 -> 33 -> 0 on the two actual
    # seam energies in one word. The writer's concealed tag is never scored.
    energies = tuple(source["source_word_energies"])
    ranks = {energy: i for i, energy in enumerate((0, 32, 33))}
    orientations = np.array(
        [(0, 1, -1)[(ranks[b] - ranks[a]) % 3] for a in energies for b in energies]
    )
    score = np.outer(orientations, orientations)
    orientation_slope = Fraction(int(np.sum(score * delta)), 2592)
    # The product reference has additive log weights, annihilated by delta.
    # Grouping the remaining log terms yields an exact log of a rational.
    exponents = Counter()
    for p, d in zip(fresh.flat, delta.flat, strict=True):
        exponents[Fraction(int(p - d), int(p))] -= int(d)
    exponents = {x: y for x, y in exponents.items() if y}
    if any(y % 18 for y in exponents.values()):
        raise ValueError("Grouping does not have whole powers at the declared scale")
    log_base = Fraction(1)
    for ratio, exponent in exponents.items():
        log_base *= ratio ** (exponent // 18)

    held = json.loads((shared_variance.HERE / "shared_variance.json").read_text())
    records = held["connected_modular_charge"]
    trials = records["shots_per_setting"]
    counts = np.asarray(records["held_counts"], dtype=np.int64)
    np.testing.assert_array_equal(counts.sum(axis=(1, 2)), [trials, trials])
    sums = [int(np.sum(row * score)) for row in counts]
    contrast = Fraction(sums[1] - sums[0], trials)
    # The score lies in [-1,1]. Union bound on both Hoeffding intervals.
    radius = 2 * sqrt(2 * log(4 / 0.01) / trials)
    return {
        "orientation_on_energy_labels": orientations.reshape(3, 3).tolist(),
        "seam_energy_labels": list(energies),
        "conditional_orientation_slope": str(orientation_slope),
        "conditional_entropy_log_base": str(log_base),
        "conditional_entropy_log_multiplier": str(Fraction(18, 2592)),
        "conditional_log_exponents_before_division_by_2592": {
            str(x): y for x, y in sorted(exponents.items())
        },
        "held_shared_minus_fresh": {
            "trials_per_arm": trials,
            "signed_score_sums_fresh_shared": sums,
            "mean_difference": float(contrast),
            "simultaneous_99pct_radius": radius,
        },
    }


def run():
    words = word_panel()
    log_base = Fraction(words["conditional_entropy_log_base"])
    multiplier = Fraction(words["conditional_entropy_log_multiplier"])
    slope = word_entropy.interval(multiplier) * word_entropy.IV.ln(
        word_entropy.interval(log_base)
    )
    ancestor = json.loads((clock.HERE / "ancestor_limit.json").read_text())
    turns = {int(face): hop for face, hop in ancestor["sector_offsets"].items()}
    first = next(i for i, turn in turns.items() if turn == 0)
    second = min(i for i, turn in turns.items() if turn == 2)
    panels = []
    for row in ancestor["midpoint_surveys"]:
        radius = row["radius"]
        total = 20 * row["frequency"] ** 2
        ell = Fraction(radius, row["eccentricity"])
        for region, inner in (("ball", -1), ("annulus", radius // 2)):
            counts = [
                clock.corner_sector(radius - turns[i])
                - clock.corner_sector(inner - turns[i])
                for i in (first, second)
            ]
            if region == "annulus" and counts != [
                row["ancestor_populations"][i] for i in (first, second)
            ]:
                raise ValueError(
                    "Sector formula disagrees with the measured enrollments"
                )
            eligible = Fraction(2 * counts[0] * counts[1], total * (total - 1))
            panels.append(
                {
                    "level": row["level"],
                    "region": region,
                    "radius": radius,
                    "sector_counts": counts,
                    "eligible": str(eligible),
                    "normalized_eligible": str(eligible / ell**4),
                }
            )
    leading = {}
    for region, face_fraction in (
        ("ball", Fraction(1, 4)),
        ("annulus", Fraction(3, 16)),
    ):
        coefficient = 2 * face_fraction**2 * Fraction(6**4, 20**2)
        leading[region] = {
            "eligible_per_fourth": str(coefficient),
            "orientation_slope_per_fourth": str(
                coefficient * Fraction(words["conditional_orientation_slope"])
            ),
            "entropy_slope_per_fourth_nats": word_entropy.endpoints(
                word_entropy.interval(coefficient) * slope
            ),
        }
    return {
        "source": "Two distinct uniform word stations with one fair preparation tag per original face. Only the correlation of the tags on the selected two faces changes.",
        "sectors": [first, second],
        "sector_offsets": [turns[first], turns[second]],
        "words": words,
        "surveys": panels,
        "limits": leading,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
