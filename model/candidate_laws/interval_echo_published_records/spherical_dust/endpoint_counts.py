"""Forecast the dust energy and equal-volume area from counted endpoint records."""
import argparse
from decimal import Decimal as D, localcontext
from hashlib import sha256
import json
import math

import numpy as np

from model.candidate_laws.interval_echo_partition_repair.published_strings.counted_tower.experiment import (
    phase_numerator, read_count, settings,
)
from .experiment import (
    HERE, NAMES, SCALES, apply, frame_inverse, inverse_gain, measure, measured_frame, regional_work,
)


STOCKS = ("blank", "ordinary", "held", "half_collar_force", "blind_collar_work")
SEED = 26112631
COEFFICIENT = 4*math.pi/15


def count_digest(levels):
    raw = [[row["design"]["inverse_radius"], stock, frame, *[bank[k] for k in ("ones", "zeros", "dyadic_one_numerator")]]
           for row in levels for stock in STOCKS for frame in NAMES
           for bank in [row["banks"][stock][frame]]]
    return sha256(json.dumps(raw, separators=(",", ":")).encode("ascii")).hexdigest()


def evaluate(data, frozen=None):
    """Regenerate each transmitter; decode every deposited zero and one when replaying."""
    levels = []
    for level, den in enumerate(SCALES, start=2):
        cfg = settings(level, banks_per_level=len(STOCKS)*len(NAMES))
        if cfg["inverse_radius"] != den:
            raise ValueError("mismatched coherent and proper-radius schedules")
        banks, summary, signals, tensors, gains = {}, {}, {}, {}, {}
        for i, stock in enumerate(STOCKS):
            survey = data["ball_surveys"][stock][level-2]
            banks[stock], signals[stock] = {}, []
            for j, (frame, obs) in enumerate(zip(NAMES, survey["observers"], strict=True)):
                if obs["frame"] != frame or survey["inverse_radius"] != den:
                    raise ValueError("rod appointment differs from the endpoint bank")
                gamma = measured_frame(obs["frame_phase_stamps"], den)[0]
                with localcontext() as ctx:
                    ctx.prec = 70
                    work = regional_work(survey["mass_tallies"], D(obs["proper_mass_weight"]), gamma,
                                         survey["pulse_uses_joint_stock"])
                    numerator = phase_numerator(work, cfg)
                shots = cfg["shots_per_bank"]
                if frozen is None:
                    rng = np.random.default_rng(np.random.SeedSequence([SEED, level, i, j]))
                    batch = 2**42
                    hits = int(sum(map(int,rng.binomial(batch, numerator/2**cfg["probability_bits"], shots//batch))))
                    bank = dict(ones=hits, zeros=shots-hits, dyadic_one_numerator=numerator)
                else:
                    bank = frozen[level-2]["banks"][stock][frame]
                    if bank["dyadic_one_numerator"] != numerator:
                        raise ValueError("deposited finite emitter differs from the native material read")
                banks[stock][frame] = bank
                signals[stock].append(read_count(bank["ones"],bank["zeros"],cfg))
            matrix = frame_inverse(survey["observers"], den)
            tensors[stock] = apply(matrix, [v[0]/(2*math.pi*COEFFICIENT) for v in signals[stock]])
            gains[stock] = inverse_gain(matrix)

        for label, upper, lower in [(s,s,"blank") for s in STOCKS[1:]] + [("collar","held","ordinary")]:
            bounds, lower_rejection, bulk = [], [], []
            eta = data["ordinary_calibration"][str(den)]["eta"]
            for j in range(10):
                top, bot = [signals[s][j] for s in (upper, lower)]
                lo, hi = top[1]-bot[2], top[2]-bot[1]
                top_row, bot_row = [data["ball_surveys"][s][level-2]["observers"][j] for s in (upper, lower)]
                area = top_row["area_over_r4"]-bot_row["area_over_r4"]
                stress = top_row["independent_stress"]-bot_row["independent_stress"]
                bounds.append(max(abs(lo+eta*area),abs(hi+eta*area)))
                lower_rejection.append(max(0,lo+eta*area,-hi-eta*area))
                bulk.append(max(abs(lo-2*math.pi*COEFFICIENT*stress),abs(hi-2*math.pi*COEFFICIENT*stress)))
            tensor_difference = [a-b for a, b in zip(tensors[upper], tensors[lower], strict=True)]
            radius = (gains[upper]+gains[lower])*cfg["entropy_difference_radius"]/(4*math.pi*COEFFICIENT)
            summary[label] = dict(entropy_stress_tensor=tensor_difference,
                                  component_sampling_radius=math.nextafter(radius, math.inf),
                                  bulk_upper=max(bulk), stationarity_upper=max(bounds),
                                  stationarity_lower=max(lower_rejection))
        levels.append(dict(design=cfg, banks=banks, comparisons=summary))
    return dict(schema="oph.dust_endpoint_counts.v1", seed=SEED, counts_sha256=count_digest(levels),
                simultaneous_coverage=0.99, signed_frame_order=list(NAMES), levels=levels)


def replay(data=None, *, write=False):
    data = measure() if data is None else data
    path = HERE/"endpoint_counts.json"
    if write:
        result = evaluate(data)
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        old = json.loads(path.read_text())
        result = evaluate(data, old["levels"])
        if any(result[k] != old[k] for k in ("schema","seed","counts_sha256","simultaneous_coverage","signed_frame_order")):
            raise ValueError("counted material or moving-frame order changed")
        for actual, saved in zip(result["levels"], old["levels"], strict=True):
            if actual["banks"] != saved["banks"] or actual["design"].keys() != saved["design"].keys():
                raise ValueError("changed dust phase banks or sampling design")
            for key, value in actual["design"].items():
                if not math.isclose(value,saved["design"][key],rel_tol=5e-14,abs_tol=0):
                    raise ValueError("changed counted dust instrument precision")
            for stock, row in actual["comparisons"].items():
                for key, value in row.items():
                    if not np.allclose(value,saved["comparisons"][stock][key],rtol=3e-11,atol=3e-13):
                        raise ValueError(f"changed counted {stock} {key}")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    result = replay(write=parser.parse_args().write)
    print(f"{len(result['levels'])*len(NAMES)*len(STOCKS)} dust endpoint banks checked.")
