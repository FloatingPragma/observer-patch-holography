"""Compare the energy-matched, probe-independent repair with disk and cosmic densities."""

import argparse
from fractions import Fraction as Q
import json

from . import pooled_fluctuations as fluctuations

pooled = fluctuations.pooled


def run():
    native = pooled.scales.budget.LAB / "landmark_arrival.json"
    measured = json.loads(native.read_bytes())
    preparation = measured["ordinary_energy_matched"]["preparations"]["entropy_stationary"]
    once = Q(preparation["uniform_proposal_acceptance"])
    coeff = [float(Q(x)) for x in measured["ordinary_coefficients"]]
    _, mu = pooled.scales.source_data()
    result = {}
    for name, acceptance in (("one_proposal", once), ("renewed_until_accepted", Q(1))):
        sky = pooled.run(native, float(acceptance))
        result[name] = dict(acceptance=str(acceptance), **sky)
        result[name]["matter_fluctuations"] = fluctuations.forecasts(
            sky["homogeneous"], coeff, mu, float(acceptance)
        )
    return result


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    path = pooled.HERE / "landmark_abundance.json"
    if args.write:
        path.write_text(raw)
    else:
        pooled.sky.neighbor.inventory.rp.galaxy.check(
            json.loads(path.read_bytes()), out
        )
    for name, row in out.items():
        print(
            name,
            max(v["within_two"] for v in row["scores"].values()),
            max(len(v["complete"]) for v in row["scores"].values()),
            min(v["predicted_sigma8"] for v in row["matter_fluctuations"].values()),
        )
