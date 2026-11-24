"""Check the long-wavelength Klein–Gordon relation on the recorded spatial phases."""

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def summarize(meter, holdout):
    mass = meter["results"]["predicted_m"]["value"]
    speed_squared = meter["results"]["predicted_c_squared"]["value"]
    assert mass > 0 and speed_squared > 0
    panels = []
    for board in holdout["boards"]:
        low = []
        all_errors = []
        for mode in board["observer"]["modes"]:
            k, measured = mode["k"], mode["omega"]
            continuum = math.sqrt(mass**2 + speed_squared * k**2)
            lattice = math.sqrt(mass**2 + 4 * speed_squared * math.sin(k / 2) ** 2)
            correction_bound = speed_squared * k**4 / (24 * mass)
            assert -1e-12 <= continuum - lattice <= correction_bound + 1e-12
            error = abs(continuum - measured)
            all_errors.append(error)
            if 0 < k <= 0.5:
                low.append((k, error, error / measured))
        panels.append(
            dict(
                cells=20 * 4 ** board["level"],
                low_modes=len(low),
                measured_k_range=[min(x[0] for x in low), max(x[0] for x in low)],
                max_low_error=max(x[1] for x in low),
                max_low_relative_error=max(x[2] for x in low),
                max_all_error=max(all_errors),
            )
        )
    return panels


if __name__ == "__main__":
    meter = json.loads((HERE / "coupling_meter.json").read_text())
    held = json.loads((HERE / "holdout.json").read_text())
    print(json.dumps(summarize(meter, held), indent=2))
