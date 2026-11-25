"""Calibrate the private-reference phase and compare it with the orbital colour centre."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier
from .measure import HERE, STARTS, bank

TURNS = (1 / 8, 1 / 4, 1 / 2, 1)
SCALES = {"native": 1.0, "half_strength": 0.5}
WEIGHTS = {"empty": 0.0, "unit": 1.0, "up": 2 / 3, "down": -1 / 3}


def probability(turn, weight, scale, quadrature):
    angle = 2 * np.pi * turn * scale * weight
    return float((1 + (np.cos(angle) if quadrature == "X" else np.sin(angle))) / 2)


def run():
    receipt = json.loads((HERE / "results.json").read_bytes())
    rng = np.random.default_rng(13116021)
    ports = np.flatnonzero(carrier.federation(0).degree[:12] == 5)[:3].tolist()
    # Exactly one of the three binary ports holds a load on every codeword.
    words = np.eye(3, dtype=int)
    service = candidate.KAPPA * (words**2).sum(axis=1)
    if len(ports) != 3 or not np.array_equal(service, [8] * 3):
        raise RuntimeError(
            "The native service does not implement the projected number coupling"
        )
    source_p = receipt["source"]["filter_acceptances"][-1]
    tri = receipt["preparations"][3]["triplet"]
    source_modes = [("vacuum", 1.0, 1.0, 0), ("standard", 1.0, 1.0, 1)]
    source_modes.extend((f"constituent_{i}", source_p, tri, 1) for i in range(3))
    results = {}
    # One shared budget for the means in 2 strengths, 5 preparations,
    # 4 waits and 2 quadratures. Condition on the pre-setting admissions.
    for variant, scale in SCALES.items():
        for tag, source, admit, word in source_modes:
            rows = []
            fits = dict.fromkeys(WEIGHTS, True)
            for turn in TURNS:
                for quad in ("X", "Y"):
                    p = probability(turn, word, scale, quad)
                    bins = bank(
                        rng,
                        [
                            1 - source,
                            source * (1 - admit),
                            source * admit * p,
                            source * admit * (1 - p),
                        ],
                    )
                    n = bins[2] + bins[3]
                    radius = float(np.sqrt(np.log(2 * 160 / 0.01) / (2 * n)))
                    observed = bins[2] / n
                    for name, weight in WEIGHTS.items():
                        fits[name] &= bool(
                            abs(probability(turn, weight, scale, quad) - observed)
                            <= radius
                        )
                    rows.append(
                        dict(
                            turn=turn,
                            wait=2 * np.pi * turn / candidate.KAPPA,
                            quadrature=quad,
                            counts=bins,
                            bright_mean=observed,
                            joint_radius=radius,
                        )
                    )
            results[variant + "/" + tag] = dict(rows=rows, covered=fits)
    return dict(
        receiver_private_ports=ports,
        one_load_service_spectrum=service.tolist(),
        starts=STARTS,
        arrays="source refused, triplet refused, bright, dark",
        controls=results,
        earlier_colour_centre_dark=receipt["material_half_turn"]["fringes"][
            "centre_one"
        ]["dark_given_triplets"],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("number.json")
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from ..open_release.measure import compare

        compare(json.loads(json.dumps(result)), json.loads(path.read_text()))
    print(
        "Independent material-number calibration and constituent phase banks verified."
    )


if __name__ == "__main__":
    main()
