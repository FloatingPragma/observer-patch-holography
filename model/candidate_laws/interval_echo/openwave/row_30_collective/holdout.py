"""Score a frozen independent material bank using disjoint coupling-meter records."""

import argparse
import json
from pathlib import Path
import numpy as np
from .experiment import HERE, K, DT, SAMPLES, grow, spectroscopy, transport
from .experiment import verify
from .formation import prepare
from model.executable_model import federation

SEED = 927016


def score(rows, obs, live, disabled, m, c2):
    e = np.array([r["omega"] for r in rows])
    k = np.array([r["k"] for r in rows])
    q2 = 2 - 2 * np.cos(k)
    lattice = np.sqrt(m * m + c2 * q2)
    par = m + c2 * q2 / (2 * m)
    cont = np.sqrt(m * m + c2 * k * k)

    def rms(x):
        y = (x - x[0]) - (e - e[0])
        return float(np.sqrt(np.mean(y[1:] ** 2)))

    controls = list(obs["same_spatial_source_changed_contrast"].values())

    def far(d):
        return max(
            (
                f["distant_max_count_difference"] / d["trials_per_arm"]
                for f in d["flights"][1:]
            ),
            default=0,
        )

    return dict(
        lattice_max=float(max(abs(lattice - e))),
        continuum_max=float(max(abs(cont - e))),
        sum_closure=obs["withheld_pair_sum_max_error"],
        shape_rms_lattice=rms(lattice),
        shape_rms_parabola=rms(par),
        shape_rms_continuum=rms(cont),
        min_native_amplitude=min(x["native_sum_amplitude"] for x in controls),
        max_changed_amplitude=max(x["changed_at_same_sum"] for x in controls),
        far_count_difference=far(live),
        far_changed_count_difference=far(disabled),
    )


def run():
    meter = json.loads((HERE / "coupling_meter.json").read_text())["results"]
    m, c2 = (meter[x]["value"] for x in ("predicted_m", "predicted_c_squared"))
    out = []
    for lev in range(3):
        commissioned = prepare(lev, SEED)
        board = federation(lev)
        route = max((p for closed, p in grow(board, SEED) if not closed), key=len)
        obs = spectroscopy(len(route), np.random.default_rng([SEED, lev, 47]))
        live, off = (transport(route, j, SEED) for j in (K, 0))
        scored = score(obs["modes"], obs, live, off, m, c2)
        cells = len(set(np.ravel(route) // 12))
        scored["operational_pass"] = bool(
            cells > 1
            and scored["lattice_max"] < 0.0048
            and scored["sum_closure"] < 0.0048
            and scored["min_native_amplitude"] > 0.40
            and scored["max_changed_amplitude"] < 0.02
            and scored["far_count_difference"] > 0.02
            and scored["far_changed_count_difference"] < 0.005
        )
        scored["shape_pass"] = bool(
            scored["shape_rms_lattice"] < 0.5 * scored["shape_rms_parabola"]
        )
        out.append(
            dict(
                level=lev,
                seed=SEED,
                cells=cells,
                dimers=len(route),
                formation=commissioned,
                score=scored,
                observer=obs,
                native_transport=live,
                changed_transport=off,
            )
        )
    return dict(
        m=m, c_squared=c2, nominal_bin_width=2 * np.pi / (SAMPLES * DT), boards=out
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    target = parser.parse_args().output
    out = run()
    if target:
        target.write_text(json.dumps(out, indent=2) + "\n")
    else:
        verify(out, json.loads((HERE / "holdout.json").read_text()))
    for b in out["boards"]:
        print(b["level"], b["dimers"], b["score"])
