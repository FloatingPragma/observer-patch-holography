"""Test exact fixed-survey entropy extremality at the measured dynamical fixed points."""

import argparse
import json
from hashlib import sha256
from os.path import relpath
from pathlib import Path

import numpy as np
import sympy as sp
from matched_einstein_limit.source_implication.directional_stress.connection import (
    force_entropy,
)

from . import persistent_density as d

HERE, c, reg, lab = d.HERE, d.c, d.reg, d.lab


def score_direction(data, score):
    native = force_entropy.c
    basis = [native.prior.sym_matrix(x) for x in native.iso.basis(True)]
    force_rows = (
        sp.Matrix(data["force_gradient_numerators"])
        / data["force_gradient_denominator"]
    )
    d, correction, _, _ = force_entropy.corrected_kernel_direction(
        *force_entropy.range_gradients(basis),
        force_rows,
        sp.Rational(data["force_gradient_error_per_row"]),
        score=score,
    )
    return np.asarray(d, dtype=float).ravel(), float(correction.b) + 1e-12


def run():
    source = HERE / "persistent_density.json"
    kernel = Path(force_entropy.__file__).with_name("force_entropy.json")
    bank = json.loads(source.read_bytes())
    data = json.loads(kernel.read_bytes())
    vec = np.array(data["direction_center_even_then_odd"])
    err = (
        float(data["true_kernel_correction_HS_upper"].strip("[]").split(",")[-1])
        + 1e-12
    )
    basis = c.dose.heads.measure.joint_basis()
    frozen = (vec, err)
    rows = {}
    for name, z in bank["schedules"].items():
        rows[name] = {}
        if "terminal_heads" not in z:
            continue
        for arm, tape in z["terminal_heads"].items():
            x, e = reg.src.observed(tape, lab.SHOTS, False)
            radius = np.linalg.norm(2 * e) + 2e-11
            recipe = c.station.preparation(x)
            margin = min(recipe["min_kept"], recipe["min_dephased"])
            floor = np.linalg.eigvalsh(x).min() - radius - 1e-11
            K = c.station.pair(basis, -c.station.log_blocks(x))
            # On each marginal-preserving HS unit tangent, only the 18
            # connected components contribute to the trace pairing.
            vec, err = frozen
            norm = np.linalg.norm(vec) + err + 1e-12
            slope = vec @ K
            error = norm * np.sqrt(8) * radius / floor + err * np.linalg.norm(K) + 2e-11
            low = abs(slope) - error
            record = {
                "reference_HS_radius": float(radius),
                "native_weight_margin": margin,
                "entropy_slope_in_frozen_direction": [
                    float(slope - error),
                    float(slope + error),
                ],
            }
            if low <= 0:
                direction, correction = score_direction(data, K)
                norm = np.linalg.norm(direction) + correction + 1e-12
                slope = direction @ K
                error = (
                    norm * np.sqrt(8) * radius / floor
                    + correction * np.linalg.norm(K)
                    + 2e-11
                )
                low = abs(slope) - error
                record["followup"] = {
                    "center_even_then_odd": direction.tolist(),
                    "direction_HS_correction": correction,
                    "entropy_slope": [float(slope - error), float(slope + error)],
                }
            if low > 0:
                for j in range(1, 60):
                    t = 2.0**-j
                    # Correct either rail decomposition using the product-Pauli
                    # expansion. Per atom at most 4 times the HS displacement
                    # is needed, including correction to the observed reference.
                    if 4 * (radius + t * norm + 1e-11) >= margin or t * norm >= floor:
                        continue
                    gain = t * low - t * t * norm**2 / (2 * (floor - t * norm)) - 1e-15
                    if gain > 0:
                        record.update(
                            sign=1 if slope > 0 else -1,
                            step=t,
                            minimum_finite_entropy_gain=float(gain),
                        )
                        break
            rows[name][arm] = record
    return {
        "works_at": rows,
        "true_HS_correction_upper": err,
        "inputs": {
            relpath(p, HERE): sha256(p.read_bytes()).hexdigest()
            for p in (
                source,
                kernel,
                HERE / "STATIONARY_SOURCE_TEST.md",
                HERE / "STATIONARY_SECTION_PROTOCOL.md",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    r = run()
    path = HERE / "stationary_source_test.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(
        {
            k: {a: v.get("minimum_finite_entropy_gain") for a, v in t.items()}
            for k, t in r["works_at"].items()
        }
    )
