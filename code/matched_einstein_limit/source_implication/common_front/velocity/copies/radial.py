"""Mean finite radial flight on repeated preparations versus the velocity-dual distance."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ....action_response import measure as action
from ....action_response.four_ticks import raw_effect, volume_ratio
from ....directional_stress.traction.certificate import four_tick_columns
from .. import measure as velocity
from .measure import HERE


def run():
    iv.dps = 32
    board = velocity.formation.federation(1)
    route = velocity.formation.region(board)
    start = 4 * route[0] + np.arange(4)
    # Each trial selects a uniformly random shell 1, 2, 3 or 4 and scores
    # whether its sharp terminal cell lies on or outside that shell.
    weights = np.repeat(
        np.minimum(velocity.cell_distances(board)[route[0] // 12], 4), 48
    ).astype(object)
    native = velocity.front.native
    coeff = native.prior.weights()[7]
    source = action.front.inputs()
    out = {}
    rho = volume_ratio(source)
    ratios = {}
    for q in (-1, 1):
        effects = []
        for odd in (False, True):
            h = native.generator(board, route, q, odd)
            a, b = four_tick_columns(h, np.arange(h.shape[0]), start, coeff)
            xa, xb = weights[:, None] * a, weights[:, None] * b
            effect = (
                sp.Matrix(a.T @ xa + b.T @ xb) + sp.I * sp.Matrix(a.T @ xb - b.T @ xa)
            ) / (4 * native.prior.FRAC**2)
            assert effect == effect.H
            effects.append(effect)
        # Even a componentwise 1e-19 bound implies a Frobenius error <2e-17
        # for this 3840x4 compression. The squared effect error is <5e-17.
        prob = raw_effect(effects, source[0], sp.Rational("6e-7"), sp.Rational("1e-16"))
        allowance = iv.mpf(["-2e-13", "2e-13"])
        jets = {
            k: j + allowance for k, j in action.source_jets(effects, source).items()
        }
        for slope in jets.values():
            assert slope.b < 0
        raw_o, raw_p = [
            raw_effect(effects, source[1][name], sp.Rational(err), sp.Rational("2e-16"))
            for name, err in (("ordinary", "12e-7"), ("unpair", "18e-7"))
        ]
        direct = jets["pairing"] / jets["ordinary"]
        # Cancel the shared compensating rotation in the matched source jets.
        common = 1 / (rho + (raw_o - rho * raw_p) / jets["pairing"])
        ratio = iv.mpf([max(direct.a, common.a), min(direct.b, common.b)])
        ratios[q] = ratio
        out[str(q)] = dict(
            collar_over_ordinary=str(ratio),
            probability=str(prob),
            probability_slopes={k: str(v) for k, v in jets.items()},
            expected_capped_velocity_distance=str(8 * prob),
            expected_distance_slopes={k: str(8 * v) for k, v in jets.items()},
        )
    assert ratios[-1].b < ratios[1].a
    error = (ratios[1].a - ratios[-1].b) / (ratios[1].a + ratios[-1].b)
    return dict(
        dwell=4,
        cap_in_cell_hops=4,
        sharp_launch=int(route[0]),
        charges=out,
        common_response_relative_error_at_least=str(error),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "radial.json").write_text(out)
    else:
        print(out, end="")
