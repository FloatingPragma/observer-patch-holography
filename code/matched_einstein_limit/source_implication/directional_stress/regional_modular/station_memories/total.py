"""Entropy of the whole closed flight on the number-preserving material control."""
import argparse
import json
from pathlib import Path

import sympy as sp
from mpmath import iv

from .. import station_region as region

c = region.c
HERE = Path(__file__).resolve().parent
SYMMETRIC = iv.mpf([-1, 1])


def run():
    iv.dps = 32
    reference, state_radius, odds = region.measure.source()
    directions, direction_radius = region.controls()["fixed_all_number_weights"]
    value = error = 0
    minimum = []
    tangent_size = []
    score_norm = 0
    for rho, d in zip(reference, directions, strict=True):
        k, arithmetic, spec = region.logarithm(rho, c.ri(0))
        centered = k - k.trace()*sp.eye(4)/4
        score_norm += region.norm(centered)**2
        lower = spec[0].a-c.ri(state_radius)
        minimum.append(lower)
        # The logarithm differential on an HS ball costs 1/lambda_min in HS.
        # Only the traceless part of the modular generator pairs with d.
        nd = region.norm(d)
        tangent_size.append(nd)
        value += c.ri(c.inner(centered, d))
        error += (nd+c.ri(direction_radius))*(c.ri(state_radius)/lower+2*arithmetic)
    error += iv.sqrt(score_norm)*c.ri(direction_radius)
    slope = value+SYMMETRIC*error
    assert slope.a > 0
    # Strong concavity of entropy along a normalized matrix chord:
    # -S'' >= ||D||_HS^2/lambda_max >= ||D||_HS^2.
    lo_norm = iv.sqrt(sum(v**2 for v in tangent_size))-c.ri(direction_radius)
    t = c.ri(sp.Rational(1,100))
    finite_upper = (t*slope-t*t*lo_norm**2/2).b
    # A blockwise sharper uniform eigenvalue ceiling follows from
    # lambda_max(R+sD)<=w_b lambda_max(rho_b)+|s| ||D_b||_HS.
    maxes = []
    for weight, rho, nd in zip((1-odds, odds), reference, tangent_size, strict=True):
        _, _, spec = region.logarithm(rho, c.ri(0))
        maxes.append(weight*(spec[-1].b+c.ri(state_radius))+t*(nd+c.ri(direction_radius)))
    strong = sum((nd-c.ri(direction_radius))**2/mx.b for nd,mx in zip(tangent_size,maxes,strict=True))
    sharper_upper = (t*slope-t*t*strong/2).b
    assert sharper_upper < finite_upper < 0
    return dict(control="fixed_all_number_weights", source_slope=str(slope),
                entropy_change_at_0p01_upper=str(sharper_upper),
                launch_entropy_constant="log(10)", fair_charge_and_launch_entropy_constant="log(20)",
                reference_conditional_eigenvalue_lower=[str(v) for v in minimum],
                tangent_HS_norm_lower=str(lo_norm.a))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    output=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        (HERE/"total.json").write_text(output)
    else:
        print(output,end="")
