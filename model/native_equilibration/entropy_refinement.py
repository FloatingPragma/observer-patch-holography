"""Uniform entropy gradients after finite preparation on native screens.

Use the existing C/E mixing theorem on each connected binary fixed-k graph:
D(rho_nR, I/d tensor rho_R) <= (d/2)(1-gamma)^n. Discard preparation
records and use the source marginal for entropy. No entropy statement on
the external reference, or on a retained/conditioned tape, is inferred.

If h>=ceil(log2 d), b=ceil(1/gamma_lower), and n=b(2h+q), each b-attempt
block contracts by at most 1/2: for 0<gamma_lower<1, Bernoulli gives
(1-gamma_lower)^b<=1/(1+b*gamma_lower)<=1/2; at gamma_lower=1 the
contraction is zero after one attempt. Thus theta=d*D(rho_n,I/d)<=2^(-q-1).
For every complete partition with counts c_i, M=d, and q_i=c_i/d,
TV(p,q)<=epsilon and min(q_i)>=1/d. Set S_count=H(p)+sum p_i log(c_i).
For positive-real count variations at fixed total d, zero-sum tangents u=p' and
w=c'/d, the exact derivative is
  S_count'=-sum u_i log(p_i/q_i)+sum w_i(p_i/q_i-1).
Thus |S_count'| <= ||u||_1 theta/(1-theta)+||w||_1 theta.
For all traceless Hermitian X, |dS_rho[X]| <= ||X||_1 theta/(1-theta).
The latter follows from (1-theta)I/d <= rho <= (1+theta)I/d and spectral
logarithms; it includes noncommuting tangents. Also log(d)-S(rho) and the
counting deficit are at most 2 theta^2/d, by log r<=r-1: the positive
and negative parts of a traceless density difference each have trace at
most epsilon, so Tr(Delta^2)<=2 epsilon^2; the classical bound is identical.

At refinement l, declare ell_l=2^-l and bound all tangent norms by 2^a_l.
Choose q_l=4+5l+a_l. The joint gradient is <=(3/32)ell_l^5, the quantum
gradient <=(1/16)ell_l^5, and either deficit <=ell_l^10/512. Consequently
each is o(ell_l^4), uniformly over all the stated tangents and partitions.
This controls only preparation's contribution. Other continuum errors,
physical radius/volume/area maps and the variation dictionary are separate.
The variation is infinitesimal: admissible two-sided intervals may shrink
with dimension. No common finite-amplitude variation interval is derived.
A subsequent common complete source/reference instrument preserves the
error <=theta/d, with all outcomes and timeouts kept and no horizon factor.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path

from . import construction as mixing
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "entropy_results.json"
ROOT = HERE.parents[1]
SOURCE_PATHS = (
    "model/native_equilibration/construction.py",
    "model/native_equilibration/results.json",
    "model/native_equilibration/entropy_refinement.py",
    "model/observer/QUANTUM_MIXING.md",
    "model/observer/quantum_mixing.py",
    "model/federated_observer/source.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/federation.py",
)


def _nonnegative_integer(value):
    if type(value) is not int or value < 0:
        raise ValueError("an original nonnegative integer is required")
    return value


def _schedule(log_dimension_upper, gap_lower, level, tangent_bits):
    h, level, a = map(_nonnegative_integer, (log_dimension_upper, level, tangent_bits))
    if type(gap_lower) is not F or not 0 < gap_lower <= 1:
        raise ValueError("an exact Fraction gap in (0,1] is required")
    q = 4+5*level+a
    block = (gap_lower.denominator+gap_lower.numerator-1)//gap_lower.numerator
    blocks = 2*h+q
    theta = F(1, 1 << (q+1))
    speed = 1 << a
    joint = speed*(theta/(1-theta)+theta)
    quantum = speed*theta/(1-theta)
    ell = F(1, 1 << level)
    if (block*gap_lower < 1 or joint > F(3, 32)*ell**5
            or quantum > F(1, 16)*ell**5 or 2*theta**2 > ell**10/512):
        raise ArithmeticError("integer entropy-refinement certificate failed")
    return {"level": level, "declared_refinement_coordinate": str(ell),
            "log2_dimension_upper": h, "tangent_norm_upper": str(speed),
            "precision_bits": q, "half_contraction_block": block,
            "blocks": blocks, "discarded_native_attempts": block*blocks,
            "gap_lower_per_attempt": str(gap_lower),
            "relative_density_error_upper": str(theta),
            "joint_counting_gradient_upper": str(joint),
            "quantum_state_gradient_upper": str(quantum),
            "entropy_deficit_upper": str(2*theta**2),
            "gradient_over_ell4_upper": str(F(3, 32)*ell),
            "entropy_deficit_over_ell4_upper": str(ell**6/512),
            "complete_future_output_error_times_dimension_upper": str(theta)}


def entropy_budget(dimension, gap_lower, level, *, tangent_bits=0):
    """All-partition/all-tangent certificate using the exact finite dimension.

    No power with exponent equal to the attempt count, floating logarithm,
    or enumeration of the d-dimensional state space is performed.
    """
    if type(dimension) is not int or dimension < 1:
        raise ValueError("a positive original integer sector dimension is required")
    row = _schedule((dimension-1).bit_length(), gap_lower, level, tangent_bits)
    if dimension == 1:
        # The sole density and sole complete partition have no nonzero tangent.
        for name in ("half_contraction_block", "blocks", "discarded_native_attempts"):
            row[name] = 0
        for name in ("relative_density_error_upper", "joint_counting_gradient_upper",
                     "quantum_state_gradient_upper", "entropy_deficit_upper",
                     "gradient_over_ell4_upper", "entropy_deficit_over_ell4_upper",
                     "complete_future_output_error_times_dimension_upper"):
            row[name] = "0"
    return row


def screen_family_budget(level, *, tangent_bits=0):
    """All-level source census, without constructing the exponentially big sector.

    The source constructor replaces each triangular face by four triangles
    with connected dual graph and retains adjacencies across old edges.
    Each twelve-port carrier is connected and native inter-seams connect
    adjacent carriers. Thus connectedness holds at all levels, with the
    source census V=240*4^l and M=630*4^l. Conservative D<=V-1 and
    binomial(V,k)<=2^V cover every supplied k. The native-attempt count is
    <=(21/16)*V^3*(2V+4+5l+a), using h=V and gamma_lower=2/(M*V*(V-1)).
    This is an allocation-free certificate, not execution or rapid mixing.
    """
    level = _nonnegative_integer(level)
    vertices, seams = 240*4**level, 630*4**level
    gap = F(2, seams*vertices*(vertices-1))
    return {"vertices": vertices, "seams": seams,
            "diameter_upper_bound": vertices-1,
            "dimension_bound": "binomial(V,k)<=2^V",
            **_schedule(vertices, gap, level, tangent_bits)}


def run():
    rows = []
    for level in (0, 1, 2):
        graph = mixing.source.federation(level)
        cert = mixing.graph_certificate(graph.ports, graph.seams)
        particles = graph.ports//3
        dimension = comb(graph.ports, particles)
        budget = entropy_budget(dimension, F(cert["quantum_gap_lower_per_attempt"]), level)
        rows.append({"particles": particles, "dimension": "binomial(V,k)",
                     **cert, **budget})
    return {"schema": "oph.native_equilibration.entropy_refinement.v1",
            "source_sha256": bound_pins({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                              for name in SOURCE_PATHS}),
            "retained_actual_screen_certificates": rows,
            "allocation_free_level_five_budget": screen_family_budget(5),
            "theorem": {
                "source": "connected binary fixed-k C/E repair, independent uniform seams and fair coins",
                "input": "arbitrary source density and untouched quantum reference within that same sector",
                "records": "discard all warm-up records; any common later complete instrument retains every outcome",
                "dimension": "d_l=binomial(V_l,k_l); the conserved sector is supplied",
                "schedule": "b_l=ceil(1/gamma_lower); h_l>=ceil(log2 d_l); n_l=b_l*(2*h_l+4+5*l+a_l)",
                "tangents": "sum p'=sum c'=0; ||p'||_1,||c'/d_l||_1,||X||_1<=2^a_l; X=X*, Tr X=0",
                "count_domain": "every complete nonempty partition, positive-real c variations at fixed total d_l",
                "quantum_domain": "source marginal on M_(d_l), not the possibly singular external-reference algebra",
                "gradient_bound": "|S_count'| <= (3/32)*ell_l^5; |dS_rho[X]| <= (1/16)*ell_l^5",
                "deficit_bound": "log(d_l)-S(rho_l), log(d_l)-S_count <= ell_l^10/512",
                "scale": "ell_l=2^-l is a declared refinement coordinate; both errors are o(ell_l^4)",
                "future_observation_bound": "complete source/reference output half-trace error <=2^(-5-5*l-a_l)/d_l",
                "cost": "discarded native attempts, with no duration or rapid-mixing claim"},
            "scope": {"exact_stationary_preparation_required": False,
                      "uniform_entropy_preparation_error_derived": True,
                      "growing_state_or_sector_enumeration_required": False,
                      "retained_or_conditioned_warmup_tape_covered": False,
                      "uniform_finite_amplitude_variation_interval_derived": False,
                      "cross_sector_weights_or_physical_vacuum_selected": False,
                      "physical_radius_area_volume_or_Einstein_branch_identified": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, sort_keys=True)+"\n"
    if args.check:
        if args.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit("native entropy-refinement receipt drift")
    else:
        with OUTPUT.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded)
    print("native entropy-refinement preparation: all checks passed")


if __name__ == "__main__":
    main()
