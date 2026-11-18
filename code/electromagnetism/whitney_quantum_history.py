"""A full-configuration potential-phase trial history with exact norm bounds.

The certificate uses global polynomial envelopes and exact Gaussian moments,
not sampled quadrature or Monte Carlo errors. Its very short time interval is
conservative and dimensionless. It is not an exact Hamiltonian trajectory or
an ordinary-physics benchmark. Existing geometry and quantization are inputs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import numpy as np

import whitney_interacting_quantum as quantum

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent/"runtime/whitney_quantum_history_receipt.json"
SCHEMA = "oph.whitney_quantum_trial_history.v1"
SCOPE = "FULL_56D_NEUTRAL_POTENTIAL_PHASE_TRIAL__GLOBAL_NORM_BOUND__DECLARED_INPUTS"
PARAMETERS = {"sigma": "1/2", "charge": "1/4", "mass_squared": "1/2",
              "quartic": "1/4", "hbar": "1", "target_norm_error": "1/10"}
PIN_PATHS = (
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "code/electromagnetism/whitney_interacting_quantum.py",
    "code/electromagnetism/verify_whitney_quantum_state.py",
    "code/electromagnetism/runtime/whitney_quantum_state_receipt.json",
    "code/electromagnetism/whitney_quantum_history.py",
    "code/electromagnetism/verify_whitney_quantum_history.py",
    "code/electromagnetism/test_whitney_quantum_history.py",
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode("ascii")


def add(*polys):
    result = {}
    for poly in polys:
        for powers, coefficient in poly.items():
            result[powers] = result.get(powers, Q(0))+coefficient
    return {key: value for key, value in result.items() if value}


def multiply(left, right):
    result = {}
    for (i, j), a in left.items():
        for (k, l), b in right.items():
            key = (i+k, j+l)
            result[key] = result.get(key, Q(0))+a*b
    return result


def gaussian_radial_moment(dimension, power, variance=Q(1, 4)):
    result = variance**power
    for j in range(power):
        result *= dimension+2*j
    return result


def gaussian_expectation(poly):
    return sum((coefficient*gaussian_radial_moment(30, i)*gaussian_radial_moment(26, j)
                for (i, j), coefficient in poly.items()), Q(0))


def polynomial_records(poly):
    return [{"powers": [i, j], "coefficient": str(value)} for (i, j), value in sorted(poly.items())]


def envelopes():
    inverse = {(0, 0): Q(209), (0, 1): Q(46656)}
    kinetic_factor = {(0, 0): Q(2*1896**2), (0, 1): Q(2*11286**2+16), (1, 0): Q(16)}
    kinetic = {key: value/4 for key, value in multiply(inverse, kinetic_factor).items()}
    potential = {(1, 0): Q(3456), (0, 1): Q(585), (1, 1): Q(729, 4), (0, 2): Q(9, 4)}
    gradient = {(1, 0): Q(4*6912**2), (0, 2): Q(2313**2),
        (2, 2): Q(4*1152**2), (0, 4): Q(81, 4), (0, 1): Q(3*1170**2),
        (2, 1): Q(3*1152**2), (0, 3): Q(243)}
    force = multiply(inverse, gradient)
    return {"inverse_metric": inverse, "initial_kinetic": kinetic,
            "potential": potential, "potential_gradient_squared": gradient,
            "force_moment": force}


def energy_bounds():
    polys = envelopes()
    kinetic = gaussian_expectation(polys["initial_kinetic"])
    potential = gaussian_expectation(polys["potential"])
    return {"kinetic_upper": kinetic, "potential_upper": potential,
            "total_energy_upper": kinetic+potential,
            "force_moment_upper": gaussian_expectation(polys["force_moment"])}


def squared_error_bound(t, bounds):
    t = abs(Q(t))
    return (bounds["total_energy_upper"]+bounds["kinetic_upper"])*t+Q(1, 6)*bounds["force_moment_upper"]*t**3


def selected_horizon(bounds):
    power = 0
    target = Q(PARAMETERS["target_norm_error"])**2
    while squared_error_bound(Q(1, 2**power), bounds) > target:
        power += 1
    return power, Q(1, 2**power)


def qs_add(a, b):
    return (a[0]+b[0], a[1]+b[1])


def qs_mul(a, b):
    return (a[0]*b[0]+5*a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def phase_configurations():
    volume, apex_gradient_squared = (Q(10), Q(10, 3)), (Q(21, 2), Q(-9, 2))
    result = []
    for name, center, boundary in (("origin", Q(0), Q(0)),
                                    ("center-real", Q(1, 2), Q(0)),
                                    ("center-imaginary", Q(1, 2), Q(0)),
                                    ("center-boundary-relative-phase", Q(1, 2), Q(1, 4))):
        q = [Q(0)]*56
        q[43 if name == "center-imaginary" else 30] = center
        for node in range(1, 13):
            q[43+node] = boundary
        # Psi=center*lambda0+i*boundary*(1-lambda0), apart from a global
        # phase in the center-imaginary case. Exact Beta(1,3) moments.
        gradient = qs_mul(volume, (apex_gradient_squared[0]*(center**2+boundary**2),
                                    apex_gradient_squared[1]*(center**2+boundary**2)))
        mass = Q(1, 2)*(center**2/Q(10)+Q(3, 5)*boundary**2)
        quartic = Q(1, 8)*(center**4/Q(35)+Q(2, 35)*center**2*boundary**2+Q(3, 7)*boundary**4)
        potential = qs_add(gradient, (volume[0]*(mass+quartic), volume[1]*(mass+quartic)))
        result.append({"id": name, "coordinates_exact": [str(x) for x in q],
                       "potential_in_Qsqrt5": [str(x) for x in potential]})
    return result


def trial_log_state(q, t, mesh=None):
    """Numerical full56D evaluator; excluded from the exact error certificate.

    Return approximate curved-measure log magnitude and potential phase.
    The certificate instead uses global analytic bounds and exact sample
    phases; these numerical coefficients carry no interval error enclosure.
    """
    if type(t) not in (int, float, Q):
        raise ValueError("finite real trial time required")
    try:
        numeric_time = float(t)
    except (ValueError, OverflowError) as error:
        raise ValueError("finite real trial time required") from error
    if not np.isfinite(numeric_time):
        raise ValueError("finite real trial time required")
    mesh = quantum.geometry(4) if mesh is None else mesh
    log_magnitude = quantum.gaussian_state_log_amplitude(q, sigma=.5, charge=.25, mesh=mesh)
    q = np.asarray(q, dtype=float)
    _, potential = quantum.coefficients(mesh.slice[:42, :30]@q[:30], q[30:43]+1j*q[43:],
        charge=.25, mass_squared=.5, quartic=.25, mesh=mesh)
    return {"log_magnitude_numeric": log_magnitude, "phase_numeric": -numeric_time*potential,
            "scope": "quadrature approximation; not used to certify the norm error"}


def build():
    bounds = energy_bounds()
    power, horizon = selected_horizon(bounds)
    configurations = phase_configurations()
    history = []
    for step in range(5):
        time = horizon*Q(step, 4)
        history.append({"time": str(time), "squared_norm_error_upper": str(squared_error_bound(time, bounds)),
            "phase_samples": [{"configuration": row["id"],
                "angle_in_Qsqrt5": [str(-time*Q(value)) for value in row["potential_in_Qsqrt5"]],
                "interpretation": "trial_to_initial_wavefunction_ratio=exp(i*angle)"} for row in configurations]})
    return {"schema": SCHEMA, "scope": SCOPE,
        "run_id": "whitney-full56d-potential-phase-certified-trial-v1",
        "source_pins": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in PIN_PATHS},
        "parameters": PARAMETERS,
        "state": {"dimension": 56, "trial_formula": "v(t,q)=exp(-i*t*V(q)/hbar)*f_sigma(q)",
            "exact_evolution": "u(t)=exp(-i*t*H/hbar)*f_sigma",
            "initial_state_parent": "code/electromagnetism/runtime/whitney_quantum_state_receipt.json",
            "measure": "dmu=rho(q)dq; rho=sqrt(det(gamma))",
            "probability_law": "|v(t,q)|^2*dmu=N(0,sigma^2 I56); independent of t",
            "neutrality": "V and f_sigma are invariant under the residual global U(1)",
            "trial_history_computed": True, "exact_Hamiltonian_history_computed": False,
            "rigorous_global_norm_bound": True, "configuration_density_moves": False,
            "classical_five_coordinate_restriction": False, "observer_history": False,
            "physical_state_preparation": False, "empirical_comparison": False,
            "continuum_QFT": False, "ordinary_physics_benchmark": False},
        "geometry_envelope": {"volume_upper": "18", "tetrahedron_volume_lower": "5/6",
            "edge_length_squared_upper": "4", "barycentric_gradient_squared_upper": "4",
            "minimum_node_multiplicity": 5, "minimum_edge_multiplicity": 2,
            "dressed_mass_lower": "1/128", "Maxwell_mass_lower": "1/81"},
        "coefficient_bounds": {"variables": "X=|a|^2,Y=|psi|^2 in orthonormal Coulomb coordinates",
            "log_density_gradient_upper": "1896+11286*sqrt(Y)",
            "polynomials": {name: polynomial_records(poly) for name, poly in envelopes().items()}},
        "integration": {"method": "exact independent scaled chi-square moments, dimensions30 and26",
            "radial_moment": "E[X^k]=sigma^(2k)*product(d+2j,j=0..k-1)",
            "global_envelopes": True, "Gaussian_tail_truncation": False,
            "Monte_Carlo_used": False, "quadrature_used_for_bound": False,
            "floating_point_used_for_bound": False},
        "bounds": {name: str(value) for name, value in bounds.items()},
        "error_certificate": {"norm": "L2(Q,dvol_gamma)",
            "bound_formula": "||u(t)-v(t)||^2 <= ((Ebar+Kbar)*abs(t)+Bbar*abs(t)^3/6)/hbar",
            "trivial_norm_error_upper": "2", "target_norm_error": "1/10",
            "dyadic_horizon_power": power, "horizon": str(horizon),
            "horizon_squared_error_upper": str(squared_error_bound(horizon, bounds)),
            "previous_dyadic_squared_error_upper": str(squared_error_bound(2*horizon, bounds)),
            "coverage": "every real t with |t|<=horizon, not only stored samples",
            "time_units": "declared dimensionless Hamiltonian time; no physical clock calibration",
            "usefulness": "very conservative tiny interval; not an ordinary-physics benchmark",
            "proof_source": "paper/tex_fragments/WHITNEY_QUANTUM_HISTORY.tex",
            "proof_label": "thm:whitney-quantum-trial-history"},
        "phase_configurations": configurations, "history": history}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    packet = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({"receipt": str(args.output), "horizon": packet["error_certificate"]["horizon"],
                      "bounds": packet["bounds"]}, sort_keys=True))
