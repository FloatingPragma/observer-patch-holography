"""A selected neutral Gaussian initial state of the interacting finite action.

Exact rational moment coefficients are distinct from illustrative numerical
log-amplitude samples. This is not a quantum time history, an observer
preparation, or a quantum state attached to the classical charged trajectory.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import numpy as np

import whitney_interacting_quantum as quantum

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/"runtime/whitney_quantum_state_receipt.json"
SCHEMA = "oph.whitney_quantum_initial_state.v1"
SCOPE = "SELECTED_NEUTRAL_GAUSSIAN_INITIAL_STATE__EXACT_MATTER_MOMENTS__NO_QUANTUM_TIME_HISTORY"
PIN_PATHS = (
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "code/electromagnetism/verify_cone_whitney_bridge.py",
    "code/electromagnetism/whitney_interacting_quantum.py",
    "code/electromagnetism/whitney_quantum_state.py",
    "code/electromagnetism/verify_whitney_quantum_state.py",
    "code/electromagnetism/test_whitney_quantum_state.py",
)
PARAMETERS = {"sigma": "1/2", "charge": "1/4", "mass_squared": "1/2", "quartic": "1/4", "hbar": "1"}


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode("ascii")


def state_definition():
    return {
        "real_configuration_dimension": 56, "radiative_coordinates": 30,
        "complex_matter_coordinates": 13,
        "coordinate_order": "30 Euclidean-orthonormal Coulomb coordinates,13 Re(psi),13 Im(psi)",
        "measure": "dmu=rho(q)dq; rho(q)=sqrt(det(gamma(q)))",
        "formula": "f_sigma(q)=rho(q)^(-1/2)*(2*pi*sigma^2)^(-14)*exp(-||q||^2/(4*sigma^2))",
        "rho_power": "-1/2", "gaussian_prefactor_power": "-14",
        "gaussian_exponent_denominator_factor": "4",
        "probability_density": "|f_sigma(q)|^2*dmu=(2*pi*sigma^2)^(-28)*exp(-||q||^2/(2*sigma^2))*dq",
        "selection": "supplied centered Gaussian width; not selected by an observer preparation",
        "action_attachment": "same supplied interacting finite Whitney action; not a wavefunction assigned to the recorded classical trajectory",
        "neutrality": "global U(1)-invariant; all real and imaginary matter coordinates have the same centered Gaussian variance",
        "normalization_scope": "analytic Gaussian/half-density identity; finite coefficient checks do not prove operator-domain claims",
        "proof_source": "paper/tex_fragments/WHITNEY_INTERACTING_QUANTUM.tex",
        "proof_label": "prop:whitney-interacting-gaussian-state",
        "analytic_domain_claim": "normalized vector in the neutral Hamiltonian operator domain, by the cited analytic proof",
        "initial_state_constructed": True, "initial_observables_provided": True,
        "quantum_time_history": False, "observer_history": False,
        "physical_comparison": False, "empirical_prediction": False,
    }


def exact_moments():
    sigma = Q(PARAMETERS["sigma"])
    # Unit volume and unit trace select algebraic coefficients, not a claim
    # that this cone has unit volume or unit transverse stiffness trace.
    coefficients = quantum.gaussian_initial_moments(sigma, Q(1), Q(1))
    l2, l4 = coefficients["matter_l2"], coefficients["matter_l4"]
    factors = {"matter_l2": l2, "matter_l4": l4,
               "mass_potential": Q(PARAMETERS["mass_squared"])*l2,
               "quartic_potential": Q(PARAMETERS["quartic"])*l4/2}
    volume = [Q(10), Q(10, 3)]
    return {
        "normalized_gaussian_norm_squared": "1",
        "one_real_coordinate_variance": str(sigma*sigma),
        "complex_gaussian_wick_factors": {"second": "2", "fourth": "8"},
        "simplex_moments_divided_by_volume": {"lambda_i_squared": "1/10",
            "lambda_i_fourth": "1/35", "lambda_i_squared_lambda_j_squared": "1/210",
            "sum_lambda_squared": "2/5", "sum_lambda_squared_squared": "6/35"},
        "volume_in_Qsqrt5": [str(x) for x in volume],
        "coefficients_per_volume": {name: str(value) for name, value in factors.items()},
        "expectations_in_Qsqrt5": {name: [str(value*x) for x in volume] for name, value in factors.items()},
        "interpretation": {"matter_l2": "E[integral_Omega |Psi_a|^2 dx]",
            "matter_l4": "E[integral_Omega |Psi_a|^4 dx]",
            "mass_potential": "E[m^2 integral_Omega |Psi_a|^2 dx]",
            "quartic_potential": "E[(g/2) integral_Omega |Psi_a|^4 dx]"},
        "derivation": "condition on a; independent circular complex Gaussian nodal matter removes every unit dressing phase; exact simplex moments then integrate the Wick factors",
        "not_evaluated": ["kinetic Hamiltonian expectation", "covariant-gradient expectation", "magnetic-energy value", "full energy expectation", "time-dependent observables"],
    }


def amplitude_samples():
    base = quantum.geometry(4)
    result = []
    for name, slot in (("origin", None), ("center-real", 30), ("center-imaginary", 43)):
        coordinates = [Q(0)]*56
        if slot is not None:
            coordinates[slot] = Q(1, 2)
        q = np.array([float(value) for value in coordinates])
        samples = []
        for order in (4, 5, 6):
            mesh = replace(quantum.geometry(order), slice=base.slice, mean_zero=base.mean_zero)
            reduction = quantum.reduced_kinetic(np.zeros(42), q[30:43]+1j*q[43:],
                charge=.25, mass_squared=.5, quartic=.25, mesh=mesh, include_potential=False)
            logdet = reduction.logdet()
            samples.append({"quadrature_order": order,
                "log_rho_numeric": float(logdet/2),
                "log_f_numeric": quantum.gaussian_state_log_amplitude(q, sigma=.5, charge=.25, mesh=mesh)})
        result.append({"id": name, "coordinates_exact": [str(value) for value in coordinates],
            "log_g_numeric": quantum.gaussian_half_density_log(q, sigma=.5), "evaluations": samples})
    return {"scope": "illustrative floating-point amplitudes from element quadrature and Schur determinants; no certified amplitude enclosure",
            "basis": "same Euclidean-orthonormal Coulomb frame at orders4,5,6; all30 radiative sample coordinates are zero",
            "samples": result}


def build():
    packet = {"schema": SCHEMA, "scope": SCOPE,
        "run_id": "whitney-neutral-gaussian-initial-state-2026-09-06-v1",
        "source_pins": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in PIN_PATHS},
        "parameters": PARAMETERS, "state": state_definition(), "exact_observables": exact_moments(),
        "numerical_amplitudes": amplitude_samples()}
    return json.loads(canonical(packet))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    packet = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({"receipt": str(args.output), "bytes": args.output.stat().st_size,
        "matter_coefficients": packet["exact_observables"]["coefficients_per_volume"]}, sort_keys=True))
