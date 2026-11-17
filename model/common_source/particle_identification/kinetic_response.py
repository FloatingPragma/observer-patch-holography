"""Inverse and frequency derivative of the complete finite CAR response.

This consumes the declared ring's full Lehmann response, not detector data.
It supplies a finite kinetic-normalization diagnostic, not a reconstruction
of a physical field, a momentum derivative, or a positive kinetic metric.
The frequency is conjugate to the programmed model parameter.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import response


SPECTRAL_TOLERANCE = 5e-12


def _validated_integral(packet):
    """Check the supplied spectral data; never substitute the CAR identity."""
    if type(packet) is not response.FiniteResponse:
        raise ValueError("a complete finite CAR response required")
    integrated = np.zeros((4, 4), complex)
    for poles, charge in ((packet.addition, 1), (packet.removal, -1)):
        if not isinstance(poles, tuple) or not poles:
            raise ValueError("both complete charged spectral channels required")
        for pole in poles:
            if type(pole) is not response.Pole:
                raise ValueError("finite charged Pole entries required")
            for value in (pole.frequency, pole.energy_cost):
                if (isinstance(value, (bool, np.bool_))
                        or not isinstance(value, (int, float, np.integer, np.floating))
                        or not math.isfinite(float(value))):
                    raise ValueError("finite real spectral frequencies and costs required")
            if (pole.charge_change != charge or pole.target_number != 2+charge
                    or type(pole.multiplicity) is not int or pole.multiplicity <= 0
                    or not math.isclose(pole.frequency, charge*pole.energy_cost,
                                        rel_tol=1e-12, abs_tol=SPECTRAL_TOLERANCE)):
                raise ValueError("charged spectral frequency convention is inconsistent")
            if np.ma.isMaskedArray(pole.residue):
                raise ValueError("complete finite spectral residues required")
            residue = np.asarray(pole.residue)
            if (residue.shape != (4, 4) or residue.dtype.kind not in "iufc"
                    or not np.isfinite(residue).all()):
                raise ValueError("complete finite spectral residues required")
            if np.max(abs(residue-residue.conj().T)) > SPECTRAL_TOLERANCE:
                raise ValueError("Hermitian spectral residues required")
            if np.linalg.eigvalsh(residue).min() < -SPECTRAL_TOLERANCE:
                raise ValueError("positive spectral residues required")
            integrated += residue
    if np.max(abs(integrated-np.eye(4))) > SPECTRAL_TOLERANCE:
        raise ValueError("actual integrated residue fails the canonical CAR sum rule")
    return integrated


def inverse_jet(packet, frequency, *, coordinates=None):
    """Return G, K=G^-1 and dK/dz from the same retained spectral data.

    For c'=A c, G'=A G A*, while K'=A^-* K A^-1 and its derivative
    transform by the same congruence. A is an invertible operator-coordinate
    change, not a tunable kinetic interaction. All four channels are retained.
    The actual pole residues must satisfy positivity and the canonical CAR
    sum within the declared numerical tolerance. These consistency checks
    do not prove that arbitrary supplied poles came from the stated Hamiltonian.
    """
    integrated = _validated_integral(packet)
    green = response.green_function(packet, frequency)
    z = complex(frequency)
    derivative = -sum((pole.residue/(z-pole.frequency)**2
                       for pole in (*packet.addition, *packet.removal)),
                      np.zeros((4, 4), complex))
    if coordinates is None:
        a = np.eye(4, dtype=complex)
    else:
        if np.ma.isMaskedArray(coordinates):
            raise ValueError("complete finite four-channel coordinates required")
        raw = np.asarray(coordinates)
        if raw.shape != (4, 4) or raw.dtype.kind not in "iufc" or not np.isfinite(raw).all():
            raise ValueError("complete finite four-channel coordinates required")
        a = np.asarray(raw, complex)
        singular = np.linalg.svd(a, compute_uv=False)
        if singular[-1] <= 4*np.finfo(float).eps*singular[0]:
            raise ValueError("numerically invertible four-channel coordinates required")
    green = a@green@a.conj().T
    derivative = a@derivative@a.conj().T
    inverse = np.linalg.solve(green, np.eye(4))
    kinetic_derivative = -inverse@derivative@inverse
    return {"green": green, "inverse": inverse,
            "inverse_frequency_derivative": kinetic_derivative,
            "CAR_integrated_residue": a@integrated@a.conj().T,
            "expected_canonical_CAR_Gram": a@a.conj().T,
            "canonical_CAR_sum_rule_residual": float(np.max(abs(integrated-np.eye(4)))),
            "spectral_validation_tolerance": SPECTRAL_TOLERANCE,
            "coordinate_condition_number": float(np.linalg.cond(a))}


def report():
    """Small numerical receipt; no finite-shot or outward error claim."""
    z = .31+.23j
    a = np.array([[2., .2j, 0., 0.], [0., .7, .1, 0.],
                  [0., 0., 1.3, -.15j], [.05, 0., 0., .9]], complex)
    ai = np.linalg.inv(a)
    rows = []
    for interaction in (False, True):
        packet = response.finite_response(interaction=interaction)
        base = inverse_jet(packet, z)
        changed = inverse_jet(packet, z, coordinates=a)
        other = inverse_jet(packet, -.19+.37j)
        step = 1e-5
        difference = (np.linalg.inv(response.green_function(packet, z+step))
                      -np.linalg.inv(response.green_function(packet, z-step)))/(2*step)
        rows.append({
            "interaction_enabled": interaction,
            "canonical_CAR_sum_rule_residual": base["canonical_CAR_sum_rule_residual"],
            "inverse_identity_residual": float(np.max(abs(base["inverse"]@base["green"]-np.eye(4)))),
            "inverse_congruence_residual": float(np.max(abs(changed["inverse"]-ai.conj().T@base["inverse"]@ai))),
            "derivative_congruence_residual": float(np.max(abs(changed["inverse_frequency_derivative"]-ai.conj().T@base["inverse_frequency_derivative"]@ai))),
            "derivative_centered_difference_residual": float(np.max(abs(difference-base["inverse_frequency_derivative"]))),
            "derivative_variation_between_two_frequencies": float(np.max(abs(other["inverse_frequency_derivative"]-base["inverse_frequency_derivative"]))),
        })
    source = Path(__file__).resolve()
    root = source.parents[3]
    sources = [source, source.with_name("response.py"), source.parent.parent/"mass_response/model.py"]
    return {"schema": "oph.finite_car_inverse_kinetic.v1",
            "source_sha256": {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            "rows": rows,
            "frequency_parameter": "complex frequency conjugate to programmed ring model parameter",
            "all_spectral_residues_retained": True,
            "actual_residue_positivity_and_CAR_sum_checked": True,
            "spectral_validation_tolerance": SPECTRAL_TOLERANCE,
            "input_is_complete_declared_Hamiltonian_response": True,
            "response_inferred_from_observer_records": False,
            "spatial_momentum_derivatives_supplied": False,
            "positive_physical_kinetic_metric_established": False,
            "physical_mass_or_species_identification": False,
            "finite_shot_or_outward_numerical_certificate": False,
            "empirical_mass_inputs": [],
            "scope": "Coordinate covariance of the inverse finite CAR response and its frequency derivative; interacting frequency dependence is retained, not absorbed into a fitted static mass."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(report(), indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
