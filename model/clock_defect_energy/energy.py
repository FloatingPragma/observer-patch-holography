"""Test a repair-activity energy ansatz before assigning physical masses.

Activity is the complete native probability of a nonzero signed current per
uniform-seam attempt. References include uniform binary66 and the four
transported faithful extensions in the retained clock inventory. These choices
are declared independently of any target mass. They do not identify the
physical vacuum, proper time, action per repair or an energy observable.

The exact kinetic tests consume independently specified action coefficients.
Their positive control is the relativistic proper-time action, not a derived
action of the native adjacent-pair candidate. No momentum response or inertial
mass for that candidate is claimed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

from .contract import binding, source_bindings, check_receipt
from ..observer import native, record_recovery as recovery
from ..native_modular_reference import history_transport, history_clock
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent


def rational(value, *, positive=False):
    if type(value) is int:
        out = F(value)
    elif type(value) is F:
        n, d = value.numerator, value.denominator
        if type(n) is not int or type(d) is not int or d <= 0 or math.gcd(abs(n), d) != 1:
            raise ValueError("canonical original rational required")
        out = F(n, d)
    else:
        raise ValueError("exact original integer or Fraction required")
    if positive and out <= 0:
        raise ValueError("positive rational required")
    return out


def activity_table(law):
    """All native attempts; no conditioning on a successful movement."""
    law = recovery.candidate(law)
    basis, seams = native.states(binary=True), recovery.seams()
    rows = []
    for state in basis:
        movement = 0
        squared_current = 0
        for edge in seams:
            for coin in (0, 1):
                output, record = recovery.recorded_step(law, state, edge, coin)
                current = recovery.record_current(law, record)
                if (output != state) != (current != 0):
                    raise ArithmeticError("movement and actual current disagree")
                movement += abs(current)
                squared_current += current*current
        denominator = 2*len(seams)
        rows.append({"state": state, "attempts": F(1),
                     "absolute_current": F(movement, denominator),
                     "squared_current": F(squared_current, denominator)})
    return tuple(rows)


def maintenance_certificate():
    tables = {law: activity_table(law) for law in ("C", "E")}
    if tables["C"] != tables["E"]:
        raise ArithmeticError("the complete C/E population costs differ")
    rows = tables["C"]
    seams = set(recovery.seams())
    cluster = tuple(i for i, row in enumerate(rows)
                    if tuple(p for p, value in enumerate(row["state"]) if value) in seams)
    background = {name: sum((row[name] for row in rows), F())/len(rows)
                  for name in ("attempts", "absolute_current", "squared_current")}
    candidate = {name: sum((rows[i][name] for i in cluster), F())/len(cluster)
                 for name in background}
    excess = {name: candidate[name]-background[name] for name in background}
    if len(cluster) != 30 or excess["attempts"] != 0 or excess["absolute_current"] >= 0:
        raise ArithmeticError("the declared localization/activity verdict changed")
    return {"carrier_dimension": len(rows), "cluster_states": len(cluster),
            "state_indices": list(cluster), "reference": "uniform binary66",
            "seam_selector": "uniform over all thirty native seams",
            "parameter": "native attempted events, including zero current",
            "background_cost_per_attempt": background,
            "cluster_cost_per_attempt": candidate,
            "signed_excess_cost_per_attempt": excess,
            "positive_action_scale_times_this_excess_is_positive": False,
            "actual_reference_vacuum_or_particle_identified": False,
            "rows": rows}


def source_cycle(law):
    """A closed LOAD path with positive activity; the exported record differs."""
    law = recovery.candidate(law)
    edge = recovery.seams()[0]
    spectator = next(i for i in range(12) if i not in edge)
    initial = tuple(int(i in (edge[0], spectator)) for i in range(12))
    state, tape, activity = initial, [], 0
    for coin in (0, 1):
        state, record = recovery.recorded_step(law, state, edge, coin)
        tape.append(record)
        activity += abs(recovery.record_current(law, record))
    if state != initial or activity != 2:
        raise ArithmeticError("native closed-load cycle changed")
    return {"initial_and_final_load": initial, "complete_records": tuple(tape),
            "conditional_coin_word_probability": F(1, 4),
            "seam_intervention": edge, "unsigned_activity": activity,
            "signed_net_current": 0,
            "source_state_function_energy_increment": 0,
            "entire_source_plus_archive_returns_to_initial_state": False,
            "conclusion": "Unsigned movement counts cannot equal increments of a source-only state energy on every branch. This does not exclude a history action or source-plus-environment energy."}


def transported_reference_activity():
    """Use predecessor target reference extensions, not posterior corners.

    Conditioning on adjacency defines the offline candidate ensemble. This
    calculation introduces no adjacency measurement or controller input.
    """
    rows = activity_table("C")
    edges = set(recovery.seams())
    cluster = tuple(i for i, row in enumerate(rows)
                    if tuple(p for p, n in enumerate(row["state"]) if n) in edges)
    references = {}
    for law in ("C", "E"):
        for branch in history_clock.branches(law):
            permutation = history_transport.port_permutation(law, branch)
            core = tuple(sorted(permutation[p] for p in (0, 1)))
            sigma = history_transport.extended_reference(law, branch, history_transport.reference())
            if core in references and references[core] != sigma:
                raise ArithmeticError("same core has inconsistent transported reference")
            references[core] = sigma
    results = []
    for core, sigma in sorted(references.items()):
        mass = sum((sigma[i] for i in cluster), F())
        activity = sum((p*row["absolute_current"] for p, row in zip(sigma, rows, strict=True)), F())
        conditional = sum((sigma[i]*rows[i]["absolute_current"] for i in cluster), F())/mass
        results.append({"core": core, "reference_cluster_mass": mass,
                        "reference_activity_per_attempt": activity,
                        "adjacency_conditioned_activity_per_attempt": conditional,
                        "signed_excess_activity_per_attempt": conditional-activity})
    return {"reference_kind": "faithful full-shell preparation-target extensions, not actual posterior corners",
            "seam_selector": "uniform over thirty seams",
            "adjacency_conditioning_is_executed_measurement": False,
            "rows": results}


def modular_cost_witness():
    """Same instantaneous costs, different exact preparation-target weights."""
    table = activity_table("C")
    sigma = history_transport.reference("local")
    fields = ("attempts", "absolute_current", "squared_current")
    for i, left in enumerate(table):
        for j, right in enumerate(table[:i]):
            if all(left[name] == right[name] for name in fields) and sigma[i] != sigma[j]:
                ratio = sigma[j]/sigma[i]
                return {"source_indices": [i, j], "source_states": [left["state"], right["state"]],
                        "equal_native_costs": {name: left[name] for name in fields},
                        "reference_weights": [sigma[i], sigma[j]],
                        "reference_weight_ratio_j_over_i": ratio,
                        "modular_energy_difference_K_i_minus_K_j": f"log({ratio})",
                        "any_function_of_this_instantaneous_cost_vector_equals_K_up_to_shift": False,
                        "modular_energy_is_physical_rest_energy_claimed": False}
    raise ArithmeticError("no unequal-reference equal-cost witness")


def kinetic_consistency(rest_energy, quadratic, wave_speed, quartic=None):
    """Compare L=-E0+A*v^2+B*v^4 with the relativistic coefficients.

    Inputs must have one common clock and coordinate dictionary. Inertia is
    extracted from the supplied kinetic coefficient (2A), independently of
    E0/c^2. Passing these two coefficients is necessary, not sufficient, for
    a complete relativistic dispersion or source-derived energy law.
    """
    energy, a, c = (rational(value, positive=True)
                    for value in (rest_energy, quadratic, wave_speed))
    mass = 2*a
    result = {"rest_energy": energy, "kinetic_mass": mass, "wave_speed": c,
              "rest_to_inertia_ratio": energy/(mass*c*c),
              "rest_inertia_match": energy == mass*c*c}
    if quartic is not None:
        b = rational(quartic)
        result.update({"quartic_coefficient": b,
                       "quartic_to_relativistic_ratio": 8*b*c**4/energy,
                       "quartic_match": 8*b*c**4 == energy})
    return result


def proper_time_coefficients(rest_energy, wave_speed):
    """Declared positive benchmark, not coefficients recovered from repair."""
    energy, c = (rational(value, positive=True) for value in (rest_energy, wave_speed))
    return energy/(2*c*c), energy/(8*c**4)


def scale_coefficients(rest_energy, quadratic, wave_speed, quartic, *,
                       energy_unit, length_unit, time_unit):
    """Convert all action coefficients together; a scalar unit choice cannot fix a failed ratio."""
    e, a, c, b = map(rational, (rest_energy, quadratic, wave_speed, quartic))
    eu, lu, tu = (rational(v, positive=True) for v in (energy_unit, length_unit, time_unit))
    return e*eu, a*eu*tu**2/lu**2, c*lu/tu, b*eu*tu**4/lu**4


def action_controls():
    e, c = F(7, 3), F(5, 2)
    a, b = proper_time_coefficients(e, c)
    controls = {"proper_time_action": kinetic_consistency(e, a, c, b),
                "same_rest_energy_twice_inertia": kinetic_consistency(e, 2*a, c, b),
                "rest_inertia_equality_but_wrong_quartic": kinetic_consistency(e, a, c, 0)}
    changed = scale_coefficients(e, 2*a, c, b, energy_unit=F(11, 7),
                                 length_unit=F(3, 5), time_unit=F(17, 13))
    scaled = kinetic_consistency(changed[0], changed[1], changed[2], changed[3])
    if scaled["rest_to_inertia_ratio"] != controls["same_rest_energy_twice_inertia"]["rest_to_inertia_ratio"]:
        raise ArithmeticError("dimensionless failure changed with units")
    return {"parameter": "declared model coordinate time and velocity",
            "action": "S=-E0*integral(sqrt(1-v^2/c^2) dt)",
            "repair_action_hypothesis": "If S=-a*N with positive action a per event and N=r*integral(d_tau) for a constant proper-time rate r, then E0=a*r. The native action and proper-time identifications remain to be derived.",
            "controls": controls, "rescaled_wrong_inertia": scaled,
            "native_particle_kinetic_response_executed": False,
            "mass_or_action_scale_fitted": False}


def _json(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def report():
    sources = ("model/clock_defect_energy/energy.py", "model/clock_defect_energy/test_energy.py",
               "model/clock_defect_energy/contract.py", "model/observer/native.py",
               "model/observer/record_recovery.py", "model/native_modular_reference/history_transport.py",
               "model/native_modular_reference/local_preparation.py",
               "model/native_modular_reference/history_clock.py", "tools/strict_json.py")
    return _json({"schema": "oph.clock_defect_energy.energy.v1", "campaign": binding(),
                  "source_sha256": bound_pins(source_bindings(sources)),
                  "maintenance": maintenance_certificate(),
                  "transported_reference_activity": transported_reference_activity(),
                  "source_cycles": {law: source_cycle(law) for law in ("C", "E")},
                  "modular_cost_witness": modular_cost_witness(),
                  "kinetic_action_controls": action_controls(),
                  "scope": {"integrated_gate_accepted": False, "physical_mass_derived": False,
                            "arbitrary_repair_energy_maps_excluded": False,
                            "gravitational_stress_constructed": False,
                            "physical_comparison_count_added": 0}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    value = report()
    path = HERE/"energy_results.json"
    if args.check:
        check_receipt(path, value)
    if args.write:
        path.write_text(json.dumps(value, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"signed_excess": value["maintenance"]["signed_excess_cost_per_attempt"],
                      "same_cost_weight_ratio": value["modular_cost_witness"]["reference_weight_ratio_j_over_i"]}, indent=2))


if __name__ == "__main__":
    main()
