"""Exact reference/law compatibility for the native Dirichlet-to-Higgs route.

For the actual faithful core-biased reference, stationarity of an unread
fair native seam mixture excludes every core-to-outer seam with positive
weight. This is an entire nonnegative weight cone, not one fitted schedule.
Its maximal22-seam graph includes the core-internal edge omitted by the
earlier completion. Full complement elimination attains family rank two,
but a conserved core-charge direction prevents rank three.

The result constrains the stated passive stationary Dirichlet attachment.
It is not a no-go for active recovery, a different joint reference, feedback
that changes the admitted seam set, or another coherent mass interaction.
No measured mass or coupling is used. The parameter is one native proposal;
its Dirichlet coefficients are not automatically physical rest energies.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix

from . import coupling
from .. import slow_geometry
from ..chiral_candidate import action
from ..hierarchy_diagnosis import native_family_dynamics
from ...native_modular_reference import local_preparation
from ...observer import native, quantum_mixing
from ...repair_scheduler import coins
from tools import strict_json


ROOT = Path(__file__).resolve().parents[3]
SEAMS = tuple(quantum_mixing.seams())
CORE = frozenset(local_preparation.CORE)
FIELD = s.QQ.algebraic_field(s.sqrt(5))


def _exact(matrix):
    return s.ImmutableMatrix(DomainMatrix.from_Matrix(s.Matrix(matrix)).convert_to(FIELD).to_Matrix())


def _dm(matrix):
    return DomainMatrix.from_Matrix(s.Matrix(matrix)).convert_to(FIELD)


def _entries(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


def _rational(value):
    if type(value) is int:
        result = F(value)
    elif type(value) is F:
        if (type(value.numerator) is not int or type(value.denominator) is not int
                or value.denominator <= 0 or math.gcd(value.numerator, value.denominator) != 1):
            raise ValueError("canonical original rational weights required")
        result = F(value.numerator, value.denominator)
    else:
        raise ValueError("exact original integer/Fraction weights required")
    if result < 0:
        raise ValueError("nonnegative native seam weights required")
    return result


def weights(value):
    if type(value) is not dict:
        raise ValueError("an explicit native-seam weight dictionary required")
    for edge in value:
        if (type(edge) is not tuple or len(edge) != 2
                or any(type(port) is not int for port in edge) or edge not in SEAMS):
            raise ValueError("original ordered native integer seams required")
    result = {edge: _rational(value.get(edge, 0)) for edge in SEAMS}
    if sum(result.values(), F()) > 1:
        raise ValueError("accepted native weights plus retained timeout must sum to one")
    return result


def crossing(edge):
    return len(set(edge)&CORE) == 1


@lru_cache(maxsize=1)
def _source():
    states = native.states(binary=True)
    index = {state: i for i, state in enumerate(states)}
    sigma = s.ImmutableMatrix(list(map(s.Rational, local_preparation.reference())))
    coordinate = s.ImmutableMatrix(12, 66, lambda i, j: states[j][i])
    rows, checked = [], 0
    source_maps = {}
    for edge in SEAMS:
        swap = s.zeros(66)
        for column, state in enumerate(states):
            output = list(state)
            output[edge[0]], output[edge[1]] = output[edge[1]], output[edge[0]]
            swap[index[tuple(output)], column] = 1
        target = (s.eye(66)+swap)/2
        for law in ("C", "E"):
            actual = s.zeros(66)
            branches = coins.coin_maps(law, edge)
            for coin, record, mapping in branches:
                for before, after in mapping:
                    actual[after, before] += s.Rational(1, 2)
                    checked += 1
            if actual != target:
                raise ArithmeticError("native recorded instrument differs from the independent swap mixture")
            source_maps[law, edge] = s.ImmutableMatrix(actual)
        defect = (sigma.T*(s.eye(66)-target)*sigma)[0, 0]
        norm_square = sum(value**2 for value in sigma-swap*sigma)
        if defect != norm_square/4 or defect != (s.Rational(157, 90000) if crossing(edge) else 0):
            raise ArithmeticError("exact native reference stationarity coefficient changed")
        rows.append({"seam": list(edge), "crosses_core_boundary": crossing(edge),
            "reference_Dirichlet_defect": str(defect), "swap_reference_difference_squared": str(norm_square),
            "admitted_by_same_reference_stationarity": defect == 0})
    admitted = tuple(edge for edge in SEAMS if not crossing(edge))
    excluded = tuple(edge for edge in SEAMS if crossing(edge))
    if len(admitted) != 22 or len(excluded) != 8 or tuple(sorted(CORE)) not in admitted:
        raise ArithmeticError("native stationary-compatible seam cone changed")
    return {"states": states, "sigma": sigma, "coordinate": coordinate,
            "population_maps": source_maps, "admitted": admitted, "excluded": excluded,
            "rows": tuple(rows), "input_coin_addresses_checked": checked}


def source():
    data = _source()
    return {**data, "population_maps": dict(data["population_maps"]),
            "rows": tuple(deepcopy(row) for row in data["rows"])}


def population_channel(value, law="C"):
    coins._law(law)
    value = weights(value)
    data = source()
    timeout = 1-sum(value.values(), F())
    result = s.Rational(timeout)*s.eye(66)
    for edge, weight in value.items():
        result += s.Rational(weight)*data["population_maps"][law, edge]
    return s.ImmutableMatrix(result)


def stationarity(value, law="C"):
    """Complete one-proposal marginal, with the identity timeout retained."""
    value = weights(value)
    data = source()
    channel = population_channel(value, law)
    defect = (s.eye(66)-channel)*data["sigma"]
    dirichlet = (data["sigma"].T*defect)[0, 0]
    cross_weight = sum((weight for edge, weight in value.items() if crossing(edge)), F())
    if dirichlet != s.Rational(157, 90000)*s.Rational(cross_weight):
        raise ArithmeticError("nonnegative seam stationarity sum identity failed")
    stationary = defect == s.zeros(66, 1)
    if stationary != (cross_weight == 0):
        raise ArithmeticError("stationarity cone does not match actual native channel")
    return {"stationary": stationary, "crossing_seam_weight": str(cross_weight),
        "reference_Dirichlet_defect": str(dirichlet),
        "timeout_probability": str(1-sum(value.values(), F())),
        "work_reference_defect_vector": [str(x) for x in defect]}


def port_generator(value):
    value = weights(value)
    result = s.zeros(12)
    for (a, b), weight in value.items():
        vector = s.zeros(12, 1); vector[a] = 1; vector[b] = -1
        result += s.Rational(weight)*vector*vector.T/2
    return s.ImmutableMatrix(result)


@lru_cache(maxsize=1)
def _family_geometry():
    basis = s.ImmutableMatrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a, b in row]
                              for row in slow_geometry.raw_frame()])
    b = _dm(basis)
    gram = _exact(basis.T*basis)[0, 0]
    constraints = b.transpose().vstack(_dm(s.ones(1, 12)))
    complement = constraints.nullspace().transpose()
    c = s.ImmutableMatrix([5 if port in CORE else -1 for port in range(12)])
    x = _exact(basis.T*c/gram)
    if (_exact(basis.T*basis) != gram*s.eye(3) or _exact(x) != s.Matrix([0, 3*s.sqrt(5)/5, 0])
            or complement.shape != (12, 8)):
        raise ArithmeticError("the actual native family projection changed")
    return {"basis": basis, "gram": gram, "complement": complement,
            "conserved_port_vector": c, "conserved_family_coordinate": x}


def family_geometry():
    geometry = _family_geometry()
    return {**geometry, "complement": deepcopy(geometry["complement"])}


def schur_family(value):
    """Full12-port complement elimination, not a closed three-band fit.

    The displayed derivative is the existing Schur theorem's
    U_harm^T U_harm/G. Its positivity is a norm identity; it is not a
    Lorentzian continuation or a physical wave-function renormalization.
    """
    value = weights(value)
    geometry = family_geometry()
    b, c, d = _dm(geometry["basis"]), geometry["complement"], _dm(port_generator(value))
    complement = c.transpose()*d*c
    if complement.det() == FIELD.zero:
        raise ValueError("this explicit Schur inverse requires a nonsingular complementary form")
    u = b-c*complement.inv()*c.transpose()*d*b
    inverse_gram = FIELD.from_sympy(1/geometry["gram"])
    static = (u.transpose()*d*u).scalarmul(inverse_gram)
    kinetic = (u.transpose()*u).scalarmul(inverse_gram)
    outside = u-b
    identity = _dm(s.eye(3))
    zero = lambda matrix: matrix.to_Matrix() == s.zeros(*matrix.shape)
    if (not zero(b.transpose()*u-_dm(geometry["gram"]*s.eye(3)))
            or not zero(d*u-b*static)
            or not zero(kinetic-identity-(outside.transpose()*outside).scalarmul(inverse_gram))):
        raise ArithmeticError("exact all-complement Schur/derivative identity failed")
    sm, zm = s.ImmutableMatrix(static.to_Matrix()), s.ImmutableMatrix(kinetic.to_Matrix())
    return {"S0": sm, "Z": zm, "U": s.ImmutableMatrix(u.to_Matrix()),
            "rank": static.rank(), "canonical": s.ImmutableMatrix((kinetic.inv()*static).to_Matrix())}


def maximal_stationary_weights():
    """One actual five-bit proposal:22 accepted addresses,10 terminal ones."""
    return {edge: F(1, 32) for edge in source()["admitted"]}


def crossing_countercontrol_weights():
    """A finite64-address proposal crosses once and therefore changes sigma."""
    result = maximal_stationary_weights()
    result[source()["excluded"][0]] = F(1, 64)
    return result


def conservation_certificate():
    """Every actual C/E quantum branch commutes with protected core charge."""
    data = source()
    charges = [sum(state[port] for port in CORE) for state in data["states"]]
    checked = 0
    for law in ("C", "E"):
        for edge in data["admitted"]:
            for coin, label, mapping in coins.coin_maps(law, edge):
                for before, after in mapping:
                    if charges[before] != charges[after]:
                        raise ArithmeticError("an admitted recorded branch changes core charge")
                    checked += 1
    geometry = family_geometry()
    c, x = geometry["conserved_port_vector"], geometry["conserved_family_coordinate"]
    if any(port_generator({edge: 1})*c != s.zeros(12, 1) for edge in data["admitted"]):
        raise ArithmeticError("a positive stationary-compatible Dirichlet term lifts the conserved direction")
    return {"all_recorded_input_addresses_checked": checked,
        "core_charge": "q=n0+n1", "centered_port_readout": "5*sum_core(n)-sum_outer(n)=6*q-2",
        "nonzero_family_coordinate": [str(v) for v in x],
        "all_C_E_branch_commutators_with_core_charge_zero_exact": True,
        "history_scope": "every finite admitted-seam history, including seam choices from prior records, independent randomness and retained timeouts",
        "claim_that_all_record_adaptive_work_marginals_are_stationary_or_linear": False,
        "positive_static_schur_family_rank_upper": 2,
        "positive_kinetic_normalization_can_remove_zero": False}


def higgs_rank_control(static):
    """Reuse the same exact M15/Higgs tensors with the maximal rank-two seed."""
    replacement = s.kronecker_product(s.eye(15), action.family_feature().inv()*static)
    kernel = s.ImmutableSparseMatrix(action.mass_kernel()*replacement)
    charge = action.electric_charge()
    entries = kernel.todok()
    rank = len(entries)
    delta = s.kronecker_product(action.SPIN_EPSILON, kernel)
    if (kernel != kernel.T or charge.T*kernel+kernel*charge != s.zeros(45)
            or delta.T != -delta or rank != 28
            or len({i for i, _ in entries}) != rank or len({j for _, j in entries}) != rank):
        raise ArithmeticError("actual rank-two Higgs/current/CAR attachment failed")
    # The established internal tensor identities are rechecked by their
    # owning function; they tensor with this same real symmetric family seed.
    identities = action.report()["identities"]["one_generation_current_Higgs_identities"]
    return {"actual_M15_Higgs_current_identities_rechecked": identities,
        "K45_rank": rank, "Delta90_rank": 2*rank, "K45_kernel_dimension": 45-rank,
        "massive_charged_Dirac_fields_counting_color": 14,
        "massless_charged_Dirac_fields_counting_color": 7,
        "massless_left_Weyl_neutrino_fields": 3,
        "seed_to_Yukawa_identification_is_a_declared_tested_route": True,
        "canonical_positive_kinetic_forms_preserve_these_ranks": True,
        "physical_charged_mass_or_vacuum_prediction": False}


def source_bindings():
    names = set(native_family_dynamics.SOURCE_FILES)|set(action.report()["source_bindings"])
    names.update({"model/common_source/full_repair_response/"+name for name in (
        "stationary_selection.py", "test_stationary_selection.py", "coupling.py", "README.md")})
    names.update(("model/observer/native.py", "model/observer/quantum_mixing.py",
                  "model/repair_scheduler/coins.py", "model/common_source/slow_geometry.py", "requirements.txt"))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def uniform_reference_control():
    """Actual connected all30-seam population law selects the uniform diagonal."""
    data = source()
    reached = {0}
    while True:
        following = set(reached)
        for edge in SEAMS:
            for coin, record, mapping in coins.coin_maps("C", edge):
                following.update(after for before, after in mapping if before in reached)
        if following == reached:
            break
        reached = following
    uniform = s.ones(66, 1)/66
    if len(reached) != 66 or any(matrix*uniform != uniform for matrix in data["population_maps"].values()):
        raise ArithmeticError("the native connected exclusion-population control changed")
    mean = data["coordinate"]*uniform
    geometry = family_geometry()
    feature = _exact(geometry["basis"].T*s.diag(*mean)*geometry["basis"]/geometry["gram"])
    if mean != s.ones(12, 1)/6 or feature != s.eye(3)/6:
        raise ArithmeticError("uniform stationary source feature is not scalar")
    return {"actual_total_two_population_graph_connected_states": len(reached),
        "all_positive_all30_weights_unique_stationary_population": "I66/66",
        "quantum_coherences_can_feed_native_output_diagonals": False,
        "one_body_port_means": "1/6", "compressed_family_feature": "I3/6",
        "same_biased_reference_or_its_clock_coefficient_inputs_preserved": False,
        "other_nonclassical_family_observables_or_different_clock_algebras_excluded": False}


def report():
    data, geometry = source(), family_geometry()
    stationary_weights = maximal_stationary_weights()
    counter_weights = crossing_countercontrol_weights()
    positive = schur_family(stationary_weights)
    counter = schur_family(counter_weights)
    x = geometry["conserved_family_coordinate"]
    if (positive["rank"] != 2 or positive["S0"]*x != s.zeros(3, 1)
            or counter["rank"] != 3):
        raise ArithmeticError("maximal stationary rank or crossing countercontrol changed")
    if port_generator(stationary_weights).rank() != 10 or port_generator(counter_weights).rank() != 11:
        raise ArithmeticError("the admitted and crossing port-kernel dimensions changed")
    c = geometry["conserved_port_vector"]
    if any((c.T*port_generator({edge: 1})*c)[0, 0] != (18 if crossing(edge) else 0) for edge in SEAMS):
        raise ArithmeticError("the exact crossing-activity energy bound changed")
    prefactor = _exact(s.Matrix([[18/(geometry["gram"]*(x.T*x)[0])]]))[0]
    if prefactor != s.Rational(5, 4)-s.sqrt(5)/4:
        raise ArithmeticError("the source-normalized variational prefactor changed")
    poles = [_exact(s.Matrix([[positive["S0"][j, j]/positive["Z"][j, j]]]))[0] for j in (0, 2)]
    if not positive["S0"].is_diagonal() or not positive["Z"].is_diagonal() or any(v.is_positive is not True for v in poles):
        raise ArithmeticError("the maximal reference-compatible two-gap witness changed")
    all_seams = {edge: F(1, 32) for edge in SEAMS}
    return {"schema": "oph.full_repair_response.stationary_selection.v1",
        "reference": {"binary_total_two_dimension": 66, "core_ports": sorted(CORE),
            "per_configuration_core_weights": ["2/225", "1/36", "2/45"],
            "faithful_exact_reference_used": True, "new_empirical_inputs": []},
        "primitive_mapping": {"laws": ["C", "E"], "seams": 30,
            "input_coin_addresses_checked": data["input_coin_addresses_checked"],
            "diagonal_channel_identity": "E_e=(I+Swap_e)/2; fair native coin and all records included",
            "seam_certificates": list(data["rows"])},
        "stationary_weight_cone": {"state_independent_nonnegative_weights_required": True,
            "mathematical_domain": "all nonnegative real seam weights; executable inputs are exact canonical rationals",
            "accepted_weight_sum_at_most_one": True, "remaining_weight_is_retained_identity_timeout": True,
            "work_marginal_reference_stationary_iff_crossing_weight_zero": True,
            "identity": "sigma^T(I-T)sigma=(1/4)*sum_e w_e||sigma-Swap_e sigma||²=(157/90000)*W_cross",
            "admitted_seams": [list(e) for e in data["admitted"]],
            "excluded_crossing_seams": [list(e) for e in data["excluded"]],
            "core_internal_edge_previously_omitted_is_admitted": [0, 1]},
        "maximal_stationary_source": {"selector_bits": 5, "accepted_addresses": 22,
            "terminal_addresses": 10, "weight_per_admitted_seam": "1/32",
            "C_stationarity": stationarity(stationary_weights, "C"),
            "E_stationarity": stationarity(stationary_weights, "E"),
            "port_Dirichlet_kernel": "span(indicator_core,indicator_outer)",
            "S0": _entries(positive["S0"]), "Z": _entries(positive["Z"]),
            "harmonic_extension": _entries(positive["U"]),
            "static_rank": positive["rank"], "canonical_nonzero_coefficients": list(map(str, poles)),
            "canonical_singular_values_have_one_exact_zero": True,
            "physical_clock_or_Lorentzian_mass_identification": False},
        "persistent_zero": conservation_certificate(),
        "same_M15_Higgs_route": higgs_rank_control(positive["S0"]),
        "nonstationary_positive_countercontrol": {"additional_crossing_seam": list(data["excluded"][0]),
            "additional_weight": "1/64", "C_stationarity": stationarity(counter_weights, "C"),
            "E_stationarity": stationarity(counter_weights, "E"),
            "family_rank": counter["rank"], "S0": _entries(counter["S0"]),
            "same_core_biased_reference_preserved": False,
            "full_rank_is_possible_when_stationarity_requirement_is_dropped": True},
        "all30_source_countercontrol": stationarity(all_seams),
        "different_reference_countercontrol": uniform_reference_control(),
        "quantitative_boundary": {"W_cross": "sum weights of eight core-to-outer seams",
            "reference_Dirichlet_defect_exact": "(157/90000)*W_cross",
            "centered_conserved_port_vector_norm_squared": 60,
            "projected_family_coordinate_norm_squared": "9/5", "family_Gram": str(geometry["gram"]),
            "source_candidate_energy_exact": "c^T D c=18*W_cross",
            "smallest_canonically_normalized_static_coefficient_upper": "((5-sqrt(5))/4)*W_cross",
            "kinetic_contract": "the native Schur derivative Z=U_harm^T U_harm/G>=I; a freely rescaled kinetic metric requires its own quantitative bound",
            "proof": "use c as an admissible extension of its nonzero family projection; harmonic minimization and Z>=I give the Rayleigh upper bound",
            "Dirichlet_defect_is_physical_entropy_production": False,
            "coefficient_bound_is_already_a_physical_mass_bound": False},
        "scope": {"classifies_entire_nonnegative_stationary_native_seam_cone_for_this_reference": True,
            "more_multi_time_current_data_changes_branchwise_conservation": False,
            "any_function_of_native_data_must_have_zero_Yukawa": False,
            "passive_Dirichlet_to_Higgs_route_can_give_three_nonzero_charged_families": False,
            "new_reference_active_recovery_or_other_coherent_mass_law_excluded": False,
            "all_OPH_particle_constructions_excluded": False,
            "observed_particle_spectrum_derived": False},
        "nearest_existing_result": "coupling.py proves the earlier protected rank1 Schur identity and positive derivative. This result classifies all reference-compatible seam weights, adds the allowed core edge and attains the sharp rank2 bound.",
        "source_bindings": source_bindings()}


def compare_exact(retained, current, path="root"):
    if type(retained) is not type(current):
        raise ValueError("stationary selection JSON type changed: "+path)
    if type(current) is dict:
        if retained.keys() != current.keys():
            raise ValueError("stationary selection fields changed: "+path)
        for key in current:
            compare_exact(retained[key], current[key], path+"."+key)
    elif type(current) is list:
        if len(retained) != len(current):
            raise ValueError("stationary selection list changed: "+path)
        for i, (a, b) in enumerate(zip(retained, current, strict=True)):
            compare_exact(a, b, path+f"[{i}]")
    elif retained != current:
        raise ValueError("stationary selection exact evidence changed: "+path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--check", type=Path)
    args = parser.parse_args()
    current = json.loads(json.dumps(report(), allow_nan=False))
    if args.check:
        compare_exact(strict_json.load(args.check), current)
        print("Recomputed all stationary-cone, native-history, Schur and Higgs identities", flush=True)
    else:
        with args.output.open("x") as handle:
            json.dump(current, handle, indent=2, allow_nan=False)
            handle.write("\n")
        print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
