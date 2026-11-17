"""One empirically calibrated M15 x family Weyl--Higgs matter action.

The same exact source currents and internal Higgs tensors as the existing
matter candidate are used with a separate, explicitly calibrated parameter
packet. All nine Yukawas and the Higgs expectation value share a published
MSbar M_Z matching convention. The source of those values, CKM fit, scheme,
and rounding remain inputs. This is an implementation sanity check, not a
derivation of masses, an independent fit test or an interacting pole solver.

Physical four-component Y_d=V_CKM diag(y_d,y_s,y_b). For all-left fields
(chi_L,chi_R^c), the directed mass coefficient is Y_d.conjugate(); its
reverse block is its transpose, not its adjoint. This distinction makes
the actual charged source current, mass diagonalization and CP convention
agree in one weak basis.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import sympy as s

from tools import strict_json
from .chiral_candidate import action as source
from .chiral_candidate import kinetic


ROOT = Path(__file__).resolve().parents[2]
PACKET = Path(__file__).with_name("calibrated_matter_parameters.json")
NEUTRINO_PACKET = Path(__file__).with_name("calibrated_neutrino_parameters.json")
CHANNELS = ("up", "down", "charged_lepton")
LEFT_TYPE = {"up": "Q", "down": "Q", "charged_lepton": "L"}
RIGHT_TYPE = {"up": "u^c", "down": "d^c", "charged_lepton": "e^c"}


def _matrix(value):
    return s.ImmutableSparseMatrix(s.SparseMatrix(value))


def _zero(value):
    return all(x == 0 or s.expand(x) == 0 or s.simplify(x) == 0 for x in value)


def _decimal(value, *, positive=True):
    if type(value) is not str:
        raise ValueError("exact decimal parameter strings required")
    try:
        result = s.Rational(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("finite exact decimal parameter required") from exc
    if (result.is_finite is not True or (positive and result <= 0)):
        raise ValueError("positive finite calibration required")
    return result


def parameters():
    """Detached explicit empirical packet; no fitted quantities are hidden."""
    packet = strict_json.load(PACKET)
    if packet.get("schema") != "oph.calibrated_matter.parameters.v1":
        raise ValueError("wrong calibrated matter parameter schema")
    if set(packet["matched_SM"]["yukawa"]) != set(CHANNELS):
        raise ValueError("all three charged Yukawa channels required")
    for row in packet["matched_SM"]["yukawa"].values():
        if type(row) is not list or len(row) != 3:
            raise ValueError("three ordered calibrated Yukawas required")
        values = tuple(map(_decimal, row))
        if not values[0] < values[1] < values[2]:
            raise ValueError("strictly ordered nonzero calibrated Yukawas required")
    for name in ("v_SM_GeV", "g1_GUT", "g2", "g3", "lambda_paper"):
        _decimal(packet["matched_SM"][name])
    for name in ("sin_theta12", "sin_theta23", "sin_theta13"):
        if not 0 < _decimal(packet["CKM"][name]) < 1:
            raise ValueError("calibrated mixing sine must lie in (0,1)")
    _decimal(packet["CKM"]["delta_radians"])
    return deepcopy(packet)


def _channel(channel):
    if type(channel) is not str or channel not in CHANNELS:
        raise ValueError("up, down or charged_lepton channel required")
    return channel


def yukawa_coordinates(channel):
    return tuple(map(_decimal, parameters()["matched_SM"]["yukawa"][_channel(channel)]))


def neutral_higgs_norm():
    """Running neutral component norm v_SM/sqrt2, in GeV."""
    return _decimal(parameters()["matched_SM"]["v_SM_GeV"])/s.sqrt(2)


def mass_coordinates(channel):
    """Calibrated running action masses, not quark or unstable lepton poles."""
    return tuple(neutral_higgs_norm()*y for y in yukawa_coordinates(channel))


def _rotations(a, b, c, delta):
    ca, cb, cc = (s.sqrt(1-x*x) for x in (a, b, c))
    phase = s.exp(s.I*delta)
    return (_matrix([[1, 0, 0], [0, cb, b], [0, -b, cb]]),
            _matrix([[cc, 0, c/phase], [0, 1, 0], [-c*phase, 0, cc]]),
            _matrix([[ca, a, 0], [-a, ca, 0], [0, 0, 1]]))


def mixing_rotations():
    packet = parameters()["CKM"]
    return _rotations(*(_decimal(packet["sin_theta"+ij]) for ij in ("12", "23", "13")),
                      _decimal(packet["delta_radians"]))


def ckm_matrix():
    r23, r13, r12 = mixing_rotations()
    return _matrix(r23*r13*r12)


def family_yukawas():
    """Physical four-component convention; left d mass eigenvectors are V."""
    values = {name: s.diag(*yukawa_coordinates(name)) for name in CHANNELS}
    values["down"] = ckm_matrix()*values["down"]
    return {name: _matrix(value) for name, value in values.items()}


def yukawa_matrices():
    """Alias for the physical family matrices, not their all-left conjugates."""
    return family_yukawas()


def _oriented_components(channel):
    """Keep the actual internal tensor only from left to conjugated-right."""
    labels = source.ordering()["all_left_components"]
    tensors = source._representation()[4][_channel(channel)]["components"]
    return tuple(_matrix(s.SparseMatrix(15, 15, {
        (i, j): value for (i, j), value in tensor.todok().items()
        if labels[i]["internal_type"] == LEFT_TYPE[channel]
        and labels[j]["internal_type"] == RIGHT_TYPE[channel]})) for tensor in tensors)


def higgs_components():
    """The two symmetric K45 coefficients per channel, in the source order."""
    yukawas = family_yukawas()
    result = {}
    for name in CHANNELS:
        rows = []
        for oriented in _oriented_components(name):
            directed = s.kronecker_product(oriented, yukawas[name].conjugate())
            rows.append(_matrix(directed+directed.T))
        result[name] = tuple(rows)
    return result


def _higgs_norm(value):
    if value is None:
        return neutral_higgs_norm()
    if isinstance(value, (bool, float)):
        raise ValueError("an exact real nonnegative Higgs norm required")
    value = s.sympify(value)
    if value.free_symbols or value.is_real is not True or value.is_nonnegative is not True or value.atoms(s.Float):
        raise ValueError("an exact real nonnegative Higgs norm required")
    return value


def mass_kernel(v=None, *, component=1, include_neutrinos=False):
    if type(component) is not int or component not in (0, 1):
        raise ValueError("Higgs component0 or1 required")
    if type(include_neutrinos) is not bool:
        raise ValueError("Boolean optional neutrino extension required")
    coefficients = higgs_components()
    value = _higgs_norm(v)
    result = value*sum((coefficients[name][component] for name in CHANNELS), s.zeros(45))
    if include_neutrinos:
        result += neutrino_higgs_kernel((value, 0) if component == 0 else (0, value))
    return _matrix(result)


def pairing_kernel(v=None, *, include_neutrinos=False):
    """Creation coefficient epsilon_spin tensor K*, antisymmetric on90 modes."""
    return _matrix(s.kronecker_product(source.SPIN_EPSILON, mass_kernel(v, include_neutrinos=include_neutrinos).conjugate()))


def mass_basis_change(*, include_neutrinos=False):
    """Weak chi=W chi_mass; source internal signs rephase only singlets."""
    if type(include_neutrinos) is not bool:
        raise ValueError("Boolean optional neutrino extension required")
    labels = source.ordering()["all_left_components"]
    blocks = [s.eye(3) for _ in range(15)]
    for name in CHANNELS:
        for (left, right), sign in _oriented_components(name)[1].todok().items():
            blocks[right] = sign*s.eye(3)
            if name == "down":
                blocks[left] = ckm_matrix()
    if include_neutrinos:
        for row in labels:
            if row["internal_type"] == "L" and row["electric_Q"] == "0":
                blocks[row["index"]] = pmns_matrix()
    return _matrix(s.diag(*blocks))


def weak_raising_current(*, families=True):
    """Actual source dGamma(E34), satisfying [T3,T+]=T+."""
    fundamental = s.zeros(5)
    fundamental[3, 4] = 1
    current = source._exterior_derivative(fundamental, source._representation()[0])
    return _matrix(s.kronecker_product(current, s.eye(3))) if families else current


def currents(*, mass_basis=False, include_neutrinos=False):
    if type(mass_basis) is not bool:
        raise ValueError("Boolean current basis required")
    generators = source.generators()
    if mass_basis:
        w = mass_basis_change(include_neutrinos=include_neutrinos)
        return tuple(_matrix(w.H*g*w) for g in generators)
    return generators


def bosonic_parameters():
    packet = parameters()["matched_SM"]
    g1, gw, gs, lp, v = (_decimal(packet[name]) for name in ("g1_GUT", "g2", "g3", "lambda_paper", "v_SM_GeV"))
    return {"g1_GUT": g1, "gY": s.sqrt(s.Rational(3, 5))*g1, "gW": gw,
            "g3": gs, "lambda_paper": lp, "lambda_SM": lp/4, "v_SM_GeV": v,
            "neutral_higgs_norm_GeV": v/s.sqrt(2)}


def bosonic_tree_control():
    """Tree comparison using the calibrated running g, lambda and loop VEV.

    This control additionally CHOOSES m2_tree=-lambda_SM*v_loop**2. SMDR's
    matched m2 instead follows its loop minimization condition. In particular
    the three zeros below belong to this supplied tree comparison, not the
    tree Hessian of the full SMDR-matched loop action. Pole matching, widths
    and the loop effective potential are separate calibration calculations.
    """
    b = bosonic_parameters()
    gy, gw, lam, v = (b[name] for name in ("gY", "gW", "lambda_SM", "v_SM_GeV"))
    coordinates = s.symbols("h0:4", real=True)
    radius = sum(x*x for x in coordinates)
    potential = -lam*v*v*radius/2+lam*radius*radius/4
    background = dict(zip(coordinates, (0, 0, v, 0), strict=True))
    hessian = s.hessian(potential, coordinates).subs(background)
    if hessian != s.diag(0, 0, 2*lam*v*v, 0):
        raise ArithmeticError("declared tree-comparison Higgs/Goldstone Hessian changed")
    neutral = s.Matrix([[gw*gw, -gw*gy], [-gw*gy, gy*gy]])*v*v/4
    photon = s.Matrix([gy, gw])
    z = s.Matrix([gw, -gy])
    if neutral*photon != s.zeros(2, 1) or not _zero(neutral*z-(gw*gw+gy*gy)*v*v*z/4):
        raise ArithmeticError("the same-action neutral gauge mass matrix changed")
    return {"potential_convention": "V_tree_comparison=-lambda_SM*v_loop^2*(Hdagger H)+lambda_SM*(Hdagger H)^2",
            "tree_comparison_m2_GeV2_exact": str(-lam*v*v),
            "calibrated_v_is_SMDR_loop_minimum_VEV": True,
            "tree_stationarity_is_an_additional_comparison_prescription": True,
            "tree_comparison_m2_is_SMDR_loop_matched_m2": False,
            "tree_Goldstone_zeros_are_claimed_for_full_SMDR_loop_action_Hessian": False,
            "scalar_Hessian_real_coordinates": _serialize(hessian),
            "real_Higgs_background_coordinates": ["0", "0", str(v), "0"],
            "Higgs_radial_mass_squared": str(2*lam*v*v), "Goldstone_zero_directions_exact": 3,
            "charged_W_mass_squared": str(gw*gw*v*v/4), "charged_W_real_directions": 2,
            "neutral_W3_B_mass_squared_matrix": _serialize(neutral),
            "Z_mass_squared": str((gw*gw+gy*gy)*v*v/4),
            "photon_null_vector_W3_B": list(map(str, photon)), "photon_zero_direction_exact": 1,
            "unbroken_color_gluon_zero_directions_exact": 8,
            "physical_pole_masses_or_widths_are_these_tree_coefficients": False}


def neutrino_parameters():
    packet = strict_json.load(NEUTRINO_PACKET)
    if packet.get("schema") != "oph.calibrated_matter.neutrino_parameters.v1":
        raise ValueError("wrong phenomenological neutrino calibration schema")
    inputs, choices = packet["oscillation_inputs"], packet["unmeasured_architecture_choices"]
    for ij in ("12", "23", "13"):
        if not 0 < _decimal(inputs["sin2_theta"+ij]) < 1:
            raise ValueError("oscillation sine squared must lie in (0,1)")
    if _decimal(choices["lightest_mass_eV"], positive=False) < 0:
        raise ValueError("nonnegative lightest neutrino mass required")
    if not 0 < _decimal(inputs["delta_m21_squared_eV2"]) < _decimal(inputs["delta_m31_squared_eV2"]):
        raise ValueError("ordered normal neutrino mass splittings required")
    return deepcopy(packet)


def neutrino_mass_coordinates():
    """Chosen low-energy phenomenological normal masses, in GeV."""
    packet = neutrino_parameters()
    inputs = packet["oscillation_inputs"]
    first = _decimal(packet["unmeasured_architecture_choices"]["lightest_mass_eV"], positive=False)
    return tuple(value/s.Integer(10)**9 for value in
                 (first, s.sqrt(first**2+_decimal(inputs["delta_m21_squared_eV2"])),
                  s.sqrt(first**2+_decimal(inputs["delta_m31_squared_eV2"]))))


def _pmns_rotations():
    packet = neutrino_parameters()["oscillation_inputs"]
    return _rotations(*(s.sqrt(_decimal(packet["sin2_theta"+ij])) for ij in ("12", "23", "13")),
                      _decimal(packet["delta_degrees"])*s.pi/180)


def pmns_matrix():
    r23, r13, r12 = _pmns_rotations()
    choices = neutrino_parameters()["unmeasured_architecture_choices"]
    phases = s.diag(1, *(s.exp(s.I*_decimal(choices[name], positive=False)/2)
                        for name in ("Majorana_phase21_radians", "Majorana_phase31_radians")))
    return _matrix(r23*r13*r12*phases)


def neutrino_family_kernel():
    """Symmetric all-left Mnu=U_PMNS* diag(mnu) U_PMNS dagger."""
    u = pmns_matrix()
    return _matrix(u.conjugate()*s.diag(*neutrino_mass_coordinates())*u.H)


def _LH_tensor():
    artifact = source._pinned("neutrino_action_operators.json", "neutrino_action_extension.json")
    return source._load_sparse(artifact["LH_invariant"])


def neutrino_higgs_kernel(higgs):
    """Optional same-Higgs Weinberg K45; coefficient has inverse-GeV units.

    This explicitly adds a calibrated dimension-five interaction. It does
    not derive its coefficient, unknown phases or scale transport.
    """
    if type(higgs) is not tuple or len(higgs) != 2:
        raise ValueError("two exact Higgs components required")
    values = tuple(s.sympify(x) for x in higgs)
    if any(x.atoms(s.Float) or x.has(s.nan, s.zoo, s.oo, -s.oo) or x.is_finite is False for x in values):
        raise ValueError("exact Higgs components required")
    ell = _LH_tensor()*s.Matrix(values)
    return _matrix(s.kronecker_product(ell*ell.T, neutrino_family_kernel()/neutral_higgs_norm()**2))


def neutrino_report():
    rotations = _pmns_rotations()
    if any(not _zero(r.H*r-s.eye(3)) for r in rotations):
        raise ArithmeticError("calibrated PMNS rotations are not exactly unitary")
    h0, h1 = s.symbols("H0 H1")
    tensor = _LH_tensor()
    ell = tensor*s.Matrix([h0, h1])
    identities = 0
    for primitive, current in zip(source._representation()[1], source.generators(families=False), strict=True):
        derivative = tensor*primitive[3:, 3:]*s.Matrix([h0, h1])
        if not _zero(current.T*ell+derivative):
            raise ArithmeticError("actual same-Higgs LH invariant failed")
        identities += 1
    masses = neutrino_mass_coordinates()
    inputs = neutrino_parameters()["oscillation_inputs"]
    if s.simplify((masses[1]**2-masses[0]**2)*10**18-_decimal(inputs["delta_m21_squared_eV2"])) != 0:
        raise ArithmeticError("calibrated solar splitting identity failed")
    if s.simplify((masses[2]**2-masses[0]**2)*10**18-_decimal(inputs["delta_m31_squared_eV2"])) != 0:
        raise ArithmeticError("calibrated atmospheric splitting identity failed")
    kernel = mass_kernel(include_neutrinos=True)
    if kernel.T != kernel or not _zero(source.electric_charge().T*kernel+kernel*source.electric_charge()):
        raise ArithmeticError("optional Majorana coefficient violates all-left charge or symmetry")
    return {"calibration": neutrino_parameters(), "same_source_current_LH_invariant_identities_exact": identities,
            "family_lift": "the checked LH invariant squares and tensors with the same symmetric calibrated kappa; no new charged field or Higgs added",
            "PMNS_rotation_unitarity_exact": True, "PMNS": _serialize(pmns_matrix()),
            "neutrino_mass_parameters_GeV_exact": list(map(str, masses)),
            "neutrino_mass_parameters_eV_numeric": [float(x*10**9) for x in masses],
            "same_action_mass_basis": "W_neutrino=U_PMNS; W^T Mnu W=diag(mnu) follows from exact unitarity",
            "weak_current_mass_basis": "charged lepton basis I; neutrino-to-lepton T+ block U_PMNS dagger",
            "both_calibrated_mass_squared_splittings_exact": True,
            "full_K45_with_optional_Weinberg_symmetric_and_electric_neutral_exact": True,
            "neutrino_family_coefficient": _serialize(neutrino_family_kernel()),
            "source_selected_neutrino_masses_or_Majorana_phases": False,
            "low_energy_to_MZ_EFT_running_or_matching_error_certified": False}


def mode_labels():
    return source.mode_labels()


def matrices(*, include_neutrinos=False):
    return {"generators15": source.generators(families=False), "generators45": currents(),
            "Q15": source.electric_charge(families=False), "Q45": source.electric_charge(),
            "Q90": source.spin_electric_charge(), "Yukawa": family_yukawas(),
            "CKM": ckm_matrix(), "K45": mass_kernel(include_neutrinos=include_neutrinos),
            "Delta90": pairing_kernel(include_neutrinos=include_neutrinos),
            "mass_basis_change45": mass_basis_change(include_neutrinos=include_neutrinos), "weak_raising45": weak_raising_current(),
            "epsilon2": source.SPIN_EPSILON}


def weyl_symbol(momentum, *, include_neutrinos=False):
    """Same-action finite free symbol in natural GeV units, c=hbar=1.

    The existing Weyl/CAR implementation receives the calibrated full
    kernel explicitly. This does not identify its coordinates or evolution
    parameter with an independently reconstructed physical spacetime.
    """
    return kinetic.bdg_symbol(momentum, kernel=mass_kernel(include_neutrinos=include_neutrinos), speed=1.)


def source_bindings():
    # Closed declared data/code closure. Existing native family coefficients
    # are not evaluated or used to calibrate this separate action.
    paths = {*source.SOURCE_FILES, "model/common_source/chiral_candidate/action.py",
             "model/common_source/chiral_candidate/kinetic.py",
             "model/common_source/calibrated_matter.py", "model/common_source/test_calibrated_matter.py",
             "model/common_source/calibrated_matter_parameters.json", "tools/strict_json.py",
             "model/common_source/calibrated_neutrino_parameters.json", "model/requirements.txt",
             "model/observer_particles/results/neutrino_action_operators.json",
             "model/observer_particles/results/neutrino_action_extension.json",
             "model/observer_particles/neutrino_action_extension.py"}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def _serialize(matrix):
    return {"shape": list(matrix.shape), "entries": [[i, j, str(value)]
            for (i, j), value in sorted(_matrix(matrix).todok().items())]}


def exact_checks():
    rotations = mixing_rotations()
    if any(not _zero(r.H*r-s.eye(3)) for r in rotations):
        raise ArithmeticError("CKM rotation is not exactly unitary")
    components, identities = higgs_components(), 0
    for name in CHANNELS:
        for primitive, current in zip(source._representation()[1], currents(), strict=True):
            h = primitive[3:, 3:]
            if source._representation()[4][name]["conjugate_Higgs"]:
                h = -h.T
            for a in (0, 1):
                k = components[name][a]
                correction = sum((h[b, a]*components[name][b] for b in (0, 1)), s.zeros(45))
                if not _zero(current.T*k+k*current+correction):
                    raise ArithmeticError("calibrated current/Higgs Ward identity failed")
                identities += 1
    k, d, q, w = mass_kernel(), pairing_kernel(), source.electric_charge(), mass_basis_change()
    if k.T != k or d.T != -d or not _zero(q.T*k+k*q):
        raise ArithmeticError("calibrated all-left charge/CAR convention failed")
    if not _zero(w.H*w-s.eye(45)):
        raise ArithmeticError("the physical mass basis is not unitary")
    diagonal = _matrix(w.T*k*w)
    expected = s.zeros(45)
    for name in CHANNELS:
        for (i, j), _ in _oriented_components(name)[1].todok().items():
            for family, mass in enumerate(mass_coordinates(name)):
                expected[3*i+family, 3*j+family] = mass
                expected[3*j+family, 3*i+family] = mass
    if not _zero(diagonal-expected):
        raise ArithmeticError("same-action mass basis failed calibrated singular-value identity")
    plus = weak_raising_current(families=False)
    t3 = source.identification_operators(families=False)["T3"]
    if not _zero(t3*plus-plus*t3-plus):
        raise ArithmeticError("source weak-current normalization changed")
    physical = w.H*weak_raising_current()*w
    quark_pairs = [(i, j) for (i, j), value in plus.todok().items()
                   if source.ordering()["all_left_components"][i]["internal_type"] == "Q"]
    if len(quark_pairs) != 3 or any(not _zero(physical[3*i:3*i+3, 3*j:3*j+3]-ckm_matrix()) for i, j in quark_pairs):
        raise ArithmeticError("the actual charged source current does not give calibrated CKM")
    return {"source_current_Higgs_identities_exact": identities,
            "CKM_rotation_unitarity_exact": True, "K45_symmetric_exact": True,
            "creation_Delta90_antisymmetric_exact": True, "neutral_electric_pairing_exact": True,
            "same_action_mass_diagonalization_exact": True, "quark_color_current_blocks_equal_CKM_exact": 3,
            "canonical_family_kinetic_metric": "I3", "charged_K45_rank_exact": 42,
            "massless_left_neutrinos_in_minimal_renormalizable_action": 3}


def report():
    packet = parameters()
    ckm = ckm_matrix()
    j = s.im(ckm[0, 0]*ckm[1, 1]*s.conjugate(ckm[0, 1]*ckm[1, 0]))
    return {"schema": "oph.calibrated_matter.v1", "calibration": packet,
            "basis": {"weak_all_left_order": source.ordering(),
                      "physical_family_Yukawa": "Yu=diag(yu,yc,yt),Yd=V_CKM*diag(yd,ys,yb),Ye=diag(ye,ymu,ytau)",
                      "all_left_directed_coefficient": "K_L,Rc=Yphysical.conjugate()*H_component; reverse block transpose",
                      "creation_pairing": "epsilon_spin tensor K.conjugate()", "family_axes_are_calibrated_species_labels": True},
            "kinetic_contract": {"metric": "canonical I45", "Weyl_principal_symbol": "-sigma.p tensor I45",
                                 "natural_units": "c=hbar=1; momentum and mass in GeV, programmed evolution parameter in inverse GeV",
                                 "same_native_clock_or_reconstructed_spacetime_calibration_derived_here": False},
            "identities": exact_checks(), "CKM": _serialize(ckm),
            "Jarlskog_exact": str(s.simplify(j)), "Jarlskog_numeric_diagnostic": float(s.N(j, 17)),
            "running_mass_parameters_GeV": {name: list(map(str, mass_coordinates(name))) for name in CHANNELS},
            "running_mass_parameters_numeric_GeV": {name: [float(x) for x in mass_coordinates(name)] for name in CHANNELS},
            "bosonic_parameters_exact": {name: str(value) for name, value in bosonic_parameters().items()},
            "bosonic_tree_control": bosonic_tree_control(),
            "optional_phenomenological_neutrino_extension": neutrino_report(),
            "mass_kernel": _serialize(mass_kernel()), "mass_basis_change": _serialize(mass_basis_change()),
            "scope": {"observed_parameters_are_calibrations": True, "same_fitted_mass_recovery_is_independent_validation": False,
                      "operator_current_Higgs_consistency_checked_in_one_action": True,
                      "finite_tree_action_replaces_SMDR_matching_or_full_quantum_SM": False,
                      "bosonic_tree_comparison_is_full_SMDR_loop_effective_potential": False,
                      "quark_running_parameters_are_free_asymptotic_poles": False,
                      "unstable_widths_or_interacting_pole_self_energies_computed": False,
                      "nonperturbative_chiral_measure_or_continuum_limit_constructed": False,
                      "source_pulses_executed_by_this_operator_module": False,
                      "physical_parameters_derived_from_native_repair": False},
            "source_bindings": source_bindings()}


def _compare(left, right, path="root"):
    if type(left) is not type(right):
        raise ValueError("calibrated matter evidence type changed: "+path)
    if type(right) is dict:
        if left.keys() != right.keys():
            raise ValueError("calibrated matter fields changed: "+path)
        for key in right:
            _compare(left[key], right[key], path+"."+key)
    elif type(right) is list:
        if len(left) != len(right):
            raise ValueError("calibrated matter list changed: "+path)
        for i, (a, b) in enumerate(zip(left, right, strict=True)):
            _compare(a, b, path+f"[{i}]")
    elif left != right:
        raise ValueError("calibrated matter evidence changed: "+path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--check", type=Path)
    args = parser.parse_args()
    current = json.loads(json.dumps(report(), allow_nan=False))
    if args.check:
        _compare(strict_json.load(args.check), current)
        print("Recomputed calibrated common-scale matter action, masses and source currents")
    else:
        with args.output.open("x") as handle:
            json.dump(current, handle, indent=2, allow_nan=False)
            handle.write("\n")


if __name__ == "__main__":
    main()
