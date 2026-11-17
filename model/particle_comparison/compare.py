"""Exact dimer attachments and an actual-data replay; no Universe membership test.

The finite Hamiltonian is certified in Q(i,sqrt(17)). Its candidate attachment
is tested against already known data. A calibration is never a prediction.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import total_ordering
import hashlib
from itertools import permutations
from pathlib import Path

import numpy as np
import sympy as sp

from tools import strict_json
from model.empirical_diagnostics import cross_scale
from model.full_fock_field.field import FieldGeometry, FullFock
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TARGETS = "model/particle_comparison/observations.json"
CHIRAL = "model/common_source/chiral_candidate/results.json"


def validate_observations(observations):
    """Admit only the named mass observables with their declared source roles.

    Numeric strings remain explicit data inputs. Changing units, replacing a
    mass limit by a g-factor limit, or dropping its assumption is not a new
    value of the same observable and must not reuse this comparison.
    """
    def fields(value, expected, numeric=()):
        if type(value) is not dict or set(value) != set(expected) | set(numeric):
            raise ValueError("physical observation schema mismatch")
        for key, wanted in expected.items():
            if type(value[key]) is not type(wanted) or value[key] != wanted:
                raise ValueError(f"physical observation contract changed: {key}")
        for key in numeric:
            if type(value[key]) is not str:
                raise ValueError("physical observations require exact numeric strings")
            try:
                number = F(value[key])
            except (ValueError, ZeroDivisionError) as exc:
                raise ValueError("finite positive rational observation required") from exc
            if number <= 0:
                raise ValueError("finite positive rational observation required")

    if (type(observations) is not dict or set(observations) !=
            {"schema", "exposure", "electron_calibration", "positron_electron_difference"}):
        raise ValueError("physical observation schema mismatch")
    fields({key:observations[key] for key in ("schema", "exposure")}, {
        "schema":"oph.particle_comparison.observations.v1",
        "exposure":"retrospective_known_data"})
    fields(observations["electron_calibration"], {
        "observable":"electron_mass_energy_equivalent", "unit":"MeV",
        "source":"CODATA 2022 recommended values",
        "citation":"https://physics.nist.gov/cuu/pdf/wall_2022.pdf",
        "role":"one common energy calibration input; not an absolute mass prediction"},
        ("value", "standard_uncertainty"))
    fields(observations["positron_electron_difference"], {
        "observable":"absolute_mass_difference_divided_by_average_mass",
        "confidence_level_percent":90, "unit":"dimensionless",
        "source":"PDG 2025 electron listing, Fee 1993 positronium spectroscopy",
        "citation":"https://pdgweb.lbl.gov/2025/listings/rpp2025-list-electron.pdf",
        "measurement_citation":"https://doi.org/10.1103/PhysRevA.48.192",
        "assumption":"positronium Rydberg constant exactly half the hydrogen one",
        "role":"conditional retrospective equality check; not a second absolute mass measurement",
        "is_g_factor_equality_measurement":False}, ("upper_bound",))
    electron = observations["electron_calibration"]
    if (F(electron["standard_uncertainty"]) >= F(electron["value"])
            or F(observations["positron_electron_difference"]["upper_bound"]) >= 2):
        raise ValueError("nontrivial positive mass calibration and normalized difference bound required")
    return observations


def load_observations(path):
    """Duplicate keys and JSON nonfinite numbers fail before schema admission."""
    return validate_observations(strict_json.loads(Path(path).read_text(encoding="utf-8")))


def rational(value):
    if type(value) not in (int, F):
        raise ValueError("exact integer or Fraction required")
    return F(value)


@total_ordering
@dataclass(frozen=True)
class Quadratic:
    """An exactly ordered element a+b sqrt(17), with rational coefficients."""
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "a", rational(self.a))
        object.__setattr__(self, "b", rational(self.b))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Quadratic) else Quadratic(value)

    def __add__(self, other):
        other = self.coerce(other)
        return Quadratic(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Quadratic(-self.a, -self.b)

    def __sub__(self, other):
        return self+-self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other)+-self

    def __mul__(self, other):
        other = self.coerce(other)
        return Quadratic(self.a*other.a+17*self.b*other.b,
                         self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        denominator = other.a**2-17*other.b**2
        if not denominator:
            raise ZeroDivisionError("zero quadratic element")
        return self*Quadratic(other.a/denominator, -other.b/denominator)

    def sign(self):
        if not self.b:
            return (self.a > 0)-(self.a < 0)
        if not self.a or (self.a > 0) == (self.b > 0):
            return 1 if self.b > 0 else -1
        comparison = self.a**2-17*self.b**2
        return (1 if self.a > 0 else -1)*(1 if comparison > 0 else -1)

    def __lt__(self, other):
        return (self-self.coerce(other)).sign() < 0

    def __abs__(self):
        return self if self.sign() >= 0 else -self

    def sympy(self):
        return sp.Rational(self.a.numerator, self.a.denominator)+sp.sqrt(17)*sp.Rational(
            self.b.numerator, self.b.denominator)

    def as_json(self):
        lower, upper = cross_scale.sqrt_bounds(F(17), 40)
        bounds = sorted((self.a+self.b*lower, self.a+self.b*upper))
        return {"rational": str(self.a), "sqrt17_coefficient": str(self.b),
                "decimal_enclosure": cross_scale._interval(bounds)}


Q = Quadratic


def spectrum():
    """Every eigenvalue, including multiplicities, grouped by number N."""
    return (
        (Q(0),),
        (Q(F(-1,16)),)*2+(Q(F(1,16)),)*2,
        (Q(0),)*3+(Q(1), Q(F(1,2), F(-1,8)), Q(F(1,2), F(1,8))),
        (Q(F(15,16)),)*2+(Q(F(17,16)),)*2,
        (Q(2),),
    )


def certify_source():
    """No fitted roots: certify source matrix characteristic polynomials exactly."""
    geometry = FieldGeometry((2,1,1), speed=.125, interaction=1.)
    space = FullFock(geometry.modes)
    matrix = 16*space.hamiltonian(geometry).toarray()
    if (not np.isfinite(matrix).all() or not np.equal(matrix.real, np.rint(matrix.real)).all()
            or not np.equal(matrix.imag, np.rint(matrix.imag)).all()):
        raise ValueError("source Hamiltonian no longer exactly Gaussian-integer/16")
    exact = sp.Matrix([[sp.Rational(int(z.real),16)+sp.I*sp.Rational(int(z.imag),16)
                        for z in row] for row in matrix])
    if exact != exact.conjugate().T:
        raise ArithmeticError("source Hamiltonian is not Hermitian")
    for i in range(16):
        for j in range(16):
            if i.bit_count() != j.bit_count() and exact[i,j] != 0:
                raise ArithmeticError("source Hamiltonian does not conserve number")
    x = sp.Symbol("x")
    sectors = []
    for n, energies in enumerate(spectrum()):
        indices = [i for i in range(16) if i.bit_count() == n]
        actual = exact.extract(indices,indices).charpoly(x).as_expr()
        expected = sp.expand(sp.prod(x-energy.sympy() for energy in energies))
        if sp.expand(actual-expected) != 0:
            raise ArithmeticError(f"source spectrum differs in sector N={n}")
        sectors.append({"number":n, "dimension":len(indices),
                        "characteristic_polynomial": str(actual),
                        "eigenvalues_with_multiplicity": [e.as_json() for e in sorted(energies)]})
    frequencies = sorted({b-a for low,high in zip(spectrum()[:-1],spectrum()[1:],strict=True)
                          for a in low for b in high})
    return {"geometry": {"sides":[2,1,1], "coordinate_speed":"1/8", "interaction":"1",
                          "unit":"dimensionless Hamiltonian model-coordinate energy"},
            "source_dimension":16, "number_sectors":sectors,
            "method":"exact characteristic polynomials of current source matrix in Gaussian rationals",
            "signed_charge_one_frequency_superset": [f.as_json() for f in frequencies],
            "includes_dark_transitions":True, "visibility_threshold":None,
            "numerical_eigenvalue_tolerance":None}, frequencies


def half_filled_candidates():
    """Ground and charged gaps of the supplied K=H-N/2, before mass data."""
    shifted = [tuple(e-Q(F(n,2)) for e in energies) for n,energies in enumerate(spectrum())]
    ground = min(e for energies in shifted for e in energies)
    locations = [n for n,energies in enumerate(shifted) for e in energies if e == ground]
    if locations != [2]:
        raise ArithmeticError("half-filled ground is not unique in N=2")
    hole, particle = min(shifted[1])-ground, min(shifted[3])-ground
    if hole != particle or hole.sign() <= 0:
        raise ArithmeticError("positive particle-hole gap equality failed")
    return {"hamiltonian":"K=H-(1/2)N", "chemical_potential":"1/2",
            "reference":"supplied zero-temperature unique N=2 ground of K",
            "reference_is_executed_beta1_mu0_Gibbs_state":False,
            "ground_energy":ground.as_json(), "ground_multiplicity":1,
            "gap":hole.as_json(), "exact_normalized_mass_difference":"0",
            "gap_unit":"inverse Hamiltonian model coordinate; hbar_model=1",
            "time_parameter":"Hamiltonian model coordinate, not calibrated observer proper time",
            "candidates":[
                {"assignment":"electron-like hole", "number_sector":1,"relative_number_charge":-1,
                 "multiplicity":shifted[1].count(min(shifted[1]))},
                {"assignment":"positron-like particle", "number_sector":3,"relative_number_charge":1,
                 "multiplicity":shifted[3].count(min(shifted[3]))}],
            "physical_electric_charge_identification_supplied":True,
            "Lorentz_spin_or_physical_particle_identity_derived":False,
            "physical_rest_mass_interpretation_derived":False}, hole


def affine_cases(frequencies, windows):
    """Exhaust m_i=kappa|omega_i-mu|, all assignments, kappa>0, mu real.

    For a,b,c assigned to electron, muon,tau let x=m_mu/m_e>1,
    y=m_tau/m_e>1. The signs of b-mu and c-mu are respectively
    sign(b-a),sign(c-a). The remaining electron sign se gives
    y=abs((c-a)/(b-a))*x+sign(c-a)*se*(1-(c-a)/(b-a)).
    This has positive slope, so rational endpoint tests are exact.
    """
    if (not isinstance(frequencies, (list,tuple)) or len(frequencies) < 3
            or any(not isinstance(f,Q) for f in frequencies)
            or len(set(frequencies)) != len(frequencies)):
        raise ValueError("at least three distinct exact frequencies required")
    if (len(windows) != 2 or any(len(w) != 2 for w in windows)
            or any(type(v) not in (int,F) for w in windows for v in w)):
        raise ValueError("two exact ratio intervals required")
    (xl,xh),(yl,yh) = tuple(tuple(map(F,w)) for w in windows)
    if not 1 < xl <= xh < yl <= yh:
        raise ValueError("ordered disjoint ratios strictly above one required")
    rows = []
    for indices in permutations(range(len(frequencies)),3):
        a,b,c = (frequencies[i] for i in indices)
        q = (c-a)/(b-a)
        st = (c-a).sign()
        for se in (-1,1):
            intercept = st*se*(1-q)
            low, high = abs(q)*xl+intercept, abs(q)*xh+intercept
            feasible = not (high < Q(yl) or Q(yh) < low)
            rows.append({"indices":list(indices), "electron_sign":se,
                         "low":low, "high":high, "feasible":feasible})
    return rows


def lepton_no_go(frequencies, windows):
    rows = affine_cases(frequencies,windows)
    feasible = [r for r in rows if r["feasible"]]
    below = [r for r in rows if r["high"] < Q(windows[1][0])]
    above = [r for r in rows if Q(windows[1][1]) < r["low"]]
    serial = [{**r,"low":r["low"].as_json(),"high":r["high"].as_json()} for r in rows]
    absolute = sorted({abs(w) for w in frequencies if w.sign()})
    return {"attachment":"m_i=kappa*abs(omega_i-mu), common kappa>0 and arbitrary common real mu",
            "target_ratios":[[str(x) for x in w] for w in windows],
            "target_windows_are_joint_confidence_region":False,
            "distinct_ordered_assignments":len(rows)//2, "electron_sign_choices":2,
            "total_cases":len(rows),"feasible_cases":len(feasible),
            "largest_rejected_upper_tau_ratio": max((r["high"] for r in below)).as_json() if below else None,
            "smallest_rejected_lower_tau_ratio": min((r["low"] for r in above)).as_json() if above else None,
            "cases_below_tau_window":len(below), "cases_above_tau_window":len(above),
            "full_case_sha256":hashlib.sha256(cross_scale.canonical(serial).encode()).hexdigest(),
            "zero_mu_absolute_frequency_dynamic_range":(absolute[-1]/absolute[0]).as_json(),
            "repeated_frequency_assignments":"excluded because equal frequencies give equal masses, contrary to disjoint target ratios",
            "status":"bounded_three_lepton_attachment_rejected" if not feasible else "bounded_attachment_has_feasible_cases",
            "scope":"fixed two-site action and charge-one transitions, including dark transitions; no sector-dependent shifts or nonlinear dictionary",
            "universe_or_OPH_program_rejected":False}


def calibrated_pair(prediction, gap, observations):
    observations = validate_observations(observations)
    if not isinstance(gap,Q) or gap.sign() <= 0 or prediction["gap"] != gap.as_json():
        raise ValueError("positive certified common particle-hole gap required")
    predicted_difference = F(prediction["exact_normalized_mass_difference"])
    if predicted_difference != 0:
        raise ValueError("this paired calibration requires the exact equal-gap certificate")
    electron = observations["electron_calibration"]
    difference = observations["positron_electron_difference"]
    centre, uncertainty = F(electron["value"]),F(electron["standard_uncertainty"])
    bound = F(difference["upper_bound"])
    if centre <= 0 or uncertainty <= 0 or bound <= 0:
        raise ValueError("positive physical calibration and uncertainty/bound required")
    scale = Q(centre)/gap
    return {"prediction":prediction, "calibration":electron,
            "calibrated_energy_scale_MeV_per_model_energy":scale.as_json(),
            "candidate_mass_energy_MeV":[electron["value"]]*2,
            "candidate_mass_uncertainty_MeV_from_shared_calibration":[electron["standard_uncertainty"]]*2,
            "calibration_errors_fully_correlated":True,
            "independently_predicted_absolute_masses":0,
            "normalized_difference_prediction":str(predicted_difference), "comparison":difference,
            "within_reported_conditional_bound":abs(predicted_difference) < bound,
            "status":"retrospective_conditional_equality_compatible_after_one_mass_calibration",
            "unique_OPH_signature":False,
            "missing_identity_contracts":["physical vacuum and chemical-potential selection",
                "physical electric charge and Lorentz spin representation",
                "same observer proper clock and physical energy calibration independent of target mass",
                "spatial continuum dispersion and isolated stable pole",
                "common topological defect, transport, interactions and mass on the same excitation"]}


def build(*, root=ROOT):
    root = Path(root)
    certificate, frequencies = certify_source()
    pair_prediction, gap = half_filled_candidates()
    prediction = {"source":certificate,"pair":pair_prediction}
    seal = hashlib.sha256(cross_scale.canonical(prediction).encode()).hexdigest()
    observations = load_observations(root/TARGETS)
    chiral = strict_json.loads((root/CHIRAL).read_text(encoding="utf-8"))
    windows = [tuple(map(F,w)) for w in chiral["separate_empirical_comparison"]["conservative_ratio_windows"]]
    if windows != [(F(206),F(207)),(F(3470),F(3485))]:
        raise ValueError("inherited broad mass comparison windows changed")
    no_go = lepton_no_go(frequencies, windows)
    pair = calibrated_pair(pair_prediction,gap,observations)
    from .dictionary import build as build_dictionary
    dictionary = build_dictionary(root=root)
    from .electromagnetism import build as build_electromagnetism
    electric_response = build_electromagnetism()
    from model.particle_response.dispersion import build as build_dispersion
    dispersion = build_dispersion()
    # This recomputes both the exact spectral prediction and actual ACT weights;
    # it does not merely read a previous cross-scale success receipt.
    replay = cross_scale.build(root=root)
    if hashlib.sha256(cross_scale.canonical(prediction).encode()).hexdigest() != seal:
        raise ArithmeticError("prediction changed after comparison")
    paths = ("model/particle_comparison/compare.py", TARGETS,
             "model/particle_comparison/test_compare.py", "model/particle_comparison/README.md",
             "model/full_fock_field/field.py",CHIRAL,
             "model/particle_comparison/dictionary.py", "model/particle_comparison/test_dictionary.py",
             "model/particle_comparison/electromagnetism.py", "model/particle_comparison/test_electromagnetism.py",
             "model/particle_response/dimer.py", "model/particle_response/spectrum.py",
             "model/gauged_field/links.py",
             "model/particle_response/dispersion.py", "model/particle_response/test_dispersion.py")
    source_bindings = bound_pins({p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    source_bindings.update(bound_pins(replay["custody"]["input_sha256"]))
    return {"schema":"oph.particle_comparison.v1", "classification":"retrospective conditional benchmark",
            "source_bindings":bound_pins(source_bindings),
            "source_spectrum":certificate,"prediction_sha256_before_comparison":seal,
            "sealing_after_prior_target_exposure_makes_prospective":False,
            "two_candidate_assignment":pair, "assigned_operator_dictionary":dictionary,
            "electromagnetic_applicability":electric_response,
            "periodic_dispersion_applicability":dispersion,
            "particle_identity_verdict":{
                "electron_positron_identification_established":False,
                "fixed_dictionary_dispersion_benchmark_satisfied":dispersion["benchmark"]["declared_dimer_anchor_and_speed_benchmark_satisfied"],
                "QED_identification":electric_response["QED_identification"],
                "finite_volume_shift_is_separate_from_same_volume_dispersion_test":True,
                "adjusting_effective_velocity_would_be_an_additional_unverified_attachment":True,
                "whole_OPH_or_continuum_excluded":False},
            "three_lepton_exhaustive_comparison":no_go,
            "actual_data_replay":{
                "full_replay_sha256":hashlib.sha256(cross_scale.canonical(replay).encode()).hexdigest(),
                "coupling":replay["coupling"],
                "ACT_posterior_lower_tail_mass":replay["cosmology"]["forward_posterior_mass_interval"],
                "ACT_posterior_mass_is_p_value":False,
                "cosmology_contract":replay["cosmology"]["physical_comparison_contract"],
                "joint_status":replay["joint_status"], "custody":replay["custody"]},
            "CP11":{"eligible_physical_dictionary":False,"eligible_physical_comparisons_added":0,
                "prospective_comparisons_added":0,"OPH_universe_membership_confirmed":False,
                "global_likelihood_or_probability":None,
                "status":"open: CP07/CP10 common physical dictionary and complete error contract absent"},
            "custody":{"input_sha256":bound_pins({p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}),
                        "exact_current_Hamiltonian_certified":True,
                        "all_assignment_cases_recomputed":True,
                        "existing_actual_data_projection_recomputed":True,
                        "physical_experiment_performed":False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write",type=Path)
    action.add_argument("--check",type=Path)
    args = parser.parse_args()
    result = build()
    if args.check:
        if strict_json.loads(args.check.read_text(encoding="utf-8")) != result:
            raise SystemExit("particle comparison receipt drift")
    else:
        with args.write.open("x",encoding="utf-8") as stream:
            stream.write(cross_scale.canonical(result))
    print(cross_scale.canonical({"pair":result["two_candidate_assignment"]["status"],
        "three_leptons":result["three_lepton_exhaustive_comparison"]["status"],
        "CP11":result["CP11"]}),end="")


if __name__ == "__main__":
    main()
