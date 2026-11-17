"""Regional vacuum entropy and independent fixed-volume geometry variation.

The angular matrix is the actual twelve-port native Laplacian, including its
constant mode. The supplied radial extension is the existing three-site
Dirichlet grid on (0,L). Its scalar Hamiltonian is

    H = p.p/(2 mu) + mu q.B.q/2,
    B = K_radial tensor I + angular_scale I tensor L_native.

Canonical site fields absorb radial-cell and angular-quadrature weights;
undoing that local rescaling does not change the regional entropy. This is
the declared unrestricted Gaussian Hamiltonian, not a claim that its vacuum
has been prepared by finite native compiler words or that the finite source
arc frequency equals sqrt(L_native) exactly.

Trace the two exterior radial layers, retaining the central layer. This is a
radial slab, not a geodesic ball. Ground-state Gaussian reduction uses the
standard harmonic-lattice covariance theorem (Peschel, cond-mat/0212631;
Plenio et al., quant-ph/0405142). No thermal occupation ansatz is substituted.

For angular stiffness c and radial stiffness kappa, its one-mode symplectic
eigenvalue obeys nu^2=(1+a/sqrt(a^2-2 kappa^2))/8, a=2 kappa+c. Let
angular_scale(z)=angular_scale(0) exp(2 u z), L(z)=L(0) exp(v z). The
action-derived metric has fixed slab proper volume exactly when
mu(z)=mu(0) exp(-2(u+v)z/3). With r=c/kappa, independently,

  d(nu^2)/dz = -(u+v) r/[2((2+r)^2-2)^(3/2)],
  d log(boundary area)/dz = -2(u+v)/3.

Every positive native angular mode therefore gives the SAME sign for the
regional vacuum entropy and boundary-area changes. No positive coefficient
of the bare area term makes their sum stationary when u+v is nonzero.
Allowing the radial coordinate length to respond does not repair this sign
within the uniform scaling class; u+v=0 only preserves the spatial geometry
and rescales the model-coordinate clock. A responsive state, a different
geometry variation, or a justified ultraviolet/subtraction dictionary is
needed before this candidate can supply generalized-entropy stationarity.
In fact the fixed site counts give proper radial spacing rate 2(u+v)/3
and proper tangential spacing rate -(u+v)/3. The comparison changes the
proper ultraviolet regulator anisotropically; keeping one isotropic proper
cutoff is an additional refinement/identification problem, not a consequence
of fixed native state count or fixed slab volume.

These are exact finite comparison identities and floating covariance checks.
They do not identify this slab with a small ball, its bare lattice entropy
with renormalized gravitational entropy, the supplied vacuum with a source-
selected physical state, or this comparison deformation with native repair.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as s

from tools import strict_json
from ..common_source import slow_geometry
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..physical_space_time.observation_process import product_mode_dictionary


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    "model/native_equilibration/vacuum_geometry.py",
    "model/native_equilibration/test_vacuum_geometry.py",
    "model/common_source/slow_geometry.py", "model/observer/native.py",
    "model/physical_space_time/observation_process.py",
    "model/relative_observer_geometry/intrinsic_geometry.py",
    "Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean",
    "evidence/native_geometric_source_20260925/spec.json",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py",
)


def _rational(value, *, positive=False):
    if type(value) not in (int, F) or (positive and value <= 0):
        raise ValueError("original exact rational parameters with the declared sign required")
    result = F(value)
    try:
        numerical = float(result)
    except OverflowError as error:
        raise ValueError("parameters must remain finite in execution precision") from error
    if (not math.isfinite(numerical) or (result != 0 and numerical == 0)
            or (positive and numerical <= 0)):
        raise ValueError("parameters must remain finite and positive in execution precision")
    return result


@lru_cache(maxsize=1)
def native_angular():
    """Bind the spectrum to all actual native endpoint/coin branches."""
    calibrated = slow_geometry.calibration()
    lap = s.Matrix(calibrated["laplacian"])
    x = s.Symbol("x")
    characteristic = s.factor(lap.charpoly(x).as_expr())
    if characteristic != x*(x-6)**5*(x*x-10*x+20)**3:
        raise ArithmeticError("the actual native angular spectrum changed")
    if lap != lap.T or lap*s.ones(12, 1) != s.zeros(12, 1):
        raise ArithmeticError("native angular form is not symmetric and conservative")
    return tuple(tuple(int(entry) for entry in row) for row in lap.tolist())


def exact_identities():
    """Symbolic covariance, geometric determinant and shape derivative checks."""
    kappa, c, r = s.symbols("kappa c r", positive=True)
    a = 2*kappa+c
    low, high = s.sqrt(a-s.sqrt(2)*kappa), s.sqrt(a+s.sqrt(2)*kappa)
    central_x, central_p = (1/low+1/high)/4, (low+high)/4
    nu_squared = (1+a/s.sqrt(a*a-2*kappa*kappa))/8
    if s.simplify((central_x*central_p-nu_squared).subs(c, r*kappa)) != 0:
        raise ArithmeticError("three-site vacuum covariance identity failed")
    dimensionless = (1+(2+r)/s.sqrt((2+r)**2-2))/8
    derivative = -1/(4*((2+r)**2-2)**s.Rational(3, 2))
    if s.simplify(s.diff(dimensionless, r)-derivative) != 0:
        raise ArithmeticError("vacuum entanglement monotonicity identity failed")
    mu, length, aa, dd = s.symbols("mu length aa dd", positive=True)
    bb, u, v, z = s.symbols("bb u v z", real=True)
    screen = s.Matrix([[aa, bb], [bb, (bb*bb+dd)/aa]])
    mu_z = mu*s.exp(-2*(u+v)*z/3)
    action = mu_z*s.diag(-1, 1, s.exp(2*u*z)*screen)
    metric = s.simplify(s.sqrt(-action.det())*action.inv())
    area = s.simplify(2*s.sqrt(metric[2:4, 2:4].det()))
    volume = s.simplify(length*s.exp(v*z)*s.sqrt(metric[1:4, 1:4].det())/4)
    if (s.simplify(area-2*mu_z) != 0
            or s.simplify(volume-length*mu**s.Rational(3, 2)*dd**s.Rational(1, 4)/4) != 0):
        raise ArithmeticError("independent action-to-area/volume identity failed")
    area_rate = s.simplify(s.diff(area, z).subs(z, 0)/(2*mu))
    if area_rate != -2*(u+v)/3:
        raise ArithmeticError("fixed-volume area response changed")
    radial_spacing = length*s.exp(v*z)*s.sqrt(metric[1, 1])/4
    radial_rate = s.simplify((s.diff(radial_spacing, z)/radial_spacing).subs(z, 0))
    tangential_rate = s.simplify((s.diff(s.sqrt(mu_z), z)/s.sqrt(mu_z)).subs(z, 0))
    if radial_rate != 2*(u+v)/3 or tangential_rate != -(u+v)/3:
        raise ArithmeticError("independent proper regulator-spacing response failed")
    return {"central_symplectic_eigenvalue_squared": str(nu_squared),
            "derivative_per_stiffness_ratio": str(derivative),
            "fixed_volume_mu_log_rate": "-2*(u+v)/3",
            "two_face_area_log_rate": str(area_rate),
            "proper_radial_lattice_spacing_log_rate": str(radial_rate),
            "proper_tangential_lattice_spacing_log_rate": str(tangential_rate),
            "proper_regulator_aspect_ratio_log_rate": str(s.simplify(radial_rate-tangential_rate)),
            "slab_proper_volume_exact": str(volume),
            "radial_coordinate_length_rate_repairs_sign_mismatch": False,
            "isometric_spatial_case": "u+v=0; area and proper radial length fixed",
            "identities_use_independent_metric_determinants": True}


def _mode(eigenvalue, kappa, scale, density, shape_rate):
    c = scale*eigenvalue
    r = c/kappa
    a = 2+r
    discriminant = math.sqrt(a*a-2)
    # Rationalize nu^2-1/4; subtracting nearly equal numbers would erase a
    # real positive occupation and falsely report a pure regional state.
    excess = 1/(4*discriminant*(a+discriminant))
    nu = math.sqrt(.25+excess)
    occupation = excess/(nu+.5)
    gap = math.log1p(1/occupation)
    frequency_minus = math.sqrt(kappa*(a-math.sqrt(2)))
    frequency_plus = math.sqrt(kappa*(a+math.sqrt(2)))
    x = (1/frequency_minus+1/frequency_plus)/(4*density)
    p = density*(frequency_minus+frequency_plus)/4
    entropy = math.log1p(occupation)+occupation*gap
    slope = -shape_rate*r*gap/(4*nu*discriminant**3)
    return {"native_angular_eigenvalue": eigenvalue, "symplectic_eigenvalue": nu,
            "occupation": occupation, "entropy_nats": entropy, "entropy_derivative": slope,
            "position_covariance": x, "momentum_covariance": p,
            "modular_position_coefficient": gap*math.sqrt(p/x),
            "modular_momentum_coefficient": gap*math.sqrt(x/p)}


def covariance(*, length=F(2), angular_scale=F(1), density=F(1)):
    """Actual 36-oscillator ground state and its 12-mode regional reduction."""
    length, scale, density = (_rational(value, positive=True)
                              for value in (length, angular_scale, density))
    kappa = float(16/(length*length))
    lap = np.asarray(native_angular(), float)
    radial = kappa*np.array([[2., -1., 0.], [-1., 2., -1.], [0., -1., 2.]])
    potential = np.kron(radial, np.eye(12))+float(scale)*np.kron(np.eye(3), lap)
    values, vectors = np.linalg.eigh(potential)
    if min(values) <= 0 or not np.isfinite(values).all():
        raise ValueError("the declared positive radial vacuum must be numerically resolvable")
    position = (vectors/(2*float(density)*np.sqrt(values)))@vectors.T
    momentum = (vectors*(float(density)*np.sqrt(values)/2))@vectors.T
    region = slice(12, 24)
    return {"potential": potential, "position": position, "momentum": momentum,
            "regional_position": position[region, region], "regional_momentum": momentum[region, region]}


def regional_modular(*, length=F(2), angular_scale=F(1), density=F(1)):
    """Modular quadratic form, independently compared to local energy below."""
    length, scale, density = (_rational(value, positive=True)
                              for value in (length, angular_scale, density))
    values, vectors = np.linalg.eigh(np.asarray(native_angular(), float))
    # The constant mode is certified exactly by the Laplacian, rather than
    # identified by a tolerance in a covariance eigensolver.
    exact = [0.]+[5-math.sqrt(5)]*3+[6.]*5+[5+math.sqrt(5)]*3
    if np.max(abs(values-exact)) > 1e-12:
        raise ArithmeticError("native spectral basis disagrees with exact multiplicities")
    modes = [_mode(lam, float(16/(length*length)), float(scale), float(density), 0.) for lam in exact]
    q = (vectors*np.array([row["modular_position_coefficient"] for row in modes]))@vectors.T
    p = (vectors*np.array([row["modular_momentum_coefficient"] for row in modes]))@vectors.T
    return q, p


def run(*, length=F(2), angular_scale=F(1), density=F(1),
        angular_log_rate=F(1), radial_length_log_rate=F(0)):
    length, scale, density = (_rational(value, positive=True)
                              for value in (length, angular_scale, density))
    u, v = map(_rational, (angular_log_rate, radial_length_log_rate))
    shape = float(u+v)
    if u+v != 0 and shape == 0:
        raise ValueError("the nonzero geometric shape rate must remain representable")
    kappa = float(16/(length*length))
    # Every radial eigenvalue uses the existing public grid dictionary.
    radial = [product_mode_dictionary(0, j, 3, length)["radial_eigenvalue"] for j in (1, 2, 3)]
    if any(s.simplify(a-b) != 0 for a, b in zip(radial,
            (s.Rational(16)/(s.Rational(length)**2)*(2-s.sqrt(2)),
             s.Rational(32)/(s.Rational(length)**2),
             s.Rational(16)/(s.Rational(length)**2)*(2+s.sqrt(2))), strict=True)):
        raise ArithmeticError("the supplied radial dictionary changed")
    exact_values = ((s.S.Zero, 1), (5-s.sqrt(5), 3), (s.Integer(6), 5), (5+s.sqrt(5), 3))
    modes = []
    energy, force = 0., 0.
    for value, multiplicity in exact_values:
        lam = float(value)
        row = _mode(lam, kappa, float(scale), float(density), shape)
        row["exact_native_eigenvalue"], row["multiplicity"] = str(value), multiplicity
        modes.append(row)
        for stiffness in radial:
            stiffness = float(stiffness)
            omega = math.sqrt(stiffness+float(scale)*lam)
            energy += multiplicity*omega/2
            force += multiplicity*(-float(v)*stiffness+float(u*scale)*lam)/(2*omega)
    cov = covariance(length=length, angular_scale=scale, density=density)
    qmod, pmod = regional_modular(length=length, angular_scale=scale, density=density)
    entropy = math.fsum(row["multiplicity"]*row["entropy_nats"] for row in modes)
    slope = math.fsum(row["multiplicity"]*row["entropy_derivative"] for row in modes)
    offdiagonal = pmod-np.diag(np.diag(pmod))
    # Icosahedral symmetry makes the thirty seam entries equal. Selecting
    # their floating maximum depends on the eigensolver's basis/rounding.
    # Fix the experiment geometrically, before evaluating its response.
    laplacian = native_angular()
    first, second = next((a, b) for a in range(12) for b in range(a+1, 12)
                         if laplacian[a][b] < 0)
    # Independent admissible state directions: add epsilon*u*u^T to the
    # momentum covariance with u=(e_first +/- e_second)/sqrt(2). These
    # classical Gaussian momentum displacements leave all qq/qp moments
    # fixed and have equal on-site pp moments. Canonical port-local stress
    # cannot distinguish them, whereas the actual modular pairing can.
    entropy_rates = [(pmod[first, first]+pmod[second, second]+sign*2*pmod[first, second])/4
                     for sign in (-1, 1)]
    report = {"schema": "oph.native_equilibration.vacuum_geometry.v1",
        "exact_identities": exact_identities(),
        "domain": {"native_ports": 12, "native_seams": 30, "radial_sites": 3,
            "topology": "three-site interval times the existing twelve-port icosahedral screen",
            "radial_boundary_condition": "Dirichlet at0,L", "radial_spacing": str(length/4),
            "coordinate_length": str(length), "angular_scale": str(scale), "action_density": str(density),
            "angular_coordinate_measure": "normalized to one; no physical area unit is inferred",
            "region": "central radial cell times complete native angular screen; two boundary faces",
            "radial_region_coordinate_interval": "(3L/8,5L/8)",
            "canonical_site_fields_absorb_quadrature_weights": True,
            "vacuum": "ground state of p.p/(2mu)+mu*q.B.q/2; complete angular zero mode retained",
            "deformation_parameter": "dimensionless z; not native attempt count or physical time",
            "angular_half_log_rate": str(u), "radial_coordinate_length_log_rate": str(v)},
        "angular_modes": modes, "regional_vacuum_entropy_nats": entropy,
        "regional_vacuum_entropy_derivative": slope,
        "two_face_area": float(2*density), "two_face_area_derivative": -4*float(density*(u+v))/3,
        "fixed_slab_proper_volume_derivative_exact": "0",
        "coordinate_clock_zero_point_energy": energy, "independent_action_shape_response": force,
        "pure_global_covariance_residual": float(np.linalg.norm(cov["position"]@cov["momentum"]-.25*np.eye(36), 2)),
        "regional_modular_momentum_offdiagonal_frobenius": float(np.linalg.norm(offdiagonal, "fro")),
        "any_port_local_energy_kernel_equals_this_modular_quadratic_form": False,
        "same_local_stress_first_law_control": {
            "ports": [first, second],
            "port_selection": "lexicographically first native seam, independent of floating modular entries",
            "state_variations": "Gaussian momentum noise along (e_a-e_b)/sqrt2 and (e_a+e_b)/sqrt2",
            "all_port_local_canonical_quadratic_stress_variations_equal": True,
            "all_position_and_mixed_covariance_variations_zero": True,
            "both_coordinate_energy_derivatives": float(1/(2*density)),
            "regional_entropy_derivatives": list(map(float, entropy_rates)),
            "entropy_derivative_difference": float(pmod[first, second]),
            "scope": "the supplied finite slab vacuum; no continuum small-ball kernel is excluded"},
        "bare_entropy_plus_positive_area_coefficient_stationary": shape == 0,
        "scope": {"fixed_volume_uniform_scale_sign_obstruction": shape != 0,
            "radial_coordinate_response_alone_can_repair_sign": False,
            "geometry_or_stress_defined_by_Einstein_equation": False,
            "source_repair_selects_this_vacuum_or_deformation": False,
            "Gaussian_vacuum_prepared_by_native_compiler_words": False,
            "one_isotropic_proper_UV_cutoff_preserved_by_nontrivial_deformation": False,
            "slab_is_a_small_geodesic_ball": False, "bare_entropy_is_renormalized_generalized_entropy": False,
            "ultraviolet_or_common_continuum_error_certified": False,
            "Einstein_equations_derived_or_all_dynamic_geometries_excluded": False},
        "source_sha256": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}}
    def finite(value):
        if isinstance(value, dict):
            return all(finite(child) for child in value.values())
        if isinstance(value, list):
            return all(finite(child) for child in value)
        return type(value) is not float or math.isfinite(value)
    if not finite(report):
        raise ValueError("vacuum entropy, action and geometry diagnostics must remain finite")
    if ((shape > 0 and not slope < 0) or (shape < 0 and not slope > 0)
            or report["pure_global_covariance_residual"] > 1e-10
            or report["regional_modular_momentum_offdiagonal_frobenius"] <= 1e-8
            or abs(pmod[first, second]) <= 1e-8):
        raise ArithmeticError("actual native vacuum sign, purity or modular locality check failed")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", type=Path)
    mode.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print("Native-spectrum vacuum, regional entropy and independent geometry response replayed.")
    elif args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end="")


if __name__ == "__main__":
    main()
