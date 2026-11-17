"""Gaussian vacuum modular charges for a proper ball with an isotropic cutoff.

This is a locally flat comparison of the existing scalar wave action. A
declared native chart coefficient is frozen, and its constant action metric
defines proper Cartesian coordinates and proper time. It is not asserted to
be the unknown limiting native tensor. In that frame the supplied regulator
is a cubic lattice of spacing a in a Dirichlet box of proper side L. The
region contains precisely the sites with |y| < R. No radial slab or angular
occupation ansatz replaces the three-dimensional vacuum partial trace.

For H=(p.p+q.V.q)/2, V=-Delta_a, the full vacuum has X=V^-1/2/2 and
P=V^1/2/2. Restrict these covariances to the ball, not V itself. If
D=X^1/2 P X^1/2=U diag(nu^2) U^T, the actual modular momentum form is

    Hp=X^1/2 U diag(log((nu+1/2)/(nu-1/2))/nu) U^T X^1/2.

A Gaussian distribution of unobserved momentum displacements changes
P_A to P_A+epsilon*v_A*v_A^T, and has entropy derivative v_A^T Hp v_A/2.
The same number is the relative entropy per squared amplitude of the
coherent displacement: the reduced coherent state is unitarily conjugate
to the reference, so its entropy is unchanged and D(rho_coh||sigma)=Delta K.
Data processing therefore makes this charge monotone under nested regions
of ONE fixed lattice vacuum. This identity is exact; reported charges are
floating executions, with no clipping of unresolved near-pure modes.

The continuum comparison is independently fixed by the conformal scalar
ball generator: 2*pi*integral_B (R^2-r^2)/(2R) T00. See Casini--Huerta,
arXiv:1007.1813 eqs.(4),(8), and Casini--Huerta--Myers, arXiv:1102.0440
eq.(2.20). Momentum-only state directions leave phi^2 unchanged at the
comparison time; hence the scalar improvement term does not enter these
two tests. This is not a claim about every stress component or massive
fields. A separate exact space-form calculation derives the fixed-volume
area coefficient from an independently specified curved metric.

Grid and box controls quantify finite discrepancies, not a proved continuum
convergence rate. A smooth native tangent/vacuum identification, one common
small-ball/modular-error limit, renormalized area density (or measured G),
and selected joint state--geometry stationarity remain separate inputs.
Neither a lattice entropy area coefficient nor an Einstein equation is fit.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.fft import ifftn
import sympy as s

from tools import strict_json
from ..observer.scalar_response.experiment import canonical
from ..physical_space_time import uniform_green
from ..relative_observer_geometry import intrinsic_geometry as geometry
from . import vacuum_geometry


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*vacuum_geometry.SOURCE_PATHS, *geometry.SOURCE_PATHS,
    *uniform_green.SOURCE_PATHS,
    "model/native_equilibration/ball_modular.py",
    "model/native_equilibration/test_ball_modular.py"))))
MIN_PURITY_GAP = 5e-14
ENTROPY_SECANT_REPLAY_TOLERANCE = 3e-6  # divides a small floating entropy difference


def _positive(value):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError("original positive exact rational lengths required")
    try:
        numerical = float(value)
    except OverflowError as error:
        raise ValueError("finite representable proper lengths required") from error
    if not math.isfinite(numerical) or not 1e-40 <= numerical <= 1e40:
        raise ValueError("proper lengths must remain in the declared1e-40..1e40 precision domain")
    return F(value)


def proper_frame():
    """An actual native point read, explicitly frozen as a comparison tensor."""
    read = geometry.native_chart_increment(level=1, attempts=2)
    screen = s.Matrix([[s.Rational(v) for v in row]
                       for row in read["averaged_screen_tensor"]])
    metric = geometry.wave_density_metric(s.diag(1, screen).tolist(), F(1, 10))
    g = metric["metric"]
    # y=E*x and tau=N*t: E.T*E=g_spatial. This identity, not a change
    # of the cutoff between directions, fixes the proper lattice frame.
    frame = g[1:4, 1:4].cholesky().T
    transform = s.diag(s.sqrt(-g[0, 0]), frame)
    if s.simplify(transform.T*s.diag(-1, 1, 1, 1)*transform-g) != s.zeros(4):
        raise ArithmeticError("the proper frame did not recover the same action metric")
    if s.simplify(transform*metric["action_density"]*transform.T/transform.det()
                  -s.diag(-1, 1, 1, 1)) != s.zeros(4):
        raise ArithmeticError("proper coordinates did not canonically normalize the complete action")
    encode = lambda matrix: [[str(s.simplify(v)) for v in row] for row in matrix.tolist()]
    return {"native_chart_read": read, "chosen_action_density": "1/10",
        "constant_comparison_metric": encode(g), "proper_coordinate_matrix": encode(transform),
        "same_action_metric_identity_exact": True,
        "proper_action_density_exact": "diag(-1,1,1,1)",
        "finite_average_identified_with_native_continuum_tensor": False,
        "proper_time_parameter": "tau=sqrt(-g_tt)*model_coordinate_t"}


def exact_ball_identities():
    """Independent ball/stress integrals and constant-curvature geometry."""
    r, R, k = s.symbols("r R k", positive=True)
    weight = (R*R-r*r)/(2*R)
    constant = s.integrate(4*s.pi*r*r*s.pi*weight, (r, 0, R))
    linear = s.integrate(4*s.pi*r**4*s.pi*weight/(3*R*R), (r, 0, R))
    if constant != 4*s.pi**2*R**4/15 or linear != 4*s.pi**2*R**4/105:
        raise ArithmeticError("independent conformal stress integrals changed")
    # ds_space^2=d rho^2+S_k(rho)^2 dOmega^2, S_k=sin(sqrt(k)rho)/sqrt(k).
    # Its exact small-k derivative follows from this metric's area element.
    sphere_radius = s.sin(s.sqrt(k)*r)/s.sqrt(k)
    area_series = s.series(4*s.pi*sphere_radius**2, k, 0, 2).removeO()
    volume_series = s.integrate(area_series, (r, 0, R))
    radius_rate = s.simplify(-s.diff(volume_series, k).subs(k, 0)/(4*s.pi*R*R))
    area_rate = s.simplify(s.diff(area_series, k).subs({k: 0, r: R})
                           +8*s.pi*R*radius_rate)
    if radius_rate != R**3/15 or area_rate != -4*s.pi*R**4/5:
        raise ArithmeticError("independent fixed-volume geodesic-ball expansion changed")
    # For the declared ultrastatic space-form metric the independently
    # computed spatial scalar curvature is 6*k and G_tau,tau=3*k.
    theta, phi, tau = s.symbols("theta phi tau", real=True)
    coordinates = (tau, r, theta, phi)
    g = s.diag(-1, 1/(1-k*r*r), r*r, r*r*s.sin(theta)**2)
    inverse = g.inv()
    @lru_cache(None)
    def gamma(i, j, h):
        return s.simplify(sum(inverse[i, n]*(s.diff(g[n, j], coordinates[h])
            +s.diff(g[n, h], coordinates[j])-s.diff(g[j, h], coordinates[n]))/2
            for n in range(4)))
    ricci = s.Matrix(4, 4, lambda i, j: s.simplify(sum(
        s.diff(gamma(n, i, j), coordinates[n])-s.diff(gamma(n, i, n), coordinates[j])
        +sum(gamma(n, n, h)*gamma(h, i, j)-gamma(n, j, h)*gamma(h, i, n)
             for h in range(4)) for n in range(4))))
    scalar = s.simplify(s.trace(inverse*ricci))
    einstein00 = s.simplify(ricci[0, 0]-g[0, 0]*scalar/2)
    if scalar != 6*k or einstein00 != 3*k:
        raise ArithmeticError("independent space-form curvature computation failed")
    return {"constant_momentum_noise_entropy_kernel": str(constant),
        "linear_momentum_noise_entropy_kernel": str(linear),
        "ball_area": "4*pi*R**2", "ball_volume": "4*pi*R**3/3",
        "space_form_metric": "-d tau**2+d r**2/(1-k*r**2)+r**2*dOmega**2",
        "space_form_scalar_curvature": str(scalar), "space_form_G_tau_tau": str(einstein00),
        "fixed_volume_geodesic_radius_derivative_per_k": str(radius_rate),
        "fixed_volume_area_derivative_per_k": str(area_rate),
        "fixed_volume_area_derivative_per_G_tau_tau": str(s.simplify(area_rate/3)),
        "geometry_defined_from_stress_or_entropy": False,
        "Newton_constant_or_renormalized_area_coefficient_derived": False}


def ball_sites(spacing=F(2, 7), box_length=F(16), radius=F(1)):
    spacing, length, radius = map(_positive, (spacing, box_length, radius))
    count = length/spacing
    if count.denominator != 1 or count % 2 or not 6 <= count <= 80 or length < 4*radius:
        raise ValueError("an even 6..80-cell Dirichlet box of side at least4R required")
    if not F(1) < radius/spacing <= F(9, 2):
        raise ValueError("the declared numerical ball uses1<R/a<=9/2")
    M = int(count)
    reach = math.ceil(radius/spacing)
    offsets = np.array(list(product(range(-reach, reach+1), repeat=3)), dtype=int)
    # Strict inclusion is exact on the original rational grid.
    numerator, denominator = (radius/spacing).numerator, (radius/spacing).denominator
    mask = denominator**2*np.sum(offsets**2, axis=1) < numerator**2
    sites = offsets[mask]+M//2
    if not len(sites) or np.any(sites <= 0) or np.any(sites >= M):
        raise ValueError("the regional algebra must lie inside the Dirichlet box")
    return sites, M, spacing, length, radius


def covariance(spacing=F(2, 7), box_length=F(16), radius=F(1)):
    """Full three-dimensional Dirichlet vacuum restricted to ball sites.

    Eight odd images of a2M-periodic FFT equal the positive Dirichlet sine
    spectrum. The periodic zero mode cancels exactly in the odd extension;
    the physical Dirichlet Hamiltonian has no zero-frequency mode.
    """
    sites, M, spacing, length, radius = ball_sites(spacing, box_length, radius)
    a = float(spacing)
    frequencies = 4/a**2*np.sin(np.pi*np.arange(2*M)/(2*M))**2
    stiffness = frequencies[:, None, None]+frequencies[None, :, None]+frequencies[None, None, :]
    inverse = np.zeros_like(stiffness)
    np.divide(.5, np.sqrt(stiffness), out=inverse, where=stiffness > 0)
    xkernel, pkernel = ifftn(inverse).real, ifftn(np.sqrt(stiffness)/2).real
    X, P = np.zeros((len(sites), len(sites))), np.zeros((len(sites), len(sites)))
    for signs in product((-1, 1), repeat=3):
        indices = (sites[:, None, :]-np.array(signs)*sites[None, :, :]) % (2*M)
        indices = tuple(indices.transpose(2, 0, 1))
        X += math.prod(signs)*xkernel[indices]
        P += math.prod(signs)*pkernel[indices]
    return (X+X.T)/2, (P+P.T)/2, sites, M


def direct_covariance(spacing=F(2, 7), box_length=F(8), radius=F(1)):
    """Independent Dirichlet spectral sum; no FFT or image extension."""
    sites, M, spacing, _, _ = ball_sites(spacing, box_length, radius)
    modes = np.array(list(product(range(1, M), repeat=3)), dtype=int)
    X, P = np.zeros((len(sites), len(sites))), np.zeros((len(sites), len(sites)))
    for start in range(0, len(modes), 2048):
        block = modes[start:start+2048]
        U = (2/M)**1.5*np.prod(np.sin(np.pi*sites[:, None, :]*block[None, :, :]/M), axis=2)
        omega = np.sqrt(np.sum(4/float(spacing)**2*np.sin(np.pi*block/(2*M))**2, axis=1))
        X += (U/(2*omega))@U.T
        P += (U*(omega/2))@U.T
    return X, P


def modular_momentum(position, momentum):
    """Recover the true regional Gaussian modular form; never clip purity."""
    matrices = []
    for original in (position, momentum):
        if np.ma.isMaskedArray(original):
            raise ValueError("unmasked original real covariance matrices required")
        value = np.asarray(original)
        if (value.dtype.kind not in "fiu" or value.ndim != 2
                or value.shape[0] != value.shape[1] or not 1 <= len(value) <= 500
                or not np.isfinite(value).all() or np.max(abs(value)) > 1e100
                or np.max(abs(value-value.T)) > 2e-13*np.max(abs(value))):
            raise ValueError("finite real symmetric covariance matrices required")
        matrices.append((value.astype(float)+value.T.astype(float))/2)
    X, P = matrices
    if X.shape != P.shape:
        raise ValueError("matching regional covariance shapes required")
    values, vectors = np.linalg.eigh(X)
    if min(values) <= 0:
        raise ValueError("positive regional position covariance required")
    root = (vectors*np.sqrt(values))@vectors.T
    D = root@P@root
    D = (D+D.T)/2
    nu2, U = np.linalg.eigh(D)
    residual = float(np.linalg.norm(D@U-U*nu2, 2))
    gap = float(min(nu2)-.25)
    if gap < MIN_PURITY_GAP or residual > gap/20:
        raise ValueError("regional vacuum purity gaps unresolved at declared execution precision")
    nu = np.sqrt(nu2)
    # nu-.5=(nu^2-.25)/(nu+.5) preserves the small positive occupation.
    occupations = (nu2-.25)/(nu+.5)
    epsilon = np.log1p(1/occupations)
    hp = root@(U*(epsilon/nu))@U.T@root
    entropy = float(np.sum(np.log1p(occupations)+occupations*epsilon))
    if not np.isfinite(hp).all() or not math.isfinite(entropy):
        raise ValueError("finite modular charges required in execution precision")
    return {"momentum_form": (hp+hp.T)/2, "entropy": entropy,
            "minimum_symplectic_gap_squared": gap, "eigen_residual_norm": residual,
            "root": root, "nu_squared": nu2, "eigenvectors": U}


def charge(position, momentum, direction):
    result = modular_momentum(position, momentum)
    if np.ma.isMaskedArray(direction):
        raise ValueError("an original unmasked real momentum direction required")
    v = np.asarray(direction)
    if (v.dtype.kind not in "fiu" or v.shape != (len(result["momentum_form"]),)
            or not np.isfinite(v).all()):
        raise ValueError("a finite original real regional momentum direction required")
    # Normalize first so an overflowing requested charge produces an
    # explicit refusal, without a numpy warning or a spurious finite value.
    scale = float(np.max(abs(v)))
    direction_unit = v/scale if scale else v.astype(float)
    value = float(direction_unit@result["momentum_form"]@direction_unit/2)*scale*scale
    if not math.isfinite(value):
        raise ValueError("the regional modular charge overflowed execution precision")
    return value


def comparison(spacing=F(2, 7), box_length=F(16), radius=F(1), *, secants=False):
    spacing, length, radius = map(_positive, (spacing, box_length, radius))
    X, P, sites, M = covariance(spacing, length, radius)
    modular = modular_momentum(X, P)
    coordinates = float(spacing)*(sites-M//2)
    f = (float(radius)**2-np.sum(coordinates**2, axis=1))/(2*float(radius))
    rows = []
    for name, profile, coefficient in (("constant", np.ones(len(sites)), 15),
            ("linear_x", coordinates[:, 0]/float(radius), 105)):
        direction = float(spacing)**1.5*profile
        actual = float(direction@modular["momentum_form"]@direction/2)
        ideal = 4*math.pi**2*float(radius)**4/coefficient
        row = {"profile": name, "modular_entropy_derivative": actual,
            "conformal_continuum_stress_integral": ideal,
            "relative_discrepancy": actual/ideal-1,
            "same_grid_conformal_quadrature": float(math.pi*float(spacing)**3*np.sum(f*profile**2))}
        if secants:
            w = modular["eigenvectors"].T@modular["root"]@direction
            negative_radius = 1/float(np.sum(w*w/(modular["nu_squared"]-.25)))
            step = 1e-6/float(radius)**4
            if step >= negative_radius/20:
                raise ValueError("requested entropy secant exceeds the resolved regional state domain")
            # Both signs are physical REGIONAL Gaussian states. Only the
            # positive sign is asserted to arise from full-vacuum noise.
            plus = modular_momentum(X, P+step*np.outer(direction, direction))
            minus = modular_momentum(X, P-step*np.outer(direction, direction))
            row["entropy_secant_control"] = {"step": step,
                "positive_noise_slope": (plus["entropy"]-modular["entropy"])/step,
                "regional_negative_noise_slope": (modular["entropy"]-minus["entropy"])/step,
                "negative_step_below_one_twentieth_numerical_positivity_radius": True}
        rows.append(row)
    return {"proper_spacing": str(spacing), "proper_box_length": str(length),
        "proper_ball_radius": str(radius), "regional_sites": len(sites),
        "full_dirichlet_sites": (M-1)**3, "minimum_frequency": 2*math.sqrt(3)/float(spacing)*math.sin(math.pi/(2*M)),
        "minimum_regional_symplectic_gap_squared": modular["minimum_symplectic_gap_squared"],
        "modular_eigen_residual_norm": modular["eigen_residual_norm"],
        "bare_region_entropy_nats": modular["entropy"],
        "pixel_union_volume": len(sites)*float(spacing)**3,
        "pixel_boundary_Hausdorff_distance_upper": math.sqrt(3)*float(spacing)/2,
        "proper_lattice_spacing_isotropic": True, "all_symplectic_modes_retained": True,
        "profiles": rows}


def run():
    controls = []
    for spacing in (F(2, 3), F(1, 2), F(2, 5), F(1, 3), F(2, 7), F(1, 4)):
        try:
            controls.append(comparison(spacing, F(16), secants=spacing == F(2, 7)))
        except ValueError as error:
            # A refusal is part of the retained precision domain, never a
            # clipped or omitted-mode successful continuum comparison.
            if "purity gaps unresolved" not in str(error):
                raise
            controls.append({"proper_spacing": str(spacing), "status": "REFUSED", "reason": str(error)})
    independent_x, independent_p = direct_covariance()
    image_x, image_p, _, _ = covariance(F(2, 7), F(8))
    maximum = max(float(np.max(abs(independent_x-image_x))), float(np.max(abs(independent_p-image_p))))
    if maximum > 2e-12:
        raise ArithmeticError("independent Dirichlet sine vacuum disagrees with image covariance")
    nested = []
    for radius in (F(4, 5), F(1), F(6, 5)):
        X, P, sites, _ = covariance(F(2, 5), F(8), radius)
        nested.append({"radius": str(radius), "sites": len(sites),
                       "constant_coherent_relative_entropy_per_amplitude_squared":
                           charge(X, P, np.full(len(sites), float(F(2, 5))**1.5))})
    if any(a["constant_coherent_relative_entropy_per_amplitude_squared"] >
           b["constant_coherent_relative_entropy_per_amplitude_squared"] for a, b in zip(nested, nested[1:])):
        raise ArithmeticError("same-vacuum nested-region relative entropy decreased")
    report = {"schema": "oph.native_equilibration.ball_modular.v1", "proper_frame": proper_frame(),
        "exact_ball_identities": exact_ball_identities(), "grid_controls": controls,
        "infrared_controls": [comparison(F(2, 5), length) for length in (F(4), F(8), F(12), F(16))],
        "independent_Dirichlet_sine_covariance_maximum_difference": maximum,
        "same_vacuum_nested_ball_relative_entropy_controls": nested,
        "physical_dictionary": {"field": "one massless real scalar; unrestricted Gaussian vacuum",
            "units": "hbar=c=1 in the declared proper frame",
            "Hamiltonian": "H_tau=(p.p+q.(-Delta_a).q)/2; q=a**(3/2)phi,p=a**(3/2)pi",
            "state_direction": "all-outcome Gaussian momentum displacement; P_A(epsilon)=P_A+epsilon*v_A*v_A.T",
            "independent_stress_increment": "delta T_tau_tau(y)=profile(y)**2/2; delta(phi**2)=0 at tau=0",
            "continuum_kernel_normalization": "2*pi, fixed independently by conformal vacuum modular theorem",
            "cutoff": "same proper spacing in all three Cartesian directions; finite cubic dispersion is not rotation invariant",
            "region": "sites strictly inside proper Euclidean ball; pixel error retained",
            "boundary": "Dirichlet at walls of a proper cubic box",
            "selected_native_vacuum_or_joint_geometry_response": False},
        "scope": {"finite_Gaussian_partial_trace_and_first_law": True,
            "coherent_relative_entropy_identity_and_nested_monotonicity_exact": True,
            "numerical_grid_and_box_comparison_only": True,
            "finite_regulator_rotation_invariant": False,
            "interval_arithmetic_or_proven_continuum_error_rate": False,
            "native_smooth_tangent_modular_identification": False,
            "renormalized_generalized_entropy_or_area_density": False,
            "joint_state_geometry_stationarity_or_Einstein_dynamics": False,
            "massive_full_matter_or_calibrated_particle_spectrum": False},
        "entropy_secant_replay_absolute_tolerance": ENTROPY_SECANT_REPLAY_TOLERANCE,
        "references": ["https://arxiv.org/abs/1007.1813", "https://arxiv.org/abs/1102.0440", "https://arxiv.org/abs/1505.04753"],
        "source_bindings": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}}
    # This also refuses any unanticipated nonfinite scalar in the receipt.
    json.dumps(report, allow_nan=False)
    return report


def compare_ball_receipt(expected, actual, path="root"):
    """Exact keys, types, pins and discrete claims; stated floating tolerance."""
    if type(expected) is not type(actual):
        raise ValueError(f"receipt type mismatch at{path}")
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError(f"receipt keys mismatch at{path}")
        if "minimum_regional_symplectic_gap_squared" in expected:
            for row in (expected, actual):
                gap = row["minimum_regional_symplectic_gap_squared"]
                residual = row["modular_eigen_residual_norm"]
                if (type(gap) is not float or type(residual) is not float
                        or not math.isfinite(gap) or not math.isfinite(residual)
                        or gap < MIN_PURITY_GAP or not 0 <= residual <= gap/20):
                    raise ValueError(f"receipt purity/precision domain failed at{path}")
        for key in expected:
            compare_ball_receipt(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f"receipt length mismatch at{path}")
        for i, (a, b) in enumerate(zip(expected, actual)):
            compare_ball_receipt(a, b, f"{path}[{i}]")
    elif isinstance(expected, float):
        key = path.rsplit(".", 1)[-1]
        nonnegative = {"modular_entropy_derivative", "bare_region_entropy_nats",
            "constant_coherent_relative_entropy_per_amplitude_squared",
            "independent_Dirichlet_sine_covariance_maximum_difference"}
        positive = {"minimum_frequency", "pixel_union_volume", "pixel_boundary_Hausdorff_distance_upper",
            "conformal_continuum_stress_integral", "same_grid_conformal_quadrature", "step",
            "positive_noise_slope", "regional_negative_noise_slope"}
        if ((key in nonnegative and (expected < 0 or actual < 0))
                or (key in positive and (expected <= 0 or actual <= 0))
                or (key == "independent_Dirichlet_sine_covariance_maximum_difference"
                    and (expected > 2e-12 or actual > 2e-12))
                or (key == "relative_discrepancy" and (expected < -1 or actual < -1))):
            raise ValueError(f"receipt physical numeric domain failed at{path}")
        tolerance = ENTROPY_SECANT_REPLAY_TOLERANCE if key in {
            "positive_noise_slope", "regional_negative_noise_slope"} else 2e-10
        if key in {"entropy_secant_replay_absolute_tolerance", "step"}:
            tolerance = 0
        if key in {"minimum_regional_symplectic_gap_squared", "modular_eigen_residual_norm"}:
            tolerance = 5e-15
        if (not math.isfinite(expected) or not math.isfinite(actual)
                or abs(expected-actual) > tolerance):
            raise ValueError(f"floating replay mismatch at{path}")
    elif expected != actual:
        raise ValueError(f"receipt mismatch at{path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    actual = run()
    if args.check:
        retained = strict_json.loads(args.check.read_text())
        compare_ball_receipt(retained, actual)
        print("Proper-ball Gaussian modular and independent geometry comparison verified.")
    elif args.write:
        args.write.write_text(canonical(actual))
    else:
        print(canonical(actual), end="")


if __name__ == "__main__":
    main()
