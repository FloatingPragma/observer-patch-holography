"""Independent exact moment and finite-density audit of a quantum initial state.

No quantum or state producer is imported. Gaussian/Wick and simplex factors
are reconstructed with Fraction arithmetic; the cone is read independently
from pinned Lean data. Illustrative a=0 metric densities are rebuilt from
simplex monomials, not element quadrature. The cited analytic paper proves
the Hilbert/operator-domain statements; this finite replay does not.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction as Q
import hashlib
from itertools import combinations, product
import json
from math import factorial, prod
from pathlib import Path
import re

import mpmath
import numpy as np
from scipy.linalg import null_space
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/"runtime/whitney_quantum_state_receipt.json"
SCOPE = "SELECTED_NEUTRAL_GAUSSIAN_INITIAL_STATE__EXACT_MATTER_MOMENTS__NO_QUANTUM_TIME_HISTORY"
PIN_PATHS = {
    "Lean/Screen/SeamCurrentEdge30Moment.lean", "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean", "code/electromagnetism/verify_cone_whitney_bridge.py",
    "code/electromagnetism/whitney_interacting_quantum.py",
    "code/electromagnetism/whitney_quantum_state.py",
    "code/electromagnetism/verify_whitney_quantum_state.py",
    "code/electromagnetism/test_whitney_quantum_state.py",
}
PARAMETERS = {"sigma": "1/2", "charge": "1/4", "mass_squared": "1/2", "quartic": "1/4", "hbar": "1"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path=OUTPUT):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result
    def constant(_):
        raise ValueError("nonfinite JSON constant")
    def finite_float(value):
        result = float(value)
        require(np.isfinite(result), "nonfinite JSON number")
        return result
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs,
                      parse_constant=constant, parse_float=finite_float)


def object_keys(value, keys, name):
    require(type(value) is dict and set(value) == keys, name)


def exact_equal(actual, expected, name):
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), name)


def rational(value):
    require(type(value) is str, "rational string required")
    try:
        result = Q(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("invalid rational") from error
    require(str(result) == value, "canonical rational")
    return result


def numeric(value, name):
    require(type(value) in (int, float) and np.isfinite(value), "finite numeric "+name)
    return float(value)


def simplex_moment(powers):
    return Q(6*prod(factorial(p) for p in powers), factorial(sum(powers)+3))


def standard_gaussian_even_moment(degree):
    """Integration-by-parts recurrence for a unit-variance real Gaussian."""
    require(type(degree) is int and degree >= 0 and degree % 2 == 0,
            "nonnegative even Gaussian degree")
    return Q(prod(range(1, degree, 2)))


def exact_geometry():
    text = (ROOT/"Lean/Screen/SeamCurrentEdge30Moment.lean").read_text(encoding="utf-8")
    text = text.split("def portVector", 1)[1].split("theorem portVector_positivePort", 1)[0]
    phi = (1+sp.sqrt(5))/2
    values = {"0": sp.Integer(0), "1": sp.Integer(1), "-1": sp.Integer(-1), "φ": phi, "-φ": -phi}
    boundary = [sp.Matrix([values[x.strip()] for x in row.split(",")])
                for row in re.findall(r"!\[([^\[\]]+)\]", text)]
    require(len(boundary) == 12 and len({tuple(v) for v in boundary}) == 12, "exact vertex census")
    adjacency = {pair for pair in combinations(range(12), 2)
                 if sp.simplify((boundary[pair[0]]-boundary[pair[1]]).dot(boundary[pair[0]]-boundary[pair[1]])) == 4}
    text = (ROOT/"Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    text = text.split("def orientedFaces", 1)[1].split("def faceEdges", 1)[0]
    faces = [tuple(map(int, row)) for row in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", text)]
    geometric_faces = {face for face in combinations(range(12), 3)
                       if all(edge in adjacency for edge in combinations(face, 2))}
    require(len(adjacency) == 30 and len(faces) == 20 and {tuple(sorted(f)) for f in faces} == geometric_faces,
            "exact face adjacency")
    text = (ROOT/"Lean/Screen/SeamCurrentCarrierQuotient.lean").read_text(encoding="utf-8")
    ends = []
    for name in ("seamLeft", "seamRight"):
        row = text.split("def "+name, 1)[1].split("![", 1)[1].split("]", 1)[0]
        ends.append([int(v) for v in row.split(",")])
    boundary_edges = list(zip(*ends, strict=True))
    require(len(boundary_edges) == 30 and {tuple(sorted(edge)) for edge in boundary_edges} == adjacency,
            "exact oriented edge census")
    volumes = [sp.simplify(abs(sp.Matrix.hstack(*(boundary[i] for i in face)).det())/6) for face in faces]
    require(all(sp.simplify(v-(3+sp.sqrt(5))/6) == 0 for v in volumes), "exact tetrahedron volumes")
    total = sp.expand(sum(volumes))
    sqrt_coefficient = total.coeff(sp.sqrt(5))
    rational_part = sp.simplify(total-sqrt_coefficient*sp.sqrt(5))
    volume = [Q(str(rational_part)), Q(str(sqrt_coefficient))]
    xyz = np.array([[0.0]*3]+[[float(x) for x in v] for v in boundary])
    edges = [(0, i) for i in range(1, 13)]+[(i+1, j+1) for i, j in boundary_edges]
    tets = [(0,)+tuple(i+1 for i in face) for face in faces]
    return volume, xyz, edges, tets


def expected_state():
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


def expected_observables(volume):
    sigma = Q(PARAMETERS["sigma"])
    # One real centered Gaussian has variance sigma^2 and fourth moment
    # 3*sigma^4. The two independent components give 3+2+3=8.
    real_second = standard_gaussian_even_moment(2)
    real_fourth = standard_gaussian_even_moment(4)
    complex_second = 2*real_second
    complex_fourth = 2*real_fourth+2*real_second**2
    second = simplex_moment((2, 0, 0, 0))
    fourth = simplex_moment((4, 0, 0, 0))
    mixed = simplex_moment((2, 2, 0, 0))
    sum_second, sum_fourth = 4*second, 4*fourth+2*6*mixed
    l2 = complex_second*sigma**2*sum_second
    l4 = complex_fourth*sigma**4*sum_fourth
    coefficients = {"matter_l2": l2, "matter_l4": l4,
        "mass_potential": Q(PARAMETERS["mass_squared"])*l2,
        "quartic_potential": Q(PARAMETERS["quartic"])*l4/2}
    return {
        "normalized_gaussian_norm_squared": "1", "one_real_coordinate_variance": str(sigma**2),
        "complex_gaussian_wick_factors": {"second": str(complex_second), "fourth": str(complex_fourth)},
        "simplex_moments_divided_by_volume": {"lambda_i_squared": str(second),
            "lambda_i_fourth": str(fourth), "lambda_i_squared_lambda_j_squared": str(mixed),
            "sum_lambda_squared": str(sum_second), "sum_lambda_squared_squared": str(sum_fourth)},
        "volume_in_Qsqrt5": [str(x) for x in volume],
        "coefficients_per_volume": {name: str(value) for name, value in coefficients.items()},
        "expectations_in_Qsqrt5": {name: [str(value*x) for x in volume] for name, value in coefficients.items()},
        "interpretation": {"matter_l2": "E[integral_Omega |Psi_a|^2 dx]",
            "matter_l4": "E[integral_Omega |Psi_a|^4 dx]",
            "mass_potential": "E[m^2 integral_Omega |Psi_a|^2 dx]",
            "quartic_potential": "E[(g/2) integral_Omega |Psi_a|^4 dx]"},
        "derivation": "condition on a; independent circular complex Gaussian nodal matter removes every unit dressing phase; exact simplex moments then integrate the Wick factors",
        "not_evaluated": ["kinetic Hamiltonian expectation", "covariant-gradient expectation", "magnetic-energy value", "full energy expectation", "time-dependent observables"],
    }


def exact_real(value):
    """Retain the supplied real component before a common array dtype exists."""
    require(not isinstance(value, (bool, np.bool_)), "finite real replay scalar")
    try:
        if isinstance(value, (int, np.integer)):
            return Q(int(value))
        if isinstance(value, Q):
            return Q(int(value.numerator), int(value.denominator))
        if isinstance(value, Decimal):
            return Q(value)
        require(isinstance(value, (float, np.floating)), "finite real replay scalar")
        return Q(*value.as_integer_ratio())
    except (ValueError, OverflowError) as error:
        raise ValueError("finite real replay scalar") from error


def real_components(value, shape):
    require(not np.ma.isMaskedArray(value), "unmasked replay data")
    raw = np.asarray(value, dtype=object)
    require(raw.shape == shape, "replay data shape")
    return np.array([exact_real(v) for v in raw.flat], dtype=object).reshape(shape)


def complex_components(value):
    require(not np.ma.isMaskedArray(value), "unmasked replay matter")
    raw = np.asarray(value, dtype=object)
    require(raw.shape == (13,), "13 complex replay matter coordinates")
    return [(exact_real(v.real), exact_real(v.imag))
            if isinstance(v, (complex, np.complexfloating)) else (exact_real(v), Q(0))
            for v in raw]


def require_original_gradient(values, edges):
    """Return its mean-zero potential after exact original edge-cycle closure.

    A numerical Coulomb residual cannot decide whether a small radiative
    field is zero. The monomial replay supports only configurations whose
    edge cochain is a vertex gradient; it does not evaluate a nonzero
    Coulomb field's magnetic energy or charged dressing phases.
    """
    values = real_components(values, (42,))
    require(len(edges) == 42 and all(len(edge) == 2 and
            all(type(i) is int and 0 <= i < 13 for i in edge) and edge[0] != edge[1]
            for edge in edges), "original gradient incidence")
    potential = [None]*13
    potential[0] = Q(0)
    for _ in range(13):
        for value, (left, right) in zip(values, edges, strict=True):
            if potential[left] is not None and potential[right] is None:
                potential[right] = potential[left]+value
            elif potential[right] is not None and potential[left] is None:
                potential[left] = potential[right]-value
        if all(value is not None for value in potential):
            break
    require(all(value is not None for value in potential), "original gradient connected incidence")
    require(all(potential[right]-potential[left] == value
                for value, (left, right) in zip(values, edges, strict=True)),
            "original edge field must be an exact vertex gradient for the a=0 replay")
    mean = sum(potential, Q(0))/13
    return [value-mean for value in potential]


def replay_context(*components):
    # Four powers of the largest scale cover the field-dependent condition
    # numbers. The full exponent span also matters: an O(1) gauge tangent
    # can cancel and leave an originally supplied 1e-100 physical velocity.
    # Numerator/denominator lengths retain exact near-one differences too,
    # which an exponent-only budget would miss. Bit lengths avoid Python's
    # decimal integer-string limit. Keep another 270 guard bits (80 digits).
    # This remains a numerical replay budget, not an interval certificate.
    lengths = [(abs(v.numerator).bit_length(), v.denominator.bit_length())
               for values in components for v in values if v]
    upper = max([0]+[numerator-denominator+1 for numerator, denominator in lengths])
    lower = min([0]+[numerator-denominator-1 for numerator, denominator in lengths])
    original_bits = max([0]+[max(pair) for pair in lengths])
    mp = mpmath.mp.clone()
    mp.prec = 270+3*upper+max(upper-lower, original_bits)
    return mp


def mp_real(mp, value):
    return mp.mpf(value.numerator)/value.denominator


def moment_system(mp, psi, xyz, edges, tets, *, potential=False):
    """Original-input simplex moments, before any binary64 Gram assembly.

    psi and xyz already belong to this private high-precision context. The
    geometry consists of the supplied vertices, not a repaired rounded Gram.
    """
    require(len(edges) == 42 and len(tets) == 20 and xyz.rows == 13 and xyz.cols == 3,
            "fixed cone replay geometry")
    charge = mp.mpf(1)/4
    metric, mass, energy = mp.matrix(68), mp.matrix(42), mp.mpf(0)
    for tet in tets:
        require(len(tet) == 4 and len(set(tet)) == 4 and
                all(type(i) is int and 0 <= i < 13 for i in tet), "tetrahedron vertices")
        points = mp.matrix([[xyz[i, j] for j in range(3)] for i in tet])
        affine = mp.matrix([[1]+list(points[i, :]) for i in range(4)])
        try:
            gradients = mp.inverse(affine)[1:, :].T
        except ZeroDivisionError as error:
            raise ValueError("nondegenerate replay tetrahedron") from error
        volume = abs(mp.det(mp.matrix([[points[i, j]-points[0, j]
                                      for j in range(3)] for i in range(1, 4)])))/6
        require(volume > 0, "positive replay tetrahedron volume")
        gradient_gram = gradients*gradients.T
        edge_columns, columns = [], []
        for e, (u, v) in enumerate(edges):
            if u in tet and v in tet:
                i, j = tet.index(u), tet.index(v)
                edge_columns.append((e, i, j))
                columns.append((e, mp.j*charge*(psi[u]-psi[v]), (i, j)))
        for i, vertex in enumerate(tet):
            columns.extend(((42+vertex, 1, (i,)), (55+vertex, mp.j, (i,))))
        for e, i, j in edge_columns:
            for f, k, l in edge_columns:
                mass[e, f] += volume/20*((1+int(i == k))*gradient_gram[j, l]
                    -(1+int(i == l))*gradient_gram[j, k]
                    -(1+int(j == k))*gradient_gram[i, l]
                    +(1+int(j == l))*gradient_gram[i, k])
        for col, value, powers in columns:
            for other, other_value, other_powers in columns:
                degree = powers+other_powers
                integral = volume*mp_real(mp, simplex_moment(tuple(degree.count(i) for i in range(4))))
                metric[col, other] += 2*mp.re(mp.conj(value)*other_value)*integral
        if potential:
            local = [psi[i] for i in tet]
            energy += volume*mp.fsum(mp.re(mp.conj(local[i])*local[j])*gradient_gram[i, j]
                                    for i, j in product(range(4), repeat=2))
            for degree, coefficient in ((2, mp.mpf(1)/2), (4, mp.mpf(1)/8)):
                for indices in product(range(4), repeat=degree):
                    term = mp.fprod(mp.conj(local[i]) if n % 2 == 0 else local[i]
                                    for n, i in enumerate(indices))
                    energy += coefficient*volume*mp.re(term)*mp_real(
                        mp, simplex_moment(tuple(indices.count(i) for i in range(4))))
    metric[:42, :42] += mass
    return metric, mass, energy


def positive_factor(mp, matrix, name):
    """A resolved Cholesky factor; determinant sign cannot establish SPD."""
    require(matrix.rows == matrix.cols and all(mp.isfinite(x) for x in matrix), name)
    scale = max((abs(x) for x in matrix), default=mp.mpf(0))
    require(scale > 0, name)
    tolerance = scale*mp.power(10, 20-mp.dps)
    require(max(abs(x) for x in matrix-matrix.T) <= tolerance, name+" symmetry")
    try:
        return mp.cholesky((matrix+matrix.T)/2, tol=tolerance)
    except ValueError as error:
        raise ValueError(name+" positive definite at replay precision") from error


def solve_factor(mp, factor, right):
    """Reuse a Cholesky factor for every right-hand side."""
    result = mp.matrix(right.rows, right.cols)
    for col in range(right.cols):
        intermediate = mp.matrix(right.rows, 1)
        for i in range(right.rows):
            intermediate[i] = (right[i, col]-mp.fsum(factor[i, j]*intermediate[j]
                                                    for j in range(i)))/factor[i, i]
        for i in range(right.rows-1, -1, -1):
            result[i, col] = (intermediate[i]-mp.fsum(factor[j, i]*result[j, col]
                                                     for j in range(i+1, right.rows)))/factor[i, i]
    return result


def reduced_moments(mp, metric, vertical, section):
    inertia, coupling = vertical.T*metric*vertical, vertical.T*metric*section
    factor = positive_factor(mp, inertia, "gauge inertia")
    eta_map = -solve_factor(mp, factor, coupling)
    gamma = section.T*metric*section+coupling.T*eta_map
    root = positive_factor(mp, gamma, "reduced kinetic metric")
    # Symmetrize only the working-precision multiplication residual already
    # checked above. No eigenvalue clipping or binary64 positivity repair.
    return (gamma+gamma.T)/2, eta_map, mp.fsum(mp.log(root[i, i]) for i in range(root.rows))


def reported(value, name):
    result = float(value)
    require(np.isfinite(result) and (value == 0 or abs((result-value)/value) <= 1e-12),
            name+" outside reliable binary64 reporting range")
    return result


def independent_log_rho(psi, xyz, edges, tets):
    """Original-input moments and high-precision constrained determinant."""
    parts = complex_components(psi)
    vertices = real_components(xyz, (13, 3))
    mp = replay_context([x for pair in parts for x in pair], vertices.flat)
    matter = [mp.mpc(mp_real(mp, re), mp_real(mp, im)) for re, im in parts]
    points = mp.matrix([[mp_real(mp, x) for x in row] for row in vertices])
    metric, mass, _ = moment_system(mp, matter, points, edges, tets)
    d = mp.matrix(42, 13)
    for edge, (i, j) in enumerate(edges):
        d[edge, i], d[edge, j] = -1, 1
    # The small fixed geometric nullspace is well-conditioned; lift its
    # supplied binary64 basis before the field-dependent matrix operations.
    transverse = null_space(np.array((d.T*mass).tolist(), dtype=float))
    require(transverse.shape == (42, 30), "Coulomb dimension")
    section = mp.matrix(68, 56)
    section[:42, :30], section[42:, 30:] = mp.matrix(transverse.tolist()), mp.eye(26)
    basis = mp.matrix([[-1]*12]+np.eye(12, dtype=int).tolist())
    vertical = mp.matrix(68, 13)
    vertical[:42, :] = d
    for i, value in enumerate(matter):
        vertical[42+i, i], vertical[55+i, i] = -mp.im(value)/4, mp.re(value)/4
    _, _, log_rho = reduced_moments(mp, metric, vertical*basis, section)
    return reported(log_rho, "log density")


def verify(packet):
    object_keys(packet, {"schema", "scope", "run_id", "source_pins", "parameters", "state", "exact_observables", "numerical_amplitudes"}, "packet schema")
    require(packet["schema"] == "oph.whitney_quantum_initial_state.v1" and packet["scope"] == SCOPE, "scope")
    require(packet["run_id"] == "whitney-neutral-gaussian-initial-state-2026-09-06-v1", "run identity")
    exact_equal(packet["parameters"], PARAMETERS, "exact state parameters")
    object_keys(packet["source_pins"], PIN_PATHS, "source pin census")
    for path in PIN_PATHS:
        require(packet["source_pins"][path] == hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), "source pin: "+path)
    exact_equal(packet["state"], expected_state(), "state normalization/density/domain contract")
    state = packet["state"]
    dimension = state["real_configuration_dimension"]
    require(dimension == state["radiative_coordinates"]+2*state["complex_matter_coordinates"], "real dimension")
    require(2*rational(state["rho_power"])+1 == 0, "half-density measure cancellation")
    require(rational(state["gaussian_prefactor_power"]) == -Q(dimension, 4), "Gaussian amplitude normalization")
    volume, xyz, edges, tets = exact_geometry()
    expected = expected_observables(volume)
    exact_equal(packet["exact_observables"], expected, "exact Gaussian/Wick/simplex observables")
    numerical = packet["numerical_amplitudes"]
    object_keys(numerical, {"scope", "basis", "samples"}, "numerical amplitude schema")
    require(numerical["scope"] == "illustrative floating-point amplitudes from element quadrature and Schur determinants; no certified amplitude enclosure", "numeric amplitude scope")
    require(numerical["basis"] == "same Euclidean-orthonormal Coulomb frame at orders4,5,6; all30 radiative sample coordinates are zero", "numeric amplitude frame")
    require(type(numerical["samples"]) is list and len(numerical["samples"]) == 3, "amplitude sample census")
    diagnostics = []
    sigma = float(Q(PARAMETERS["sigma"]))
    for row, (name, slot) in zip(numerical["samples"], (("origin", None), ("center-real", 30), ("center-imaginary", 43)), strict=True):
        object_keys(row, {"id", "coordinates_exact", "log_g_numeric", "evaluations"}, "amplitude sample schema")
        require(row["id"] == name, "amplitude sample identity")
        target = ["0"]*56
        if slot is not None:
            target[slot] = "1/2"
        exact_equal(row["coordinates_exact"], target, "sample exact coordinates")
        q = np.array([float(rational(x)) for x in row["coordinates_exact"]])
        log_g = -dimension/4*np.log(2*np.pi*sigma*sigma)-q@q/(4*sigma*sigma)
        require(abs(numeric(row["log_g_numeric"], "log g")-log_g) < 2e-12, "Gaussian log amplitude")
        log_rho = independent_log_rho(q[30:43]+1j*q[43:], xyz, edges, tets)
        log_f = log_g-log_rho/2
        require(type(row["evaluations"]) is list and len(row["evaluations"]) == 3, "quadrature sample census")
        differences = []
        for result, order in zip(row["evaluations"], (4, 5, 6), strict=True):
            object_keys(result, {"quadrature_order", "log_rho_numeric", "log_f_numeric"}, "quadrature result schema")
            require(type(result["quadrature_order"]) is int and result["quadrature_order"] == order, "exact quadrature order")
            rho_difference = abs(numeric(result["log_rho_numeric"], "log rho")-log_rho)
            f_difference = abs(numeric(result["log_f_numeric"], "log f")-log_f)
            require(rho_difference < 2e-10 and f_difference < 2e-10, "independent simplex metric/amplitude")
            differences.append(max(rho_difference, f_difference))
        diagnostics.append({"id": name, "log_rho_numeric": log_rho, "log_f_numeric": float(log_f),
                            "max_quadrature_replay_difference": max(differences)})
    require(abs(diagnostics[1]["log_f_numeric"]-diagnostics[2]["log_f_numeric"]) < 2e-12,
            "global-phase-equivalent amplitude samples")
    return {"accepted": True, "scope": SCOPE, "gaussian_sigma": "1/2",
        "real_configuration_dimension": 56, "gaussian_norm_squared": "1",
        "matter_l2_coefficient": expected["coefficients_per_volume"]["matter_l2"],
        "matter_l4_coefficient": expected["coefficients_per_volume"]["matter_l4"],
        "volume_in_Qsqrt5": expected["volume_in_Qsqrt5"],
        "initial_state_constructed": True, "initial_observables_provided": True,
        "quantum_time_history": False, "observer_history": False, "physical_comparison": False,
        "analytic_operator_domain_proved_by_numeric_replay": False,
        "numeric_diagnostics": diagnostics}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), sort_keys=True))
