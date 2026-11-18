"""Independent exact replay of the full56D trial-history certificate.

No history producer or quantum coefficient evaluator is imported. Geometry is
reconstructed in Q(sqrt(5)); scalar envelopes, Gaussian radial moments and
sampled phases use exact arithmetic. The analytic norm comparison is proved
in the cited paper, not by floating-point checks or a five-coordinate model.
"""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import combinations
import json
from math import factorial, isfinite, isqrt
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent/"runtime/whitney_quantum_history_receipt.json"
SCOPE = "FULL_56D_NEUTRAL_POTENTIAL_PHASE_TRIAL__GLOBAL_NORM_BOUND__DECLARED_INPUTS"
PINS = {
    "Lean/Screen/SeamCurrentEdge30Moment.lean", "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "code/electromagnetism/whitney_interacting_quantum.py",
    "code/electromagnetism/verify_whitney_quantum_state.py",
    "code/electromagnetism/runtime/whitney_quantum_state_receipt.json",
    "code/electromagnetism/whitney_quantum_history.py",
    "code/electromagnetism/verify_whitney_quantum_history.py",
    "code/electromagnetism/test_whitney_quantum_history.py",
}
PARAMETERS = {"sigma": "1/2", "charge": "1/4", "mass_squared": "1/2",
              "quartic": "1/4", "hbar": "1", "target_norm_error": "1/10"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def equal(actual, expected, message):
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), message)


def rational(value):
    require(type(value) is str, "canonical rational string required")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("invalid rational") from error
    require(str(answer) == value, "noncanonical rational")
    return answer


def load(path=OUTPUT):
    def pairs(items):
        obj = {}
        for key, value in items:
            require(key not in obj, "duplicate JSON key")
            obj[key] = value
        return obj
    def floating(value):
        result = float(value)
        require(isfinite(result), "nonfinite JSON number")
        return result
    def constant(_):
        raise ValueError("nonfinite JSON constant")
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs,
                      parse_float=floating, parse_constant=constant)


class Golden:
    """Exact a+b*sqrt(5), including comparison without decimal rounding."""
    def __init__(self, a=0, b=0):
        self.a, self.b = Q(a), Q(b)

    @staticmethod
    def cast(value):
        return value if isinstance(value, Golden) else Golden(value)

    def __add__(self, value):
        value = self.cast(value)
        return Golden(self.a+value.a, self.b+value.b)

    __radd__ = __add__

    def __neg__(self):
        return Golden(-self.a, -self.b)

    def __sub__(self, value):
        return self+-self.cast(value)

    def __rsub__(self, value):
        return self.cast(value)+-self

    def __mul__(self, value):
        value = self.cast(value)
        return Golden(self.a*value.a+5*self.b*value.b, self.a*value.b+self.b*value.a)

    __rmul__ = __mul__

    def __truediv__(self, value):
        value = self.cast(value)
        denominator = value.a**2-5*value.b**2
        require(denominator != 0, "zero exact field divisor")
        return self*Golden(value.a/denominator, -value.b/denominator)

    def sign(self):
        if self.b == 0:
            return (self.a > 0)-(self.a < 0)
        if self.a == 0 or (self.a > 0) == (self.b > 0):
            return 1 if self.b > 0 else -1
        difference = self.a**2-5*self.b**2
        return ((difference > 0)-(difference < 0))*(1 if self.a > 0 else -1)

    def pair(self):
        return [str(self.a), str(self.b)]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b, strict=True)), Golden())


def subtract(a, b):
    return [x-y for x, y in zip(a, b, strict=True)]


def determinant3(rows):
    a, b, c = rows
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


def inverse(matrix):
    n = len(matrix)
    work = [row.copy()+[Golden(int(i == j)) for j in range(n)] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((i for i in range(col, n) if work[i][col].sign()), None)
        require(pivot is not None, "singular exact simplex")
        work[col], work[pivot] = work[pivot], work[col]
        divisor = work[col][col]
        work[col] = [value/divisor for value in work[col]]
        for i in range(n):
            if i != col:
                factor = work[i][col]
                work[i] = [x-factor*y for x, y in zip(work[i], work[col], strict=True)]
    return [row[n:] for row in work]


def exact_geometry():
    source = (ROOT/"Lean/Screen/SeamCurrentEdge30Moment.lean").read_text(encoding="utf-8")
    source = source.split("def portVector", 1)[1].split("theorem portVector_positivePort", 1)[0]
    phi = Golden(Q(1, 2), Q(1, 2))
    vocabulary = {"0": Golden(), "1": Golden(1), "-1": Golden(-1), "φ": phi, "-φ": -phi}
    vertices = [[Golden()]*3]+[[vocabulary[x.strip()] for x in row.split(",")]
        for row in re.findall(r"!\[([^\[\]]+)\]", source)]
    require(len(vertices) == 13 and len({tuple(tuple(x.pair()) for x in row) for row in vertices}) == 13,
            "exact cone vertex census")
    source = (ROOT/"Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    source = source.split("def orientedFaces", 1)[1].split("def faceEdges", 1)[0]
    faces = [tuple(int(x)+1 for x in row) for row in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", source)]
    adjacency = {edge for edge in combinations(range(1, 13), 2)
        if (dot(subtract(vertices[edge[0]], vertices[edge[1]]), subtract(vertices[edge[0]], vertices[edge[1]]))-4).sign() == 0}
    geometric_faces = {face for face in combinations(range(1, 13), 3)
        if all(edge in adjacency for edge in combinations(face, 2))}
    require(len(faces) == 20 and len(adjacency) == 30 and {tuple(sorted(f)) for f in faces} == geometric_faces,
            "exact face and adjacency census")
    source = (ROOT/"Lean/Screen/SeamCurrentCarrierQuotient.lean").read_text(encoding="utf-8")
    ends = []
    for name in ("seamLeft", "seamRight"):
        row = source.split("def "+name, 1)[1].split("![", 1)[1].split("]", 1)[0]
        ends.append([int(x)+1 for x in row.split(",")])
    boundary_edges = list(zip(*ends, strict=True))
    require(len(boundary_edges) == 30 and {tuple(sorted(e)) for e in boundary_edges} == adjacency,
            "exact edge census")
    edges = [(0, i) for i in range(1, 13)]+boundary_edges
    for i, j in edges:
        diff = subtract(vertices[i], vertices[j])
        require((Golden(4)-dot(diff, diff)).sign() >= 0, "edge length envelope")
    elements = []
    for face in faces:
        tet = (0,)+face
        coords = [vertices[i] for i in tet]
        det = determinant3(coords[1:])
        volume = (det if det.sign() > 0 else -det)/6
        require((volume-Q(5, 6)).sign() >= 0, "simplex volume lower bound")
        inv = inverse([[Golden(1)]+row for row in coords])
        gradients = [[inv[j+1][i] for j in range(3)] for i in range(4)]
        require(all((Golden(4)-dot(g, g)).sign() >= 0 for g in gradients), "barycentric gradient envelope")
        elements.append((tet, volume, gradients))
    total_volume = sum((row[1] for row in elements), Golden())
    require((Golden(18)-total_volume).sign() >= 0, "cone volume upper bound")
    require(min(sum(i in tet for tet, _, _ in elements) for i in range(13)) >= 5, "node multiplicity")
    require(min(sum(i in tet and j in tet for tet, _, _ in elements) for i, j in edges) >= 2, "edge multiplicity")
    # Local affine-field and dominant-corner inequalities are proved in the
    # paper. These exact scalar comparisons bind their geometric constants.
    require(2*Q(5, 6)/(30*4) >= Q(1, 81), "Maxwell coercivity arithmetic")
    require(5*Q(5, 6)/(80*4) >= Q(1, 128), "dressed-mass coercivity arithmetic")
    return elements, total_volume


def sqrt_rational(value):
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    require(a*a == value.numerator and b*b == value.denominator, "exact square-root constant")
    return Q(a, b)


def record_polynomial(poly):
    return [{"powers": [x, y], "coefficient": str(value)} for (x, y), value in sorted(poly.items())]


def product(left, right):
    output = {}
    for first, a in left.items():
        for second, b in right.items():
            key = tuple(x+y for x, y in zip(first, second, strict=True))
            output[key] = output.get(key, Q(0))+a*b
    return output


def moment(dimension, power):
    # A different recurrence organization from the producer: integrate the
    # radial Gaussian density recursively before applying variance1/4.
    value = Q(1)
    for j in range(power, 0, -1):
        value *= Q(dimension+2*j-2, 4)
    return value


def expectation(poly):
    return sum((value*moment(30, x)*moment(26, y) for (x, y), value in poly.items()), Q(0))


def independent_bounds():
    volume, e, mass_squared, quartic = Q(18), Q(1, 4), Q(1, 2), Q(1, 4)
    maxwell_inverse, dressed_inverse, gradient = Q(81), Q(128), Q(2)
    velocity_a = sqrt_rational(maxwell_inverse)
    scalar0 = sqrt_rational(dressed_inverse/2)
    scalar1 = sqrt_rational(dressed_inverse*volume*maxwell_inverse*e**2)
    relative_G = sqrt_rational(8*volume*e**2)
    relative_I = sqrt_rational(8*volume*maxwell_inverse*e**2)
    ell0 = (68*relative_G*(scalar0+velocity_a)+12*relative_I)/2
    ell1 = (68*relative_G*(scalar1+e*velocity_a)+12*relative_I*e)/2
    inverse_metric = {(0, 0): maxwell_inverse+dressed_inverse,
                      (0, 1): 4*dressed_inverse*volume*e**2*maxwell_inverse}
    kinetic = product(inverse_metric, {(0, 0): ell0**2/2,
                                      (0, 1): ell1**2/2+4, (1, 0): Q(4)})
    curl_squared = 6*(2*gradient**2)**2
    magnetic_force = volume*curl_squared
    # Both affine-in-|a| factors in the gradient estimate are <=4(1+|a|).
    require(2*gradient <= 4 and e*9 <= 4 and e*(4+9) <= 4 and e**2*9 <= 4,
            "spatial derivative factor bounds")
    base = 4*volume*16
    ca_y = base+2*mass_squared*volume*e
    ca_y2 = 2*quartic*volume*e
    cp_y = base+2*mass_squared*volume
    cp_y3 = 2*quartic*volume
    grad_poly = {(1, 0): 4*magnetic_force**2, (0, 2): 4*ca_y**2,
        (2, 2): 4*base**2, (0, 4): 4*ca_y2**2, (0, 1): 3*cp_y**2,
        (2, 1): 3*base**2, (0, 3): 3*cp_y3**2}
    potential = {(1, 0): magnetic_force/2,
        (0, 1): 2*volume*(2*gradient)**2+mass_squared*volume,
        (1, 1): 2*volume*(9*e)**2, (0, 2): quartic*volume/2}
    force = product(inverse_metric, grad_poly)
    # Independently expand the kinetic expectation through low radial moments.
    a0, a1 = inverse_metric[(0, 0)], inverse_metric[(0, 1)]
    kinetic_upper = (ell0**2*(a0+a1*moment(26, 1))/2
        +ell1**2*(a0*moment(26, 1)+a1*moment(26, 2))/2
        +4*(a0*(moment(30, 1)+moment(26, 1))
             +a1*(moment(30, 1)*moment(26, 1)+moment(26, 2))))
    require(kinetic_upper == expectation(kinetic), "independent kinetic expansion")
    potential_upper = expectation(potential)
    return ({"inverse_metric": inverse_metric, "initial_kinetic": kinetic,
             "potential": potential, "potential_gradient_squared": grad_poly, "force_moment": force},
            {"kinetic_upper": kinetic_upper, "potential_upper": potential_upper,
             "total_energy_upper": kinetic_upper+potential_upper,
             "force_moment_upper": expectation(force)}, ell0, ell1)


def simplex_average(powers):
    answer = Q(6, factorial(sum(powers)+3))
    for power in powers:
        answer *= factorial(power)
    return answer


def sampled_potential(q, elements):
    require(len(q) == 56 and all(x == 0 for x in q[:30]), "exact sample has zero connection")
    answer = Golden()
    for tet, volume, gradients in elements:
        real = [q[30+i] for i in tet]
        imag = [q[43+i] for i in tet]
        real_gradient = [sum((value*grad[c] for value, grad in zip(real, gradients, strict=True)), Golden()) for c in range(3)]
        imag_gradient = [sum((value*grad[c] for value, grad in zip(imag, gradients, strict=True)), Golden()) for c in range(3)]
        scalar = {}
        for i in range(4):
            for j in range(4):
                powers = tuple(int(k == i)+int(k == j) for k in range(4))
                scalar[powers] = scalar.get(powers, Q(0))+real[i]*real[j]+imag[i]*imag[j]
        fourth = {}
        for left, a in scalar.items():
            for right, b in scalar.items():
                powers = tuple(x+y for x, y in zip(left, right, strict=True))
                fourth[powers] = fourth.get(powers, Q(0))+a*b
        mass = sum((c*simplex_average(p) for p, c in scalar.items()), Q(0))/2
        interaction = sum((c*simplex_average(p) for p, c in fourth.items()), Q(0))/8
        answer += volume*(dot(real_gradient, real_gradient)+dot(imag_gradient, imag_gradient)+mass+interaction)
    require(answer.sign() >= 0, "sample potential positivity")
    return answer


def parent_initial_state():
    path = ROOT/"code/electromagnetism/verify_whitney_quantum_state.py"
    spec = importlib.util.spec_from_file_location("_oph_initial_state_parent_for_history", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.verify(module.load(ROOT/"code/electromagnetism/runtime/whitney_quantum_state_receipt.json"))
    require(result["accepted"] is True and result["gaussian_sigma"] == "1/2"
            and result["real_configuration_dimension"] == 56 and result["quantum_time_history"] is False,
            "independent initial-state parent")


def verify(packet):
    fields = {"schema", "scope", "run_id", "source_pins", "parameters", "state", "geometry_envelope",
              "coefficient_bounds", "integration", "bounds", "error_certificate", "phase_configurations", "history"}
    require(type(packet) is dict and set(packet) == fields, "history schema")
    require(packet["schema"] == "oph.whitney_quantum_trial_history.v1" and packet["scope"] == SCOPE, "history scope")
    require(packet["run_id"] == "whitney-full56d-potential-phase-certified-trial-v1", "history identity")
    equal(packet["parameters"], PARAMETERS, "fixed declared parameters")
    require(type(packet["source_pins"]) is dict and set(packet["source_pins"]) == PINS, "source pin census")
    for path in PINS:
        require(packet["source_pins"][path] == hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), "source pin: "+path)
    equal(packet["state"], {"dimension": 56, "trial_formula": "v(t,q)=exp(-i*t*V(q)/hbar)*f_sigma(q)",
        "exact_evolution": "u(t)=exp(-i*t*H/hbar)*f_sigma",
        "initial_state_parent": "code/electromagnetism/runtime/whitney_quantum_state_receipt.json",
        "measure": "dmu=rho(q)dq; rho=sqrt(det(gamma))",
        "probability_law": "|v(t,q)|^2*dmu=N(0,sigma^2 I56); independent of t",
        "neutrality": "V and f_sigma are invariant under the residual global U(1)",
        "trial_history_computed": True, "exact_Hamiltonian_history_computed": False,
        "rigorous_global_norm_bound": True, "configuration_density_moves": False,
        "classical_five_coordinate_restriction": False, "observer_history": False,
        "physical_state_preparation": False, "empirical_comparison": False,
        "continuum_QFT": False, "ordinary_physics_benchmark": False}, "state and interpretation boundary")
    equal(packet["geometry_envelope"], {"volume_upper": "18", "tetrahedron_volume_lower": "5/6",
        "edge_length_squared_upper": "4", "barycentric_gradient_squared_upper": "4",
        "minimum_node_multiplicity": 5, "minimum_edge_multiplicity": 2,
        "dressed_mass_lower": "1/128", "Maxwell_mass_lower": "1/81"}, "geometry envelope")
    elements, volume = exact_geometry()
    polys, bounds, ell0, ell1 = independent_bounds()
    equal(packet["coefficient_bounds"], {
        "variables": "X=|a|^2,Y=|psi|^2 in orthonormal Coulomb coordinates",
        "log_density_gradient_upper": str(ell0)+"+"+str(ell1)+"*sqrt(Y)",
        "polynomials": {name: record_polynomial(poly) for name, poly in polys.items()}}, "derived coefficient envelopes")
    equal(packet["integration"], {"method": "exact independent scaled chi-square moments, dimensions30 and26",
        "radial_moment": "E[X^k]=sigma^(2k)*product(d+2j,j=0..k-1)",
        "global_envelopes": True, "Gaussian_tail_truncation": False,
        "Monte_Carlo_used": False, "quadrature_used_for_bound": False,
        "floating_point_used_for_bound": False}, "exact global integration scope")
    equal(packet["bounds"], {name: str(value) for name, value in bounds.items()}, "exact energy and force bounds")
    def error(time):
        return time*(bounds["total_energy_upper"]+bounds["kinetic_upper"]+time*time*bounds["force_moment_upper"]/6)
    target = Q(1, 100)
    horizon, power = Q(1), 0
    while error(horizon) > target:
        horizon /= 2
        power += 1
    equal(packet["error_certificate"], {"norm": "L2(Q,dvol_gamma)",
        "bound_formula": "||u(t)-v(t)||^2 <= ((Ebar+Kbar)*abs(t)+Bbar*abs(t)^3/6)/hbar",
        "trivial_norm_error_upper": "2", "target_norm_error": "1/10",
        "dyadic_horizon_power": power, "horizon": str(horizon),
        "horizon_squared_error_upper": str(error(horizon)),
        "previous_dyadic_squared_error_upper": str(error(2*horizon)),
        "coverage": "every real t with |t|<=horizon, not only stored samples",
        "time_units": "declared dimensionless Hamiltonian time; no physical clock calibration",
        "usefulness": "very conservative tiny interval; not an ordinary-physics benchmark",
        "proof_source": "paper/tex_fragments/WHITNEY_QUANTUM_HISTORY.tex",
        "proof_label": "thm:whitney-quantum-trial-history"}, "global error certificate and horizon")
    require(error(horizon) <= target < error(2*horizon), "maximal dyadic horizon")
    require(type(packet["phase_configurations"]) is list and len(packet["phase_configurations"]) == 4,
            "phase configuration census")
    potentials = {}
    names = ("origin", "center-real", "center-imaginary", "center-boundary-relative-phase")
    for row, name in zip(packet["phase_configurations"], names, strict=True):
        require(type(row) is dict and set(row) == {"id", "coordinates_exact", "potential_in_Qsqrt5"}
                and row["id"] == name, "phase configuration schema")
        expected = ["0"]*56
        if name != "origin":
            expected[43 if name == "center-imaginary" else 30] = "1/2"
        if name == "center-boundary-relative-phase":
            expected[44:56] = ["1/4"]*12
        equal(row["coordinates_exact"], expected, "exact phase coordinates")
        potential = sampled_potential([rational(value) for value in expected], elements)
        equal(row["potential_in_Qsqrt5"], potential.pair(), "independent exact sample potential")
        potentials[name] = potential
    require(type(packet["history"]) is list and len(packet["history"]) == 5, "complete trial sample schedule")
    for index, row in enumerate(packet["history"]):
        time = horizon*Q(index, 4)
        equal(row, {"time": str(time), "squared_norm_error_upper": str(error(time)),
            "phase_samples": [{"configuration": name, "angle_in_Qsqrt5": (potentials[name]*-time).pair(),
                "interpretation": "trial_to_initial_wavefunction_ratio=exp(i*angle)"} for name in names]},
            "exact history phase or sample error")
    parent_initial_state()
    return {"accepted": True, "scope": SCOPE, "real_configuration_dimension": 56,
        "state_samples": 5, "phase_configurations": 4, "horizon": str(horizon),
        "target_norm_error": "1/10", "squared_norm_error_upper": str(error(horizon)),
        "kinetic_upper": str(bounds["kinetic_upper"]), "force_moment_upper": str(bounds["force_moment_upper"]),
        "volume_in_Qsqrt5": volume.pair(), "global_time_coverage": True,
        "trial_history_computed": True, "exact_Hamiltonian_history_computed": False,
        "configuration_density_moves": False, "observer_history": False,
        "physical_state_preparation": False, "empirical_comparison": False,
        "ordinary_physics_benchmark": False, "analytic_proof_formalized_in_Lean": False,
        "numeric_quadrature_used_for_bound": False}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), sort_keys=True))
