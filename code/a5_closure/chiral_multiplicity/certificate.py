"""Exact family-index and weak-strength checks on the declared exterior menu."""

from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import sys

from sympy import I, Matrix, diag, eye, kronecker_product, sqrt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import super_tannakian_matter_lift_certificate as matter  # noqa: E402
from multiplicity_window_certificate import ckm_parameter_count  # noqa: E402


COMPONENTS = matter._scan_components(-2, 3)
# Selected even matter indices in the order Q, u^c, d^c, L, e^c.
DIRECT = (3, 2, 8, 9, 4)
MIRROR = (6, 5, 0, 1, 7)


def columns(indices):
    result = []
    for i in indices:
        values = matter._scan_anomalies([COMPONENTS[i]])
        result.append([values[k] for k in ("grav", "su3", "su2", "u1_cubed")])
    return Matrix(result).T


def inspect(counts):
    material = [c for c, n in zip(COMPONENTS, counts, strict=True) for _ in range(n)]
    anomalies = matter._scan_anomalies(material)
    cube = sum(c["weak"] * {"1": 0, "3": 1, "3bar": -1}[c["color"]] for c in material)
    q2, raising, dimension = Fraction(0), 0, 0
    # T+ maps the lower to the upper member of one normalized doublet.
    for c in material:
        colors = 1 if c["color"] == "1" else 3
        t3 = (-3, 3) if c["weak"] == 2 else (0,)
        q2 += colors * sum(Fraction(c["q"] + t, 6) ** 2 for t in t3)
        raising += colors * (c["weak"] == 2)
        dimension += colors * c["weak"]
    # Presence of q=±6 uncoloured weak singlet restricts a central kernel
    # element to a sixth root in U(1). This applies to every test material.
    assert any(
        c["color"] == "1" and c["weak"] == 1 and abs(c["q"]) == 6 for c in material
    )
    center = [
        (a, b, k)
        for a, b, k in product(range(3), range(2), range(6))
        if all(
            (
                2 * a * {"1": 0, "3": 1, "3bar": -1}[c["color"]]
                + 3 * b * (c["weak"] - 1)
                + k * c["q"]
            )
            % 6
            == 0
            for c in material
        )
    ]
    assert not any(anomalies.values()) and cube == 0
    assert set(center) == {(k % 3, k % 2, k) for k in range(6)}
    net = [counts[i] - counts[j] for i, j in zip(DIRECT, MIRROR, strict=True)]
    assert len(set(net)) == 1
    return dict(
        index=net[0],
        dimensions=dimension,
        weak_raising_strength=raising,
        charge_square_trace=str(q2),
        bare_weak_doublets=raising,
        witten_parity=raising % 2,
        mixed_and_cubic_anomalies=anomalies,
        pure_color_cube=cube,
        faithful_quotient_kernel_order=len(center),
    )


def generations(left, right=0):
    m = [0] * 10
    for i, j in zip(DIRECT, MIRROR, strict=True):
        m[i], m[j] = left, right
    return m


def mass_and_spin(n):
    # One-Higgs positive diagonal Yukawa comparator, in units v/sqrt(2)=1.
    # The spinor factor is supplied independently; these are colourless Qem=-1 states.
    sx = kronecker_product(eye(n), -I * Matrix([[0, 1], [1, 0]]))
    sy = kronecker_product(eye(n), Matrix([[0, -1], [1, 0]]))
    assert sx**2 == sy**2 == -eye(2 * n) and sx * sy == -sy * sx
    for masses in [list(range(1, n + 1)), [1] * n]:
        mass = kronecker_product(diag(*masses), eye(2))
        assert mass * sx == sx * mass and mass * sy == sy * mass
    return dict(
        charge=-1,
        rotation_order_sign=-1,
        allowed_sector_mass_diagonals=[list(range(1, n + 1)), [1] * n],
    )


def run():
    a = columns(DIRECT)
    assert columns(MIRROR) == -a
    assert a.rank() == 4 and a.nullspace() == [Matrix.ones(5, 1)]
    # Every deletion of one component from a generation violates the finite anomaly gate.
    deletion_signatures = {}
    for name, i in zip(("Q", "u_c", "d_c", "L", "e_c"), DIRECT, strict=True):
        missing = generations(1)
        missing[i] = 0
        deletion_signatures[name] = matter._scan_anomalies(
            [c for c, n in zip(COMPONENTS, missing, strict=True) if n]
        )
        assert any(deletion_signatures[name].values())
    band = json.loads(
        (HERE.parent / "manifests/family_band_attachment_reference.json").read_text()
    )
    dims = band["spectral_resolution"]["band_dimensions"]
    admissible = {
        dims[name] for name, size in band["band_action_kernels"].items() if size == 1
    }
    selected = dims[band["selection"]["minimizer"]]
    assert band["selection"]["strict"] and admissible == {3, 5} and selected == 3
    replicas = {}
    for n in range(1, 7):
        record = inspect(generations(n))
        assert record["weak_raising_strength"] == 4 * n and Fraction(
            record["charge_square_trace"]
        ) == Fraction(16 * n, 3)
        replicas[n] = dict(
            **record,
            **mass_and_spin(n),
            intrinsic_ckm_phases=ckm_parameter_count(n)["phases"],
            one_loop_weak_coefficient=str(Fraction(43 - 8 * n, 6)),
            admissible_complete_faithful_port_band=n in admissible,
            minimum_cost_family_band=n == selected,
        )
    # An explicit CP-capable three-generation block embeds in n=3,4,5.
    omega = (-1 + I * sqrt(3)) / 2
    fourier = Matrix(
        [[omega ** ((i * j) % 3) / sqrt(3) for j in range(3)] for i in range(3)]
    )
    assert (fourier * fourier.H - eye(3)).applyfunc(
        lambda z: z.simplify()
    ) == Matrix.zeros(3)
    invariant = (
        (
            fourier[0, 0]
            * fourier[1, 1]
            * fourier[0, 1].conjugate()
            * fourier[1, 0].conjugate()
        )
        .as_real_imag()[1]
        .simplify()
    )
    assert invariant**2 == Fraction(1, 108)
    mixed = inspect(generations(2, 1))
    for k in ("dimensions", "charge_square_trace", "weak_raising_strength"):
        assert mixed[k] == replicas[3][k]
    extra = generations(3)
    extra[DIRECT[-1]] += 1
    extra[MIRROR[-1]] += 1
    return dict(
        field_order=["Q", "u_c", "d_c", "L", "e_c"],
        direct=list(DIRECT),
        mirror=list(MIRROR),
        four_anomaly_rows=a.tolist(),
        kernel_rank=1,
        one_missing_component=deletion_signatures,
        replicas=replicas,
        cp_block_j_squared=str(invariant**2),
        two_generations_plus_mirror=mixed,
        three_plus_vectorlike_charged_singlet=inspect(extra),
    )


if __name__ == "__main__":
    text = json.dumps(run(), default=lambda i: int(i), indent=2) + "\n"
    if sys.argv[1:] == ["--write"]:
        (HERE / "results.json").write_text(text)
    else:
        print(text, end="")
