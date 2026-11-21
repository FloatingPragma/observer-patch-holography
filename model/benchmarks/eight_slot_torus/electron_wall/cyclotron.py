"""Two frequency transitions of the wall walk on a uniform magnetic torus."""
import argparse
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from .__main__ import HERE
from .walk import dimer, rest_properties

PROTOCOL = json.loads((HERE / "cyclotron.json").read_text())
WEAK_PERIODS = json.loads((HERE / "weak_field.json").read_text())["periods"]


def enrolled(rows, periods):
    expected = [(n, q) for n in periods for q in PROTOCOL["charges"]]
    return bool(expected) and [(r["period"], r["charge"]) for r in rows] == expected


def magnetic_block(n, length, charge, grade):
    """12 times the discriminant in a mirror-spin grade, with a real spin basis."""
    t = sparse.coo_matrix((np.ones(n), ((np.arange(n) + 1) % n, np.arange(n))), shape=(n, n)).tocsr()
    unit = sparse.eye(n)
    angle = 2 * np.pi * charge * np.arange(n) / n
    # a = 1 + (1-cos kx) + (1-cos(ky-q A)), with ky=kz=0.
    a = sparse.diags(3 - np.cos(angle)) - (t + t.T) / 2
    wall_reflection = np.fliplr(np.eye(length))
    sigma_z = np.diag([1., -1.])
    skew = np.array([[0., 1.], [-1., 0.]])
    sigma_x = np.array([[0., 1.], [1., 0.]])
    def kron(left, right):
        return sparse.kron(left, right, format="csr")
    return (-kron((t - t.T) / 2, np.kron(np.eye(length), skew))
            - kron(sparse.diags(np.sin(angle)), np.kron(np.eye(length), sigma_x))
            + grade * kron(a, np.kron(wall_reflection, sigma_z))
            + grade * kron(unit, np.kron(dimer(length, 0) @ wall_reflection, sigma_z)))


def oscillator_sources(n, length, grade):
    b = 2 * np.pi / n
    centered = (np.arange(n) + n // 2) % n - n // 2
    gaussian = np.exp(-b * centered ** 2 / 2)
    wall = np.linalg.svd(dimer(length, 1))[0][:, -1]
    for envelope in (gaussian, centered * gaussian):
        probe = np.zeros((n, length, 2))
        probe[..., 0 if grade == 1 else 1] = envelope[:, None] * wall
        yield probe.ravel() / np.linalg.norm(probe)


def modes(n, length, charge, grade, count):
    matrix = magnetic_block(n, length, charge, grade)
    mass = rest_properties(length)["singular_mass"]
    seed = np.random.default_rng(171).standard_normal(matrix.shape[0])
    roots, vectors = eigsh(matrix, k=count, sigma=mass, v0=seed, tol=1e-12)
    index = np.argsort(roots)
    roots, vectors = roots[index], vectors[:, index]
    residual = np.linalg.norm(matrix @ vectors - vectors * roots, axis=0).max()
    if not np.isfinite(residual) or residual > 1e-9:
        raise ArithmeticError("unresolved eigenvectors")
    weights, winners = [], []
    for probe in oscillator_sources(n, length, grade):
        raw = vectors.T @ probe
        prob = raw ** 2
        weights.append([*prob.tolist(), float(max(0., 1 - prob.sum()))])
        winners.append(int(np.argmax(prob)))
    if winners[0] == winners[1] or min(ws[i] for ws, i in zip(weights, winners)) <= .5:
        raise ArithmeticError("orbital detector could not distinguish two dominant lines")
    rates = np.arcsin(roots / 12)
    return dict(rates=rates.tolist(), weights=weights, assigned=winners,
                residual=float(residual))


def experiment(spec):
    results = []
    for n in spec["periods"]:
        for q in spec["charges"]:
            grades = [modes(n, spec["private_length"], q, g, spec["levels_per_grade"]) for g in (1, -1)]
            selected = [np.array(g["rates"])[g["assigned"]] for g in grades]
            orbital = np.array([rates[1] - rates[0] for rates in selected])
            spin_gap = q * (selected[1][0] - selected[0][0])
            if np.any(orbital <= 0) or spin_gap <= 0:
                raise ArithmeticError("wrong order of the observed lines")
            results.append(dict(period=n, charge=q, grades=grades, spin_gap=float(spin_gap),
                                orbital=orbital.tolist(), g_each=(2 * spin_gap / orbital).tolist(),
                                g=float(2 * spin_gap / np.mean(orbital))))
    return results


def check_reversal(results):
    for i in range(0, len(results), 2):
        first, second = results[i:i + 2]
        if abs(first["g"] - second["g"]) > 1e-8:
            raise ArithmeticError("opposite charges give different magnetic factors")
        for a, b in zip(first["grades"], reversed(second["grades"])):
            if np.max(np.abs(np.array(a["rates"]) - b["rates"])) > 1e-9:
                raise ArithmeticError("charge reversal did not exchange spectra")


def initial_gate(results, target):
    if not enrolled(results, PROTOCOL["periods"]):
        return False
    for charge in PROTOCOL["charges"]:
        errors = [abs(r["g"] - target) for r in results if r["charge"] == charge]
        if not (errors[-1] < .05 and all(y < x for x, y in zip(errors, errors[1:]))):
            return False
    return True


def weak_gate(results, initial, target):
    if not enrolled(initial, PROTOCOL["periods"]) or not enrolled(results, WEAK_PERIODS):
        return False
    for q in (-1, 1):
        previous = next(r for r in reversed(initial) if r["charge"] == q)
        for i, row in enumerate(r for r in results if r["charge"] == q):
            before = np.array([previous["g"], *previous["g_each"]])
            after = np.array([row["g"], *row["g_each"]])
            extrapolated = 2 * row["g"] - previous["g"]
            if row["period"] != 2 * previous["period"] or not np.all(abs(after-target) < abs(before-target)):
                return False
            if i > 0 and (extrapolated >= 2 or abs(extrapolated - target) >= .01):
                return False
            previous = row
    return True


def compare(stored, results):
    if not isinstance(stored, list) or len(stored) != len(results):
        raise ValueError("incorrect acquisition size")
    def match(left, right):
        array = np.asarray(left)
        if array.shape != np.shape(right) or array.dtype.kind not in "fiu" or not np.isfinite(array).all():
            raise ValueError("invalid numerical record")
        if np.asarray(right).dtype.kind in "iu":
            if array.dtype.kind not in "iu" or not np.array_equal(array, right):
                raise ValueError("invalid integer label in magnetic receipt")
            return
        np.testing.assert_allclose(array, right, atol=1e-9, rtol=1e-9)
    for entry, fresh in zip(stored, results, strict=True):
        if entry.keys() != fresh.keys() or len(entry["grades"]) != 2:
            raise ValueError("invalid line record")
        for name in ("period", "charge", "spin_gap", "orbital", "g_each", "g"):
            match(entry[name], fresh[name])
        for old, current in zip(entry["grades"], fresh["grades"], strict=True):
            if old.keys() != current.keys():
                raise ValueError("invalid grade")
            for name in current:
                match(old[name], current[name])


def magnetic_series(weak_field, *, acquire=False):
    periods = WEAK_PERIODS if weak_field else PROTOCOL["periods"]
    results = experiment(dict(PROTOCOL, periods=periods))
    if not acquire:
        filename = "weak_lines.json" if weak_field else "cyclotron_lines.json"
        compare(json.loads((HERE / filename).read_text()), results)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=str, help="write every line, both probes and their complement")
    parser.add_argument("--weak-field", action="store_true", help="run the independent halved-field series")
    parser.add_argument("--require-gate", action="store_true", help="return failure if the scientific gate fails")
    args = parser.parse_args()
    results = magnetic_series(args.weak_field, acquire=args.acquire is not None)
    if args.acquire:
        with open(args.acquire, "x", encoding="utf-8") as output:
            json.dump(results, output, separators=(",", ":"), allow_nan=False)
            output.write("\n")
    for row in results:
        print(row["period"], row["charge"], row["g_each"], row["g"])
    check_reversal(results)
    target = rest_properties(PROTOCOL["private_length"])["g"]
    if args.weak_field:
        initial = magnetic_series(False)
        accepted = weak_gate(results, initial, target)
    else:
        accepted = initial_gate(results, target)
    print("Scientific gate:", "PASS" if accepted else "FAIL")
    if args.require_gate and not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
