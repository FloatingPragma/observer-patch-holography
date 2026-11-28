"""Exact spin-half spectral continuation under one supplied ring twist."""

import argparse
import json
from math import comb
from pathlib import Path

import sympy as s

from .. import certificate, ring


X, Z, U = s.symbols("x z u")


def forward(up, down):
    words, index = ring.sector(4, up, down)
    hop = s.zeros(len(words))
    onsite = s.zeros(len(words))
    translate = s.zeros(len(words))
    phases = []
    for column, word in enumerate(words):
        loads = [((word >> (2 * x)) & 3).bit_count() for x in range(4)]
        onsite[column, column] = 16 * loads.count(2)
        phases.append(s.I ** (sum(x * n for x, n in enumerate(loads)) % 4))
        target, phase = ring.shifted(word, 4)
        translate[index[target], column] = phase
        for old in range(8):
            new = (old + 2) % 8
            if word & (1 << old) and not word & (1 << new):
                missing, first = ring.car(word, old)
                target, second = ring.car(missing, new)
                hop[index[target], column] = -first * second
    gauge = s.diag(*phases)
    assert onsite + hop + hop.T == certificate.matrix(ring.hamiltonian(4, up, down))
    assert gauge * hop * gauge.H == s.I * hop
    assert translate * gauge.H == s.I ** ((up + down) % 4) * gauge.H * translate
    assert all(translate * h == h * translate for h in (onsite, hop))
    return onsite, hop


def frame(up, down, k):
    orbits = ring.orbits(4, up, down, k)
    basis = s.zeros(len(ring.sector(4, up, down)[0]), len(orbits))
    for j, orbit in enumerate(orbits):
        assert len(orbit) == 4
        for step, (i, phase) in enumerate(orbit):
            basis[i, j] = phase * (-s.I) ** (k * step % 4) / 2
    assert basis.H * basis == s.eye(len(orbits))
    return basis


def polynomial(up, down, k):
    onsite, hop = forward(up, down)
    top_onsite, top_hop = forward(up + 1, down - 1)
    raising = certificate.matrix(ring.spin_plus(4, up, down))
    assert all(
        raising * lower == upper * raising
        for lower, upper in (
            (onsite, top_onsite),
            (hop, top_hop),
            (hop.T, top_hop.T),
        )
    )
    b = frame(up, down, k)
    top_b = frame(up + 1, down - 1, k)
    assert b.cols == 6 and top_b.cols == 1 and (raising * b).rank() == 1

    def make(d, f):
        return d + Z * f + f.T / Z

    full = (b.H * make(onsite, hop) * b).applyfunc(s.expand)
    spin_three = s.expand((top_b.H * make(top_onsite, top_hop) * top_b)[0])
    factor = s.cancel((X * s.eye(6) - full).det(method="domain-ge") / (X - spin_three))
    half_angle = (1 + s.I * U) / (1 - s.I * U)
    numerator = s.Poly(s.cancel(factor.subs(Z, half_angle) * (1 + U**2) ** 5), X, U)
    assert numerator.domain == s.ZZ
    return numerator.as_expr(), full


def discriminant_factor(numerator):
    determinant = s.discriminant(numerator, X)
    known = 2**44 * (1 - U**2) ** 4 * (1 + U**2) ** 24
    residual = s.Poly(s.cancel(determinant / known), U, domain=s.ZZ)
    assert residual.degree() == 24
    powers = [residual.nth(j) for j in range(25)]
    bernstein = [
        sum(powers[j] * s.Rational(comb(k, j), comb(24, j)) for j in range(k + 1))
        for k in range(25)
    ]
    assert min(bernstein) > 0
    return int(min(bernstein))


def endpoint_slope(full, sign, k):
    # Interpolation over the exact distinct energies of the six-cell block.
    shift = 8 * (1 - sign)
    middle, radical, isolated = (
        (24, 2 * s.sqrt(21), (18, 30)) if k == 3 else (24, 2 * s.sqrt(13), (14, 34))
    )
    low = middle - radical - shift
    rest = (middle + radical - shift, *(v - shift for v in isolated))
    stop = full.subs(Z, s.I)
    identity = s.eye(6)
    actual_charpoly = stop.charpoly(X).as_expr()
    assert (
        s.simplify(
            actual_charpoly
            - (X - low) ** 2 * (X - rest[0]) ** 2 * (X - rest[1]) * (X - rest[2])
        )
        == 0
    )
    projector = identity
    for eig in rest:
        projector = (projector * (stop - eig * identity) / (low - eig)).applyfunc(
            s.simplify
        )
    assert s.trace(projector) == 2
    derivative = (s.I * Z * full.diff(Z)).subs(Z, s.I)
    compressed = (projector * derivative * projector).applyfunc(s.simplify)
    assert s.simplify(s.trace(compressed)) == 0
    slope_squared = s.simplify(s.trace(compressed**2) / 2)
    assert slope_squared > 0
    return low / 16, slope_squared


def static_residues(sign):
    up, down = (3, 2) if sign > 0 else (2, 1)
    neel = sum(1 << (2 * x + x % 2) for x in range(4))
    seed = s.zeros(36, 1)
    seed[ring.sector(4, 2, 2)[1][neel]] = 1
    ground = certificate.projected(2, 2, 0, seed)
    crank = certificate.matrix(ring.field(4, 2, 2, 0, int(sign < 0), sign > 0))
    input_vector = crank * ground
    residues = []
    for index in (0, 2):
        projected = certificate.projected(up, down, index, input_vector)
        fraction = projected.dot(projected) / (2 * ground.dot(ground))
        error = s.Rational(1, 10**15)
        residues.append(certificate.interval(fraction - error, fraction + error, 11))
    return residues


def run():
    ring.source()
    results = []
    for sign in (1, -1):
        up, down = (3, 2) if sign > 0 else (2, 1)
        levels = certificate.sector(up, down)[1]
        residues = static_residues(sign)
        for k in (3, 1):
            numerator, full = polynomial(up, down, k)
            at_start = s.Poly(numerator.subs(U, 0), X)
            assert at_start.rem(levels[0].factor).is_zero
            assert at_start.count_roots(-s.oo, levels[0].hi) == 1
            bound = discriminant_factor(numerator)
            high, slope = endpoint_slope(full, sign, k)
            if k == 3:
                assert s.simplify(slope - (40 + 8 * s.sqrt(21)) / 21) == 0
                assert levels[2].lo < 16 * high < levels[2].hi
            else:
                assert 16 * high > levels[2].hi
            final_half = s.Poly(numerator.subs(U, 1), X)
            left = s.floor(16 * high * 10**8) / 10**8
            right = left + s.Rational(1, 10**8)
            assert final_half.count_roots(-s.oo, left) == 0
            assert final_half.count_roots(left, right) == 1
            results.append(
                {
                    "charge": sign,
                    "initial_translation_character": k,
                    "after_large_gauge_return_character": (k + sign) % 4,
                    "minimum_Bernstein_coefficient_of_discriminant_residual": bound,
                    "ending_unshifted_H_energy": str(high),
                    "initial_unshifted_H_energy_bounds": certificate.interval(
                        levels[0].lo / 16, levels[0].hi / 16, 10
                    ),
                    "terminal_multiplicity_at_fixed_spin_projection": 2,
                    "squared_16H_slope_in_bond_phase": str(slope),
                    "ends_in_second_bright_line": k == 3,
                    "lower_bright_start_local_CAR_residue": residues[0],
                    "ending_local_CAR_residue_of_one_slope_direction": residues[1]
                    if k == 3
                    else None,
                }
            )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2) + "\n"
    receipt = Path(__file__).with_name("results.json")
    if args.write:
        receipt.write_text(output, encoding="utf-8")
    else:
        assert receipt.read_text(encoding="utf-8") == output
        print(
            "Verified four spin-half contact-ring flux endpoint maps and CAR residues."
        )


if __name__ == "__main__":
    main()
