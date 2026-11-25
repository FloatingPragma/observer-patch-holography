"""Cyclotomic coefficient witnesses for the three balanced Fourier intensities."""

from argparse import ArgumentParser
from itertools import combinations
from pathlib import Path
import json

from sympy import Poly, cyclotomic_poly, symbols

T, Z = symbols("T Z")
PHI = Poly(cyclotomic_poly(24, T), T, domain="ZZ")


def reduce_coefficient(value):
    """T is a primitive 24th root; coefficients here are exact integers."""
    return Poly(value, T, domain="ZZ").rem(PHI).as_expr()


def channel(k):
    # In the standard complex embedding, omega=T^8 and sqrt(2)=T^3+T^(-3).
    q = T ** ((8 * k) % 24) * Z**2 + (T**3 + T**21) * Z + T ** ((-8 * k) % 24)
    squared = Poly(q**2, Z)
    return [reduce_coefficient(squared.nth(m)) for m in range(5)]


def minor(a, b):
    """Determinant of the constant and fourth-power coefficients."""
    return reduce_coefficient(a[0] * b[4] - b[0] * a[4])


def certificate():
    assert PHI.as_expr() == T**8 - T**4 + 1
    assert reduce_coefficient((T**3 + T**21) ** 2) == 2
    assert reduce_coefficient(1 + T**8 + T**16) == 0
    channels = [channel(k) for k in range(3)]
    assert all(p[2] == 4 and p[4] != 0 for p in channels)
    omega_difference = reduce_coefficient(T**8 - T**16)
    pairs = {}
    for i, j in combinations(range(3), 2):
        det = minor(channels[i], channels[j])
        assert det in (omega_difference, -omega_difference) and det != 0
        pairs[f"{i},{j}"] = str(det)
    top_light = [x - y for x, y in zip(channels[0], channels[1], strict=True)]
    mid_light = [x - y for x, y in zip(channels[2], channels[1], strict=True)]
    gap = minor(top_light, mid_light)
    assert gap in (3 * omega_difference, -3 * omega_difference) and gap != 0
    return dict(
        primitive_root_order=24,
        minimal_polynomial=str(PHI.as_expr()),
        variable="Z=exp(2i/9)",
        channel_coefficients_ascending=[[str(c) for c in p] for p in channels],
        pair_minors=pairs,
        common_offset_invariant_gap_minor=str(gap),
    )


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    out = json.dumps(certificate(), indent=2) + "\n"
    if parser.parse_args().write:
        Path(__file__).with_name("results.json").write_text(out)
    else:
        print(out, end="")


if __name__ == "__main__":
    main()
