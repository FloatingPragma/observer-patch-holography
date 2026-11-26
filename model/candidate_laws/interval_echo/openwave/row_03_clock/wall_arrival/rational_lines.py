"""Exact reflection polynomials and Rabin irreducibility certificates."""

from sympy import Poly, factorint, symbols

X = symbols("x")
# A witness prime for each enrolled reflection block, found after acquisition.
PRIMES = {14: 3, 30: 11, 50: 3}


def reflection_polynomial(m):
    if m % 2:
        raise ValueError("These certificates use a weak central bond")
    p, q = Poly(1, X), Poly(X, X)
    for k in range(2, m + 1):
        off = 16 if k % 2 == 0 else 1
        diagonal = 1 if k == m else 0
        p, q = q, Poly(X - diagonal, X) * q - off**2 * p
    return q


def certificate(m):
    poly = reflection_polynomial(m)
    assert poly.degree() == m and poly.all_coeffs()[:2] == [1, -1]
    prime = PRIMES[m]
    f = poly.set_modulus(prime)
    x = Poly(X, X, modulus=prime)
    required = {m // int(q) for q in factorint(m)}
    xp = x
    checks = []
    for i in range(1, m + 1):
        xp = (xp**prime).rem(f)
        if i in required:
            checks.append((i, (xp - x).gcd(f).degree()))
    return prime, checks, xp == x


if __name__ == "__main__":
    for n in PRIMES:
        print(n, certificate(n))
