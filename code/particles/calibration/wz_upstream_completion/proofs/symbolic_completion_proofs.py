#!/usr/bin/env python3
"""Exact symbolic checks for the theorem stack. This is not a self-energy engine."""
from sympy import symbols, Rational, simplify, Matrix

# Anomaly arithmetic
y = {
    "Q": Rational(1, 6), "uc": Rational(-2, 3), "dc": Rational(1, 3),
    "L": Rational(-1, 2), "ec": Rational(1, 1),
}
assert simplify(2 * y["Q"] + y["uc"] + y["dc"]) == 0
assert simplify(3 * y["Q"] + y["L"]) == 0
assert simplify(6 * y["Q"] + 3 * y["uc"] + 3 * y["dc"] + 2 * y["L"] + y["ec"]) == 0
assert simplify(6 * y["Q"]**3 + 3 * y["uc"]**3 + 3 * y["dc"]**3 + 2 * y["L"]**3 + y["ec"]**3) == 0

# Pure-SM one-loop gauge beta coefficients from the declared census.
# Nonabelian: b=-(11/3)C2(G)+(2/3)sum_Weyl T(R)+(1/3)sum_complex T(R).
b3 = -Rational(11, 3) * 3 + Rational(2, 3) * 6
b2 = -Rational(11, 3) * 2 + Rational(2, 3) * 6 + Rational(1, 3) * Rational(1, 2)
# U(1): component multiplicities are explicit.
sumY2_per_gen = (
    6 * Rational(1, 6)**2 + 3 * Rational(2, 3)**2 + 3 * Rational(1, 3)**2
    + 2 * Rational(1, 2)**2 + 1
)
bY = Rational(2, 3) * (3 * sumY2_per_gen) + Rational(1, 3) * (2 * Rational(1, 2)**2)
assert b3 == -7
assert b2 == -Rational(19, 6)
assert bY == Rational(41, 6)

# Neutral tree mass matrix
g, gp, v = symbols("g gp v", nonzero=True)
M = v**2 / 4 * Matrix([[g**2, -g * gp], [-g * gp, gp**2]])
assert simplify(M.det()) == 0
assert simplify(M.trace() - v**2 * (g**2 + gp**2) / 4) == 0

# FJ parameter counterterm identity from the defining equation m2+lambda*vF^2=0.
dm2, dl, m2, lam, dv = symbols("dm2 dl m2 lam dv", nonzero=True)
# The first-order variation is dm2 + dl*v^2 + 2*lambda*v*dv = 0,
# with m2=-lambda*v^2. Solve algebraically for dv/v.
v2 = -m2 / lam
dv_over_v = simplify(-(dm2 + dl * v2) / (2 * lam * v2))
assert simplify(dv_over_v - Rational(1, 2) * (dm2 / m2 - dl / lam)) == 0

# One-loop finite-coordinate orientation.  If pL=pF+h*dp and the two charts
# describe the same observable, then OF1(pF)=OL1(pF)+dp*dO0/dp.
pF, dp, c0, c1, c2, d0, d1, eps_loop = symbols(
    "pF dp c0 c1 c2 d0 d1 eps_loop"
)
O0 = lambda p: c0 + c1 * p + c2 * p**2
OL1 = lambda p: d0 + d1 * p
pL = pF + eps_loop * dp
OF1 = OL1(pF) + dp * (c1 + 2 * c2 * pF)
OF = O0(pF) + eps_loop * OF1
OL = O0(pL) + eps_loop * OL1(pL)
assert simplify((OF - OL).expand().coeff(eps_loop, 0)) == 0
assert simplify((OF - OL).expand().coeff(eps_loop, 1)) == 0

# Charged pole series for Gamma=s-w+h*P1(s)+h^2*P2(s).
h, w = symbols("h w")
P1, P2, P1p, a1, a2 = symbols("P1 P2 P1p a1 a2")
expr = h * a1 + h**2 * a2 + h * P1 + h**2 * (a1 * P1p + P2)
sol1 = -P1
sol2 = -P2 - sol1 * P1p
assert simplify(expr.coeff(h, 1).subs(a1, sol1)) == 0
assert simplify(expr.coeff(h, 2).subs({a1: sol1, a2: sol2})) == 0

# Neutral determinant series at the massive root.
z = symbols("z", nonzero=True)
A1, Z1, Z2, X1, Y1, Z1p, b1, b2c = symbols("A1 Z1 Z2 X1 Y1 Z1p b1 b2c")
coef1 = z * (b1 + Z1)
coef2 = z * (b2c + b1 * Z1p + Z2) + (b1 + A1) * (b1 + Z1) - X1 * Y1
b1sol = -Z1
b2sol = -Z2 - b1sol * Z1p + X1 * Y1 / z
assert simplify(coef1.subs(b1, b1sol)) == 0
assert simplify(coef2.subs({b1: b1sol, b2c: b2sol})) == 0

# Square-root series.
m, s1, s2 = symbols("m s1 s2", nonzero=True)
x = h * s1 / m**2 + h**2 * s2 / m**2
root_series = m * (1 + x / 2 - x**2 / 8)
expected = m + h * s1 / (2 * m) + h**2 * (s2 / (2 * m) - s1**2 / (8 * m**3))
assert simplify(root_series.expand().coeff(h, 0) - expected.expand().coeff(h, 0)) == 0
assert simplify(root_series.expand().coeff(h, 1) - expected.expand().coeff(h, 1)) == 0
assert simplify(root_series.expand().coeff(h, 2) - expected.expand().coeff(h, 2)) == 0

print("PASS anomaly arithmetic")
print("PASS pure-SM one-loop gauge beta coefficients")
print("PASS electroweak neutral mass eigenvalues")
print("PASS FJ parameter-counterterm identity")
print("PASS FJ one-loop finite-coordinate orientation")
print("PASS charged strict pole series through two loops")
print("PASS neutral strict pole series and A-Z mixing order")
print("PASS strict square-root series")
