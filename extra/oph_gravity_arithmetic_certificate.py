from decimal import Decimal, getcontext

getcontext().prec = 120

pi = Decimal(
    "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
)
c = Decimal("299792458")
h = Decimal("6.62607015e-34")
hbar = h / (2 * pi)

B_star = Decimal(
    "3.6078739146803215760518414801601725476877083072171853821242070198789871988886e70"
)
ell2 = (3 * pi) / B_star
ell = ell2.sqrt()
G = (c**3) * ell2 / hbar
nu_cs = Decimal("9192631770")
gamma_cs = ell * nu_cs / c

print("B_star            =", B_star)
print("ell_star^2        =", ell2)
print("ell_star          =", ell)
print("gamma_Cs          =", gamma_cs)
print("hbar              =", hbar)
print("G_OPH             =", G)
print("check B_star*ell2 =", B_star * ell2)
print("check 3*pi        =", 3 * pi)
print("relative residual =", (B_star * ell2 - 3 * pi) / (3 * pi))
