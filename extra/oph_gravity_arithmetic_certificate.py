from decimal import Decimal, getcontext

getcontext().prec = 120

pi = Decimal(
    "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
)
c = Decimal("299792458")
h = Decimal("6.62607015e-34")
hbar = h / (2 * pi)

nu_cs = Decimal("9192631770")
gamma_cs = Decimal(
    "4.9559743365484194636657782433696431879484319705825e-34"
)
epsilon_cs = 2 * pi * gamma_cs
ell = gamma_cs * c / nu_cs
ell2 = ell * ell
B_star = (3 * pi) / ell2
G = (c**5) * gamma_cs * gamma_cs / (hbar * nu_cs * nu_cs)

print("epsilon_Cs        =", epsilon_cs)
print("gamma_Cs          =", gamma_cs)
print("B_star            =", B_star)
print("ell_star^2        =", ell2)
print("ell_star          =", ell)
print("hbar              =", hbar)
print("G_OPH             =", G)
print("check B_star*ell2 =", B_star * ell2)
print("check 3*pi        =", 3 * pi)
print("relative residual =", (B_star * ell2 - 3 * pi) / (3 * pi))
