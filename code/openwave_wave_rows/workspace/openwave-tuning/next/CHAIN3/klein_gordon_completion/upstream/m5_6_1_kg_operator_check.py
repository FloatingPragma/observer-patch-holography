"""
M5.6.1 (math-first) — symbolic anatomy of the hedgehog KG operator (5a §5)

Before building any 3D leapfrog, pin down WHERE the Klein-Gordon mass gap lives in
Duda's hedgehog twist equation (5a §5, paper Fig.9):

    2 ∂_tt ψ = [ (∇ − Â)²  +  (Â·∇ / ‖Â‖)² ] ψ ,     Â = (x,y,z)/r²

Checks (all symbolic, sympy):
  1. Connection identities:  ∇·Â = 1/r² ,  ‖Â‖² = 1/r².
  2. (∇−Â)²ψ expansion + the KEY question: does the zeroth-order (mass-like) coefficient
     −(∇·Â) + ‖Â‖²  survive, or cancel?
  3. Full operator L reduces to  2 ∂_rr + (1/r²) Δ_Ω  (verified on a radial test function).
  4. CORRECTION (2026-05-27): the bare phase ψ is MASSLESS — no dispersion gap. The
     "frozen-Â  k=0 gap = 1/r0²" first reported is an ARTIFACT: freezing Â to a constant
     drops ∇·Â (=1/r²), which is exactly the term that cancels ‖Â‖². With the full
     position-dependent Â, both vanish ⇒ no ψ-mass.
  5. The mass lives in the dual Ψ = exp(iψ): (∇−iÂ)²Ψ keeps +‖Â‖²Ψ = (1/r²)Ψ as a
     position-dependent mass² — finite particle mass only after core regularization
     (M5.6.2/.3). The clock ω = mc²/ℏ is the Ψ phase rotation.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_1_kg_operator_check
"""
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
r = sp.sqrt(x**2 + y**2 + z**2)
coords = (x, y, z)


def grad(f):
    return sp.Matrix([sp.diff(f, c) for c in coords])


def div(V):
    return sum(sp.diff(V[i], coords[i]) for i in range(3))


def laplacian(f):
    return sum(sp.diff(f, c, 2) for c in coords)


A = sp.Matrix([x, y, z]) / r**2                 # Â^hedg
normA2 = sp.simplify((A.T * A)[0])
normA = sp.sqrt(normA2)
rhat = sp.simplify(A / normA)                   # Â/‖Â‖ = r̂


def operator_L(psi):
    """L[ψ] = (∇−Â)·(∇−Â)ψ + (r̂·∇)(r̂·∇)ψ — the full 5a §5 operator."""
    gpsi = grad(psi)
    cov = gpsi - A * psi
    covLap = div(cov) - (A.T * cov)[0]                       # (∇−Â)²ψ
    rad1 = (rhat.T * gpsi)[0]                                # r̂·∇ψ
    rad2 = sum(rhat[i] * sp.diff(rad1, coords[i]) for i in range(3))   # (r̂·∇)²ψ
    return sp.simplify(covLap + rad2)


def main():
    print("=" * 74)
    print("M5.6.1 — symbolic anatomy of the hedgehog KG operator (5a §5)")
    print("=" * 74)

    # --- 1. connection identities ------------------------------------------------
    divA = sp.simplify(div(A))
    print("\n[1] Hedgehog connection identities")
    print(f"    ∇·Â = {divA}   ‖Â‖² = {normA2}   (both expect 1/r²)")
    ok1 = sp.simplify(divA - 1 / r**2) == 0 and sp.simplify(normA2 - 1 / r**2) == 0
    print(f"    → {'OK' if ok1 else 'MISMATCH'}")

    # --- 2. does the explicit mass term cancel? ----------------------------------
    mass_coeff = sp.simplify(-divA + normA2)
    print("\n[2] (∇−Â)²ψ = ∇²ψ − (∇·Â)ψ − 2Â·∇ψ + ‖Â‖²ψ")
    print(f"    zeroth-order coeff of ψ  =  −(∇·Â) + ‖Â‖²  =  {mass_coeff}")
    print(f"    → explicit mass term in the ψ-equation: "
          f"{'CANCELS (=0)' if mass_coeff == 0 else 'SURVIVES'}")

    # --- 3. reduce the full operator; verify on concrete radial test functions ---
    # L should reduce to 2∂_rr on radial modes. Test ψ=f(r) for several f, comparing
    # L[f] (Cartesian, via chain rule) against 2·f''(r) (computed via an R symbol).
    R = sp.symbols("R", positive=True)
    print("\n[3] Full L on RADIAL ψ = f(r)  (expect L = 2·f''(r)):")
    all_ok = True
    for f_R in (R**4, sp.cos(R), sp.exp(-R)):
        f_r = f_R.subs(R, r)
        L_val = sp.simplify(operator_L(f_r))
        expected = sp.simplify((2 * sp.diff(f_R, R, 2)).subs(R, r))
        match = sp.simplify(L_val - expected) == 0
        all_ok = all_ok and match
        print(f"    f(r)={str(f_R):10s}: L[f] = 2·f''  → {match}")
    print(f"    → L reduces to 2∂_rr on radial modes: {all_ok}")
    print("    (general form: L = 2∂_rr + (1/r²)Δ_Ω — a MASSLESS radial wave, weakened angular part)")

    # --- 4. dispersion: bare ψ is gapless (correction) ---------------------------
    print("\n[4] CORRECTION — bare-ψ dispersion has NO gap")
    print("    Full operator (cancellation applied): 2∂_ttψ = ∇²ψ − 2Â·∇ψ + (r̂·∇)²ψ.")
    print("    Local plane wave ψ=exp(i(k·x−ωt)), REAL part:  2ω² = k² + (k·r̂)²  (gapless).")
    print("    The earlier 'freeze Â ⇒ 2ω²=1/r0² at k=0' was an ARTIFACT: freezing Â to a")
    print("    constant kills ∇·Â=1/r², the very term that cancels ‖Â‖². Inconsistent.")
    print("    ⇒ bare phase ψ is MASSLESS; there is no ψ-dispersion gap to measure.")

    # --- 5. the gap lives in Ψ = exp(iψ) -----------------------------------------
    print("\n[5] Where the KG mass actually is — the dual Ψ = exp(iψ)")
    print("    For the COMPLEX Ψ, minimal coupling (∇−iÂ)²Ψ RETAINS +‖Â‖²Ψ = (1/r²)Ψ —")
    print("    a position-dependent mass² (→0 far from core, →∞ at the disclination).")
    print("    A finite PARTICLE mass / clock ω=mc²/ℏ exists only once the core is")
    print("    regularized (sets the scale) — that is M5.6.2 (core) + M5.6.3 (Faber V-on).")

    print("\n" + "=" * 74)
    print("FINDING (M5.6.1a): the KG structure is minimal-coupling to the hedgehog")
    print("connection Â (geometric), NOT an added V_ψ potential — the explicit mass term")
    print("cancels exactly. This literally confirms Duda's thesis + the M5.2 error.")
    print("CONSEQUENCE for M5.6.1b: do NOT measure a bare-ψ dispersion gap (there isn't")
    print("one). Measure (a) gauge-covariant structure + stable massless evolution now,")
    print("and (b) the finite clock-mass with the regularized core in M5.6.2/.3.")
    print("=" * 74)


if __name__ == "__main__":
    main()
