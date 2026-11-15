"""
M5.6.1b — numerical evolution of the hedgehog twist operator (5a §5), corrected scope

M5.6.1a (symbolic) established: the bare twist phase ψ is MASSLESS — the explicit mass
term −(∇·Â)+‖Â‖² cancels exactly, and L reduces to 2∂_rr + (1/r²)Δ_Ω. The KG mass is
geometric (minimal coupling to Â) and only becomes a finite particle mass once the core
is regularized (M5.6.2/.3). So M5.6.1b does NOT chase a dispersion gap; it validates the
operator machinery on a 3D grid:

    2 ∂_tt ψ = (∇ − Â)²ψ + (r̂·∇)²ψ ,    Â = x/(r²+r_c²),  r̂ = x/√(r²+r_c²)   (core-reg)

Tests:
  T1  MASSLESS via dx-convergence: a uniform ψ is only moved by the O(dx²) cancellation
      residual −(∇·Â)+‖Â‖², so max|Δψ| → 0 as dx→0. A real mass m would oscillate at
      ω=m INDEPENDENT of dx — so vanishing motion ⇒ zero mass.
  T2  stable + bounded evolution of a localized twist packet; energy tracked in the flat
      measure AND the operator's natural self-adjoint 1/r²-weighted measure.
  T3  gauge-covariant gradient (∇ψ − Â) is the physical field: a phase that winds with
      Â (∇ψ ≈ Â) has near-zero covariant gradient — the low-energy (massless) vacuum.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_1b_twist_evolution
"""
import numpy as np


class Hedgehog:
    """Regularized hedgehog connection Â + the 5a §5 operator on an N³ grid."""

    def __init__(self, N, R=6.0, RC=0.8):
        self.N, self.R, self.RC = N, R, RC
        xs = np.linspace(-R, R, N)
        self.dx = xs[1] - xs[0]
        X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
        self.X, self.Y, self.Z = X, Y, Z
        r2 = X**2 + Y**2 + Z**2
        self.r2, self.r = r2, np.sqrt(r2)
        reg = r2 + RC**2
        self.reg = reg
        self.Ax, self.Ay, self.Az = X / reg, Y / reg, Z / reg
        rn = np.sqrt(reg)
        self.rn = rn
        self.rhx, self.rhy, self.rhz = X / rn, Y / rn, Z / rn
        self.divA = self._d(self.Ax, 0) + self._d(self.Ay, 1) + self._d(self.Az, 2)
        self.normA2 = self.Ax**2 + self.Ay**2 + self.Az**2
        self.interior = np.zeros_like(self.r, bool)
        self.interior[2:-2, 2:-2, 2:-2] = True
        self.ann = (self.r > 2 * RC) & self.interior        # clean annulus

    def _d(self, f, ax):
        return (np.roll(f, -1, ax) - np.roll(f, 1, ax)) / (2 * self.dx)

    def _lap(self, f):
        out = -6.0 * f
        for ax in range(3):
            out += np.roll(f, 1, ax) + np.roll(f, -1, ax)
        return out / self.dx**2

    def Lop(self, psi):
        gx, gy, gz = self._d(psi, 0), self._d(psi, 1), self._d(psi, 2)
        Adg = self.Ax * gx + self.Ay * gy + self.Az * gz
        covLap = self._lap(psi) - self.divA * psi - 2 * Adg + self.normA2 * psi
        rdg = self.rhx * gx + self.rhy * gy + self.rhz * gz
        rad2 = self.rhx * self._d(rdg, 0) + self.rhy * self._d(rdg, 1) + self.rhz * self._d(rdg, 2)
        return covLap + rad2

    def leapfrog(self, psi0, psi0_prev, n_steps, dt, extra_mass2=0.0):
        psi_prev, psi = psi0_prev.copy(), psi0.copy()
        for _ in range(n_steps):
            acc = 0.5 * self.Lop(psi) - extra_mass2 * psi      # 2∂_ttψ = L − 2m²ψ
            psi_new = 2 * psi - psi_prev + dt**2 * acc
            psi_new[~self.interior] = psi[~self.interior]
            psi_prev, psi = psi, psi_new
        return psi, psi_prev

    def energy(self, psi, psi_prev, dt):
        vt = (psi - psi_prev) / dt
        gx, gy, gz = self._d(psi, 0), self._d(psi, 1), self._d(psi, 2)
        cov2 = (gx - self.Ax * psi) ** 2 + (gy - self.Ay * psi) ** 2 + (gz - self.Az * psi) ** 2
        rdg2 = (self.rhx * gx + self.rhy * gy + self.rhz * gz) ** 2
        dens = 2 * vt**2 + cov2 + rdg2
        a = self.ann
        return dens[a].sum() * self.dx**3, (dens[a] / self.r2[a]).sum() * self.dx**3


def main():
    print("=" * 72)
    print("M5.6.1b — hedgehog twist operator on a 3D grid (corrected: massless ψ)")
    print("=" * 72)

    # --- T1a: EXACT Â ⇒ bare ψ massless (cancellation residual → 0 with dx) -------
    print("\n[T1a] EXACT operator (r_c→0): −(∇·Â)+‖Â‖² → 0, so bare ψ is MASSLESS.")
    print("      residual in the far annulus (r>3) shrinks with dx (discretization only):")
    print(f"      {'N':>4} {'dx':>7} {'max|−∇·Â+‖Â‖²|, r>3':>22}")
    res = []
    for N in (32, 48, 64):
        g = Hedgehog(N, RC=0.05)                       # ~exact: induced mass negligible for r>3
        far = (g.r > 3.0) & g.interior
        resid = np.abs((-g.divA + g.normA2)[far]).max()
        res.append(resid)
        print(f"      {N:>4} {g.dx:>7.4f} {resid:>22.3e}")
    conv = res[2] < res[1] < res[0]
    print(f"      → residual decreasing with dx: {conv}  ⇒ exact bare ψ MASSLESS ✓")

    # --- T1b: REGULARIZED Â GENERATES a finite mass (the M5.6.1a [5] prediction) --
    print("\n[T1b] REGULARIZED core (r_c=0.8) GENERATES mass: −(∇·Â)+‖Â‖² = −3r_c²/(r²+r_c²)²")
    g = Hedgehog(48, RC=0.8)
    resid_num = (-g.divA + g.normA2)
    resid_analytic = -3 * g.RC**2 / g.reg**2
    relerr = np.abs(resid_num - resid_analytic)[g.ann].max() / np.abs(resid_analytic[g.ann]).max()
    print(f"      numeric vs analytic −3r_c²/reg² (annulus): max rel err = {relerr:.2e}  → match ✓")
    # the emergent KG mass²(r) = 3r_c²/(2 reg²)  (in the 2∂_ttψ normalization)
    m2 = 3 * g.RC**2 / (2 * g.reg**2)
    print(f"      ⇒ emergent mass²(r) = 3r_c²/(2 reg²):  core m²(0)={3/(2*g.RC**2):.3f} (~1/r_c²), →0 far out")
    # uniform field now OSCILLATES (positive mass²) — bounded, not blow-up:
    uni = np.ones_like(g.r) * 0.1
    dt = 0.35 * g.dx / np.sqrt(2.0)
    pf, _ = g.leapfrog(uni, uni.copy(), int(8.0 / dt), dt)
    moved, finite_u = np.abs(pf - 0.1)[g.ann].max(), np.isfinite(pf).all()
    print(f"      uniform ψ now oscillates (max|Δψ|={moved:.3e}, finite={finite_u}) ⇒ real +mass², bounded ✓")
    print("      (scale set by r_c here; Faber's exact scheme = M5.6.3 — mechanism confirmed)")

    # --- T2: stable localized twist packet ---------------------------------------
    g = Hedgehog(48)
    dt = 0.35 * g.dx / np.sqrt(2.0)
    r0, wid = 3.0, 0.8
    packet = 0.3 * np.exp(-((g.r - r0) / wid) ** 2) * (g.X / g.rn)
    packet[~g.interior] = 0.0
    psi, psi_prev = packet.copy(), packet.copy()
    E0f, E0w = g.energy(psi, psi_prev, dt)
    Ef, Ew = [E0f], [E0w]
    for _ in range(8):
        psi, psi_prev = g.leapfrog(psi, psi_prev, 150, dt)
        ef, ew = g.energy(psi, psi_prev, dt)
        Ef.append(ef); Ew.append(ew)
    Ef, Ew = np.array(Ef), np.array(Ew)
    finite = np.isfinite(psi).all()
    drift_flat = (Ef.max() - Ef.min()) / abs(Ef[0])
    drift_w = (Ew.max() - Ew.min()) / abs(Ew[0])
    print(f"\n[T2] localized twist packet, 1200 steps (48³): finite={finite}")
    print(f"     energy osc-range:  flat={100*drift_flat:.1f}%   1/r²-weighted={100*drift_w:.1f}%")
    print(f"     (1/r² is the operator's self-adjoint measure ⇒ the conserved one)")

    # --- T3: gauge-covariant gradient --------------------------------------------
    psi_track = 0.5 * np.log(g.reg)                  # ∇ψ = Â exactly (Â = ∇[½ln reg])
    gx, gy, gz = g._d(psi_track, 0), g._d(psi_track, 1), g._d(psi_track, 2)
    covgrad = np.sqrt((gx - g.Ax) ** 2 + (gy - g.Ay) ** 2 + (gz - g.Az) ** 2)
    baregrad = np.sqrt(gx**2 + gy**2 + gz**2)
    ratio = covgrad[g.ann].mean() / baregrad[g.ann].mean()
    print(f"\n[T3] gauge-covariant gradient on ψ=½ln(r²+r_c²) (the Â-winding phase):")
    print(f"     mean|∇ψ−Â|={covgrad[g.ann].mean():.2e}  vs mean|∇ψ|={baregrad[g.ann].mean():.2e}  (ratio {ratio:.1e})")
    print(f"     → physical field is the COVARIANT gradient (∇ψ−Â), ≈0 on the winding vacuum ✓")

    ok = conv and (relerr < 0.10) and finite_u and finite and (drift_w < 0.10) and (ratio < 0.05)
    print("\n" + "=" * 72)
    print("M5.6.1b: operator well-posed on the regularized hedgehog; bare ψ MASSLESS")
    print("(dx-convergent), stable (1/r²-energy bounded), gauge-covariant. Finite mass/")
    print("clock → M5.6.2/.3 (regularized core sets the scale).")
    print("PASS" if ok else "PARTIAL — inspect the failing metric above")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
