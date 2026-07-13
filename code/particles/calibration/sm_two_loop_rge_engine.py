"""Standard Model two-loop renormalization group equations in MS-bar.

Couplings evolved: g1 (hypercharge, GUT normalization g'^2 = (3/5) g1^2),
g2 (SU(2)_L), g3 (SU(3)_c), yt (top Yukawa), lam (Higgs quartic).
All Yukawa couplings other than yt are set to zero (yb = ytau = 0).

Convention: V = lam * (H^dag H)^2, so mH_tree^2 = 2 * lam * v^2 and
lam(mt) ~ 0.126 for Mh = 125 GeV. Sources that use V = (lam'/2)(phi^dag phi)^2
(e.g. Luo and Xiao) relate to this convention by lam = lam'/2, with every
beta coefficient converted accordingly.

Beta functions are d/dln(mu). Sources:

  [1] Buttazzo, Degrassi, Giardino, Giudice, Sala, Salvio, Strumia,
      arXiv:1307.3536, appendix B. Printed there as d(g_i^2)/dln(mubar^2),
      d(yt^2)/dln(mubar^2) and d(lam)/dln(mubar^2); converted here via
      beta_g = RHS/g, beta_yt = RHS/yt, beta_lam = 2*RHS. Their potential
      is V = -m^2/2 |H|^2 + lam |H|^4, identical to the convention above.

  [2] Luo and Xiao, arXiv:hep-ph/0207271 (Phys.Rev.Lett. 90 (2003) 011601),
      complete two-loop SM RGEs. Their quartic is (lam_LX/2)(phi^dag phi)^2,
      i.e. lam_LX = 2*lam. After conversion, their beta_yt and beta_lam
      two-loop expressions (with H^dag H -> yt^2, F_D = F_L = 0, n_g = 3)
      agree with [1] term by term. Their gauge two-loop coefficients match
      the Machacek-Vaughn values quoted in [1].

One-loop reference set (reproduced exactly by beta_1loop):

  16pi^2 b_g1  = (41/10) g1^3
  16pi^2 b_g2  = (-19/6) g2^3
  16pi^2 b_g3  = -7 g3^3
  16pi^2 b_yt  = yt*( (9/2) yt^2 - (17/20) g1^2 - (9/4) g2^2 - 8 g3^2 )
  16pi^2 b_lam = 24 lam^2 - 6 yt^4
                 + (3/8)*( 2 g2^4 + (g2^2 + (3/5) g1^2)^2 )
                 + lam*( 12 yt^2 - 9 g2^2 - (9/5) g1^2 )

Two-loop terms (cross-checked between [1] and [2]):

  (16pi^2)^2 b_g1^(2) = g1^3 * ( (199/50) g1^2 + (27/10) g2^2
                                 + (44/5) g3^2 - (17/10) yt^2 )
  (16pi^2)^2 b_g2^(2) = g2^3 * ( (9/10) g1^2 + (35/6) g2^2
                                 + 12 g3^2 - (3/2) yt^2 )
  (16pi^2)^2 b_g3^(2) = g3^3 * ( (11/10) g1^2 + (9/2) g2^2
                                 - 26 g3^2 - 2 yt^2 )
  (16pi^2)^2 b_yt^(2) = yt * ( -12 yt^4 - 12 lam yt^2 + 6 lam^2
        + yt^2*( (393/80) g1^2 + (225/16) g2^2 + 36 g3^2 )
        + (1187/600) g1^4 - (9/20) g1^2 g2^2 + (19/15) g1^2 g3^2
        - (23/4) g2^4 + 9 g2^2 g3^2 - 108 g3^4 )
  (16pi^2)^2 b_lam^(2) = -312 lam^3 - 144 lam^2 yt^2
        + lam^2*( 108 g2^2 + (108/5) g1^2 )
        + lam yt^2*( -3 yt^2 + 80 g3^2 + (45/2) g2^2 + (17/2) g1^2 )
        + lam*( -(73/8) g2^4 + (1887/200) g1^4 + (117/20) g1^2 g2^2 )
        + yt^4*( 30 yt^2 - 32 g3^2 - (8/5) g1^2 )
        + yt^2*( -(9/4) g2^4 - (171/100) g1^4 + (63/10) g1^2 g2^2 )
        + (305/16) g2^6 - (3411/2000) g1^6
        - (289/80) g2^4 g1^2 - (1677/400) g2^2 g1^4

Validation benchmark (Buttazzo et al., Mh = 125.15, Mt = 173.34,
alpha_s(MZ) = 0.1184): NNLO-matched couplings at mu = Mt integrated to
mu = M_Pl = 1.2209e19 GeV and compared to their published 3-loop endpoint.
Run `python sm_two_loop_rge.py` to execute validate().
"""

import math

PI = math.pi
K1 = 1.0 / (16.0 * PI ** 2)
K2 = K1 * K1

MT = 173.34
MPL = 1.2209e19

# NNLO-matched couplings at mu = Mt (Buttazzo et al., Mh = 125.15,
# Mt = 173.34, alpha_s(MZ) = 0.1184). The paper quotes the hypercharge
# coupling as g_Y(Mt) = 0.35830 (their eq. 60); the state variable here is
# the GUT-normalized g1 = sqrt(5/3) * g_Y. Their Planck-scale table
# (eq. 64) quotes g_1 in GUT normalization directly.
GY_AT_MT = 0.35830
COUPLINGS_AT_MT = (math.sqrt(5.0 / 3.0) * GY_AT_MT,
                   0.64779, 1.1666, 0.93690, 0.12604)

# Published 3-loop endpoint at mu = M_Pl from the same reference
# (g1 GUT-normalized).
BENCHMARK_AT_MPL = (0.6154, 0.5055, 0.4873, 0.3825, -0.0113)


def beta_1loop(g1, g2, g3, yt, lam):
    """One-loop beta functions d/dln(mu) for (g1, g2, g3, yt, lam)."""
    g1s, g2s, g3s, yts = g1 * g1, g2 * g2, g3 * g3, yt * yt
    b_g1 = K1 * (41.0 / 10.0) * g1 ** 3
    b_g2 = K1 * (-19.0 / 6.0) * g2 ** 3
    b_g3 = K1 * (-7.0) * g3 ** 3
    b_yt = K1 * yt * (4.5 * yts - (17.0 / 20.0) * g1s
                      - (9.0 / 4.0) * g2s - 8.0 * g3s)
    b_lam = K1 * (24.0 * lam * lam - 6.0 * yts * yts
                  + (3.0 / 8.0) * (2.0 * g2s ** 2
                                   + (g2s + 0.6 * g1s) ** 2)
                  + lam * (12.0 * yts - 9.0 * g2s - 1.8 * g1s))
    return (b_g1, b_g2, b_g3, b_yt, b_lam)


def beta_2loop(g1, g2, g3, yt, lam):
    """One-loop plus two-loop beta functions d/dln(mu).

    Two-loop coefficients transcribed from Buttazzo et al. arXiv:1307.3536
    appendix B and cross-checked against Luo-Xiao hep-ph/0207271 after the
    lam_LX = 2*lam conversion.
    """
    g1s, g2s, g3s, yts = g1 * g1, g2 * g2, g3 * g3, yt * yt
    b1 = beta_1loop(g1, g2, g3, yt, lam)

    b2_g1 = g1 ** 3 * ((199.0 / 50.0) * g1s + (27.0 / 10.0) * g2s
                       + (44.0 / 5.0) * g3s - (17.0 / 10.0) * yts)
    b2_g2 = g2 ** 3 * ((9.0 / 10.0) * g1s + (35.0 / 6.0) * g2s
                       + 12.0 * g3s - (3.0 / 2.0) * yts)
    b2_g3 = g3 ** 3 * ((11.0 / 10.0) * g1s + (9.0 / 2.0) * g2s
                       - 26.0 * g3s - 2.0 * yts)
    b2_yt = yt * (-12.0 * yts * yts - 12.0 * lam * yts + 6.0 * lam * lam
                  + yts * ((393.0 / 80.0) * g1s + (225.0 / 16.0) * g2s
                           + 36.0 * g3s)
                  + (1187.0 / 600.0) * g1s ** 2
                  - (9.0 / 20.0) * g1s * g2s
                  + (19.0 / 15.0) * g1s * g3s
                  - (23.0 / 4.0) * g2s ** 2
                  + 9.0 * g2s * g3s
                  - 108.0 * g3s ** 2)
    b2_lam = (-312.0 * lam ** 3 - 144.0 * lam * lam * yts
              + lam * lam * (108.0 * g2s + (108.0 / 5.0) * g1s)
              + lam * yts * (-3.0 * yts + 80.0 * g3s
                             + (45.0 / 2.0) * g2s + (17.0 / 2.0) * g1s)
              + lam * (-(73.0 / 8.0) * g2s ** 2
                       + (1887.0 / 200.0) * g1s ** 2
                       + (117.0 / 20.0) * g1s * g2s)
              + yts * yts * (30.0 * yts - 32.0 * g3s - (8.0 / 5.0) * g1s)
              + yts * (-(9.0 / 4.0) * g2s ** 2
                       - (171.0 / 100.0) * g1s ** 2
                       + (63.0 / 10.0) * g1s * g2s)
              + (305.0 / 16.0) * g2s ** 3
              - (3411.0 / 2000.0) * g1s ** 3
              - (289.0 / 80.0) * g2s ** 2 * g1s
              - (1677.0 / 400.0) * g2s * g1s ** 2)

    return (b1[0] + K2 * b2_g1,
            b1[1] + K2 * b2_g2,
            b1[2] + K2 * b2_g3,
            b1[3] + K2 * b2_yt,
            b1[4] + K2 * b2_lam)


def run(couplings, lnmu0, lnmu1, n_steps=8000, loops=2,
        return_trajectory=False):
    """Integrate the RGEs from lnmu0 to lnmu1 with fixed-step RK4.

    couplings: iterable (g1, g2, g3, yt, lam) at lnmu0.
    loops: 1 or 2, selecting beta_1loop or beta_2loop.
    return_trajectory: if True, return (endpoint, trajectory) where
    trajectory is a list of (lnmu, couplings) samples including both ends.
    """
    beta = beta_2loop if loops == 2 else beta_1loop
    y = tuple(float(c) for c in couplings)
    h = (lnmu1 - lnmu0) / n_steps
    traj = [(lnmu0, y)] if return_trajectory else None
    x = lnmu0
    for _ in range(n_steps):
        k1 = beta(*y)
        k2 = beta(*(y[i] + 0.5 * h * k1[i] for i in range(5)))
        k3 = beta(*(y[i] + 0.5 * h * k2[i] for i in range(5)))
        k4 = beta(*(y[i] + h * k3[i] for i in range(5)))
        y = tuple(y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i]
                                      + 2.0 * k3[i] + k4[i])
                  for i in range(5))
        x += h
        if return_trajectory:
            traj.append((x, y))
    if return_trajectory:
        return y, traj
    return y


def _lambda_zero_crossing(traj):
    """Scale mu (GeV) of the first lam sign change along a trajectory.

    Linear interpolation in lnmu between the bracketing samples. Returns
    None if lam does not change sign.
    """
    for (x0, y0), (x1, y1) in zip(traj, traj[1:]):
        l0, l1 = y0[4], y1[4]
        if l0 == 0.0:
            return math.exp(x0)
        if l0 * l1 < 0.0:
            x = x0 + (x1 - x0) * l0 / (l0 - l1)
            return math.exp(x)
    return None


def validate():
    """Certification of the transcription against published benchmarks.

    Integrates the two-loop system from mu = Mt = 173.34 GeV with the
    NNLO-matched Buttazzo et al. inputs up to M_Pl = 1.2209e19 GeV and
    compares to their published 3-loop endpoint. Also verifies that
    beta_1loop matches the reference one-loop set at a test point and
    that lam crosses zero between 1e9 and 1e12 GeV.
    """
    lnmu0 = math.log(MT)
    lnmu1 = math.log(MPL)
    end, traj = run(COUPLINGS_AT_MT, lnmu0, lnmu1, n_steps=8000, loops=2,
                    return_trajectory=True)

    names = ("g1", "g2", "g3", "yt", "lam")
    residuals = {n: end[i] - BENCHMARK_AT_MPL[i] for i, n in enumerate(names)}
    tolerances = {"g1": 0.005, "g2": 0.005, "g3": 0.005,
                  "yt": 0.010, "lam": 0.004}
    endpoint_ok = {n: abs(residuals[n]) < tolerances[n] for n in names}

    mu_zero = _lambda_zero_crossing(traj)
    zero_ok = mu_zero is not None and 1e9 < mu_zero < 1e12

    # One-loop reference check at an asymmetric test point.
    g1, g2, g3, yt, lam = 0.44, 0.61, 1.05, 0.87, 0.11
    got = beta_1loop(g1, g2, g3, yt, lam)
    ref = (
        K1 * (41.0 / 10.0) * g1 ** 3,
        K1 * (-19.0 / 6.0) * g2 ** 3,
        K1 * (-7.0) * g3 ** 3,
        K1 * yt * ((9.0 / 2.0) * yt ** 2 - (17.0 / 20.0) * g1 ** 2
                   - (9.0 / 4.0) * g2 ** 2 - 8.0 * g3 ** 2),
        K1 * (24.0 * lam ** 2 - 6.0 * yt ** 4
              + (3.0 / 8.0) * (2.0 * g2 ** 4
                               + (g2 ** 2 + 0.6 * g1 ** 2) ** 2)
              + lam * (12.0 * yt ** 2 - 9.0 * g2 ** 2
                       - (9.0 / 5.0) * g1 ** 2)),
    )
    oneloop_ok = all(abs(a - b) <= 1e-15 * max(1.0, abs(b))
                     for a, b in zip(got, ref))

    return {
        "sources": {
            "gauge_2loop": "Buttazzo et al. 1307.3536 app. B; matches "
                           "Machacek-Vaughn coefficients",
            "yt_2loop": "Buttazzo et al. 1307.3536 app. B; cross-checked "
                        "vs Luo-Xiao hep-ph/0207271 (lam_LX = 2*lam)",
            "lam_2loop": "Buttazzo et al. 1307.3536 app. B; cross-checked "
                         "vs Luo-Xiao hep-ph/0207271 (lam_LX = 2*lam)",
            "benchmark": "Buttazzo et al. 1307.3536 couplings at Mt and "
                         "M_Pl (Mh=125.15, Mt=173.34, alpha_s=0.1184)",
        },
        "inputs_at_mt": dict(zip(names, COUPLINGS_AT_MT)),
        "endpoint_at_mpl": dict(zip(names, end)),
        "benchmark_at_mpl": dict(zip(names, BENCHMARK_AT_MPL)),
        "residuals": residuals,
        "tolerances": tolerances,
        "endpoint_ok": endpoint_ok,
        "lambda_zero_crossing_gev": mu_zero,
        "lambda_zero_crossing_ok": zero_ok,
        "oneloop_reference_ok": oneloop_ok,
        "all_ok": all(endpoint_ok.values()) and zero_ok and oneloop_ok,
    }


if __name__ == "__main__":
    result = validate()
    for key, value in result.items():
        print(f"{key}: {value}")
    if not result["all_ok"]:
        raise SystemExit(1)
