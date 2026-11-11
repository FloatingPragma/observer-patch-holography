"""Placement of the transition on the coupling chain's scale, and the cosmic time of that temperature."""
import math

HBAR_GEV_S = 6.582119569e-25      # CODATA
E_PLANCK_GEV = 1.22089e19         # display unit of the chain
E_CELL_GEV = 9.5599e18            # cell energy of the chain
ALPHA_S_MZ_CHAIN = 0.11834        # coupling chain, declared forms
LAMBDA3_CHAIN_GEV = 0.3400        # four-loop running, three-loop decoupling of ALPHA_S_MZ_CHAIN
LAMBDA3_CHAIN_CORNERS_GEV = (0.3251, 0.3547)

T_PC_MEV = (156.5, 1.5)           # HotQCD 2019
LAMBDA3_MEV = (338.0, 12.0)       # FLAG 2021
G_ABOVE = 61.75                   # gluons, u d s quarks, e, mu, three neutrinos, photon
G_BELOW = 17.25                   # pions, e, mu, three neutrinos, photon


def ratio(T=T_PC_MEV, lam=LAMBDA3_MEV):
    r = T[0] / lam[0]
    return r, r * math.hypot(T[1] / T[0], lam[1] / lam[0])


def placement(T=T_PC_MEV, lam=LAMBDA3_MEV):
    r, dr = ratio(T, lam)
    t = r * LAMBDA3_CHAIN_GEV * 1e3
    corners = [r * c * 1e3 for c in LAMBDA3_CHAIN_CORNERS_GEV]
    d_corner = max(abs(c - t) for c in corners)
    return dict(ratio=r, ratio_err=dr, T_MeV=t, err_ratio_MeV=dr * LAMBDA3_CHAIN_GEV * 1e3,
                err_T_import_MeV=t * T[1] / T[0], err_Lambda_import_MeV=t * lam[1] / lam[0],
                err_chain_corners_MeV=d_corner, err_total_MeV=math.hypot(dr * LAMBDA3_CHAIN_GEV * 1e3, d_corner),
                chain_corners_MeV=corners, Lambda_chain_over_Lambda_lattice=LAMBDA3_CHAIN_GEV * 1e3 / lam[0],
                T_over_E_cell=t * 1e-3 / E_CELL_GEV)


def cosmic_time_seconds(T_MeV, g_star):
    """Radiation era: t = sqrt(45 / (16 pi^3 g*)) E_Planck hbar / T^2."""
    T = T_MeV * 1e-3
    return math.sqrt(45.0 / (16.0 * math.pi ** 3 * g_star)) * E_PLANCK_GEV / T ** 2 * HBAR_GEV_S


def tick_seconds():
    return HBAR_GEV_S / E_CELL_GEV


def g_star(bosons, fermions):
    return bosons + 7.0 * fermions / 8.0


def table():
    p = placement()
    t_hi, t_lo = cosmic_time_seconds(p['T_MeV'], G_ABOVE), cosmic_time_seconds(p['T_MeV'], G_BELOW)
    alt = {'Borsanyi et al. 2020 with FLAG 2021': placement((158.0, 0.6), LAMBDA3_MEV),
           'HotQCD 2019 with FLAG 2024': placement(T_PC_MEV, (338.0, 10.0)),
           'HotQCD 2019 with ALPHA 2017': placement(T_PC_MEV, (341.0, 12.0))}
    return dict(inputs=dict(alpha_s_mZ=ALPHA_S_MZ_CHAIN, Lambda3_GeV=LAMBDA3_CHAIN_GEV,
                            Lambda3_corners_GeV=list(LAMBDA3_CHAIN_CORNERS_GEV), T_pc_MeV=list(T_PC_MEV),
                            Lambda3_lattice_MeV=list(LAMBDA3_MEV), E_cell_GeV=E_CELL_GEV, E_Planck_GeV=E_PLANCK_GEV),
                primary=p, alternatives={k: dict(ratio=v['ratio'], T_MeV=v['T_MeV'], err_ratio_MeV=v['err_ratio_MeV'])
                                         for k, v in alt.items()},
                g_above=G_ABOVE, g_below=G_BELOW, tick_seconds=tick_seconds(),
                t_seconds_g_above=t_hi, t_seconds_g_below=t_lo,
                ticks_g_above=t_hi / tick_seconds(), ticks_g_below=t_lo / tick_seconds(),
                pure_gauge=dict(Tc_over_sqrt_sigma_Boyd=[0.629, 0.003],
                                Tc_over_sqrt_sigma_LTW=[0.596 + 0.453 / 9.0, math.hypot(0.004, 0.030 / 9.0)],
                                Tc_r0=[0.7498, 0.0050], r0_Lambda0=[0.602, 0.048],
                                Tc_over_Lambda0=[0.7498 / 0.602, 0.7498 / 0.602 * math.hypot(0.0050 / 0.7498, 0.048 / 0.602)]))
