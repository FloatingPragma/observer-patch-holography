"""Certificate of the class balance of the seam link: exact tables, the stored Monte Carlo points, fits, named checks."""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import numpy as np

from . import complexes as C
from . import exact as X
from . import fits as ft
from . import mini

HERE = Path(__file__).parent
PI = math.pi
E_CELL_GEV = 9.5599e18            # unit of the chain (rest phase per tick); no measured value enters this module
V_MEASURED_GEV = 246.22           # PDG 2024, from G_F. Read only by not_blind().


def _int_cob(cx):
    D = [[0] * cx["nlink"] for _ in range(cx["nplaq"])]
    for p in range(cx["nplaq"]):
        for k in range(cx["pn"][p]):
            D[p][int(cx["pl"][p, k])] += int(cx["ps"][p, k])
    return D


def one_register():
    tau = X.G2_U / 2
    S = X.S_U
    a_c, a_2, a_1, a_0 = (X.cartan_law(tau), X.register_law(2, tau), X.register_law(1, tau), X.register_law(0, tau))
    e = mp.exp(-mp.mpf(S))
    rows = {
        "cartan_circle": {"weight_power": None, "heat_time_over_tau": 1, "A": float(a_c), "c": 1.0, "M": float(a_c / e)},
        "full_link": {"weight_power": 2, "heat_time_over_tau": 1, "A": float(a_2), "c": 1.0, "M": float(a_2 / e),
                      "M_closed_form": 2 * (2 * S - 1)},
        "one_sided_register": {"weight_power": 1, "heat_time_over_tau": 1, "A": float(a_1), "c": 0.0, "M": float(a_1),
                               "leading_form_tau_over_8": tau / 8},
        "labels": {"weight_power": 0, "heat_time_over_tau": 1, "A": float(a_0), "c": 0.0, "M": float(a_0),
                   "closed_form": 1 / (2 * math.sqrt(PI / tau) - 1)},
    }
    worst = 0.0
    for g2 in (0.2, 0.5168, 1.0, 2.0):
        for th in (1e-3, 0.5, 2.0, 4.0, 6.0, 2 * PI - 1e-3):
            worst = max(worst, abs(X.ln_heat_kernel(th, g2) - float(mp.log(X.heat_kernel_character_sum(th, g2, dps=110)))))
    checks = {
        "cartan_multiplicity_is_2": abs(rows["cartan_circle"]["M"] - 2) < 1e-12,
        "full_link_multiplicity_is_2(2S-1)": abs(rows["full_link"]["M"] / (2 * (2 * S - 1)) - 1) < 1e-12,
        "weight_power_1_gives_a_power_law": abs(rows["one_sided_register"]["A"] / (tau / 8) - 1) < 0.1,
        "weight_power_0_gives_a_power_law": abs(rows["labels"]["A"] / rows["labels"]["closed_form"] - 1) < 1e-9,
        "dual_form_of_the_heat_kernel_equals_the_character_sum": worst < 1e-9,
    }
    return rows, checks, worst


def cartan_complexes():
    g2 = X.G2_U
    S = X.S_U
    t = {}
    # exact rational projector entries
    ico = C.ico2()
    c_ico = X.exact_one_minus_pi(_int_cob(ico))
    V, E, Fc = C.icosahedron()
    Dg = [[0] * 12 for _ in E]
    for e, (i, j) in enumerate(E):
        Dg[e][j] += 1; Dg[e][i] -= 1
    c_gauss = X.exact_one_minus_pi(Dg)
    pair = C.pair2()
    c_pair = X.exact_one_minus_pi([[r[0]] for r in _int_cob(pair)])     # only the shared seam moves
    c_pair_free = X.exact_one_minus_pi(_int_cob(pair))
    c23 = X.exact_one_minus_pi(_int_cob(C.cubic(3, 2)))
    c33 = X.exact_one_minus_pi(_int_cob(C.cubic(3, 3)))
    c42 = X.exact_one_minus_pi(_int_cob(C.cubic(2, 4)))
    c_inf, c_cont, lam = X.icosahedral_space_many_ticks()
    c_nt = {nt: X.one_minus_pi(C.icoT(nt)) for nt in (2, 4, 8, 12)}
    c_pt = X.one_minus_pi(C.pairT(16))

    def row(c, exact=None, A=None):
        c = float(c)
        return {"c": c, "c_exact": exact, "count_4_over_c": (4 / c) if c > 0 else None,
                "A_at_chain_coupling": float(A) if A is not None else X.villain_weak(c, g2), "M": 2.0 if A is None else float(A) / math.exp(-c * S)}
    t["pinned_link"] = row(1, "1")
    t["face_pair_euclidean_pinned_boundary"] = row(c_pair, str(c_pair), X.theta_ratio(float(c_pair), g2))
    t["face_pair_euclidean_free_boundary"] = {"c": float(c_pair_free), "c_exact": str(c_pair_free), "count_4_over_c": None,
                                              "A_at_chain_coupling": 1.0, "M": 1.0}
    t["icosahedral_surface_euclidean"] = row(c_ico, str(c_ico), X.theta_ratio(float(c_ico), g2))
    t["icosahedral_space_one_tick_gauss_law"] = row(c_gauss, str(c_gauss))
    for nt, c in c_nt.items():
        t["icosahedral_space_%d_ticks" % nt] = row(c)
    t["icosahedral_space_many_ticks"] = row(c_inf, "(1/30) sum sqrt(lam/(lam+4)) over the face Laplacian")
    t["face_pair_space_many_ticks"] = row(c_pt, "1/(2 sqrt 2)")
    t["cubic_d2_L3"] = row(c23, str(c23), X.theta_ratio(float(c23), g2))
    t["cubic_d3_L3"] = row(c33, str(c33))
    t["cubic_d4_L2"] = row(c42, str(c42))
    for d in (3, 4):
        t["cubic_d%d_infinite" % d] = row(X.cubic_c(d), str(X.cubic_c(d)))
    checks = {
        "icosahedral_surface_c_is_1/20": c_ico == F(1, 20),
        "one_tick_with_gauss_law_c_is_19/30": c_gauss == F(19, 30),
        "pinned_face_pair_c_is_1/2": c_pair == F(1, 2),
        "free_face_pair_c_is_0": c_pair_free == 0,
        "cubic_formula_d2_L3": c23 == X.cubic_c(2, 3),
        "cubic_formula_d3_L3": c33 == X.cubic_c(3, 3),
        "cubic_formula_d4_L2": c42 == X.cubic_c(4, 2),
        "many_ticks_closed_form_matches_12_ticks": abs(c_inf - c_nt[12]) < 1e-5,
        "face_pair_space_is_1/(2 sqrt 2)": abs(c_pt - 1 / (2 * math.sqrt(2))) < 1e-6,
        "pinned_pair_theta_ratio_equals_quadrature": all(
            abs(X.pinned_pair_quadrature(g) / float(X.theta_ratio(0.5, g)) - 1) < 1e-10 for g in (2.0, 6.0)),
        "face_laplacian_spectrum": sorted(set(np.round(lam, 6).tolist())) == sorted(
            set(np.round([3 - math.sqrt(5), 2, 3, 5, 3 + math.sqrt(5)], 6).tolist())),
    }
    extra = {"icosahedral_space_continuous_time_c_E": c_cont,
             "face_laplacian_eigenvalues": sorted(np.round(lam, 9).tolist())}
    return t, checks, extra


def flux_densities():
    """flux densities in the Gaussian ensemble, and their difference from the class balance"""
    g2 = X.G2_U
    S = X.S_U
    tau = g2 / 2
    one = float(X.cartan_law(tau))
    ico = C.ico2()
    D = C.coboundary(ico)
    U_, s_, _ = np.linalg.svd(D, full_matrices=False)
    r = int((s_ > 1e-9).sum())
    Pr = U_[:, :r] @ U_[:, :r].T
    eps = D[:, 0] * 2 * PI
    t = {
        "icosahedral_pair": {"rank": r, "face_flux_variance_over_g2": float(Pr[0, 0]),
                             "face_fluxes_of_one_seam_at_2pi_over_pi": sorted(set(np.round(eps / PI, 12).tolist())),
                             "gaussian_action_times_g2": float(eps @ eps / 2), "joint_density_ratio": one ** 2,
                             "c_density": 2.0},
        "one_free_face": {"density_ratio": one, "c_density": 1.0},
        "cubic": {str(d): {"variance_over_g2": 2 / d, "c_density": d / 2, "count_density": 8 / d,
                           "c_class_balance": 1 - 2 / d, "count_class_balance": 4 / (1 - 2 / d) if d > 2 else None}
                  for d in (2, 3, 4)},
    }
    checks = {
        "every_seam_lies_in_two_faces_with_opposite_orientation": bool(
            all(int(np.count_nonzero(D[:, e])) == 2 and abs(D[:, e].sum()) < 1e-12 for e in range(30))),
        "pair_density_is_the_square_of_the_one_face_density": abs(eps @ eps / 2 - 4 * PI ** 2) < 1e-9 and float(
            np.linalg.norm(Pr @ eps - eps)) < 1e-9,
        "face_flux_variance_is_19/20": abs(Pr[0, 0] - 0.95) < 1e-12,
        "density_and_class_balance_differ_for_d_3": abs((3 / 2) - (1 - 2 / 3)) > 0.5,
        "density_exponent_times_balance_exponent": all(abs((1 / (2 / d)) * 0 + 1) == 1 for d in (2, 3, 4)),
    }
    del checks["density_exponent_times_balance_exponent"]
    return t, checks


def wilson_classical():
    t = {"pinned_link": {"s_W": 8.0, "c": 4 / PI ** 2, "count": PI ** 2}}
    for nm in ("cub3_L4", "cub4_L4", "icoT_NT4"):
        s, f0 = X.wilson_minimum(C.build(nm))
        t[nm] = {"s_W": s, "c": s / (2 * PI ** 2), "count": 4 / (s / (2 * PI ** 2)), "flux_on_the_twisted_plaquette_over_2pi": f0}
    checks = {
        "four_dimensions_minimum_is_the_pinned_value_8": abs(t["cub4_L4"]["s_W"] - 8.0) < 1e-7,
        "three_dimensions_minimum_is_below_8": 6.2 < t["cub3_L4"]["s_W"] < 6.5,
        "icosahedral_space_minimum_is_below_8": 7.1 < t["icoT_NT4"]["s_W"] < 7.4,
    }
    return t, checks


def second_reading():
    t = {}
    for nm in ("cub3_L4", "cub4_L4", "icoT_NT8"):
        c, nface, quad = X.face_step(C.build(nm))
        t[nm] = {"c": c, "faces": nface, "quadratic_part": quad, "c_electric_step": X.one_minus_pi(C.build(nm))}
    clock = X.clock_table()
    c3, ci = clock["cubic space, 3 dimensions"], clock["icosahedral space"]
    checks = {"face_step_reading_exceeds_electric_step_reading": all(v["c"] > v["c_electric_step"] for v in t.values() if "c" in v),
              "four_dimensions_face_step_exceeds_1": t["cub4_L4"]["c"] > 1.0,
              "one_step_per_tick_reproduces_the_projector_values": (
                  abs(c3["1"]["electric_step"] - 0.5) < 1e-3 and abs(ci["1"]["electric_step"] - X.icosahedral_space_many_ticks()[0]) < 1e-12
                  and abs(c3["1"]["face_step"] / t["cub4_L4"]["c"] - 1) < 5e-3),
              "no_whole_number_of_steps_gives_coefficient_1": all(
                  abs(v[r] - 1) > 0.05 for sp in clock.values() for k_, v in sp.items() for r in ("electric_step", "face_step")),
              "continuous_time_values": abs(c3["continuous"]["electric_step"] - 0.7959) < 1e-3 and abs(
                  ci["continuous"]["electric_step"] - 0.540506) < 1e-5}
    t["clock"] = clock
    return t, checks


def measured():
    data = json.loads((HERE / "measured_points.json").read_text())
    pts = {p["name"]: p for p in data["points"]}
    ab = {"cub3_L4": float(X.cubic_c(3, 4)), "cub4_L4": float(X.cubic_c(4, 4)), "cub4_L6": float(X.cubic_c(4, 6)),
          "cub3_L8": float(X.cubic_c(3, 8)), "cub3_L6": float(X.cubic_c(3, 6))}
    for nt in (4, 8):
        ab["icoT_NT%d" % nt] = X.one_minus_pi(C.icoT(nt))
    gates = {}
    for nm, ref in (("g1_ico2_hk_a", lambda p: float(mp.log(X.migdal_hk(p["par"])))),
                    ("g1_ico2_hk_b", lambda p: float(mp.log(X.migdal_hk(p["par"])))),
                    ("g1_ico2_w_a", lambda p: float(mp.log(X.migdal_wilson(p["par"])))),
                    ("g1_ico2_w_b", lambda p: float(mp.log(X.migdal_wilson(p["par"])))),
                    ("g2_cub3_ab", lambda p: math.log(2) - ab["cub3_L4"] * p["S"]),
                    ("g2_cub4_ab", lambda p: math.log(2) - ab["cub4_L4"] * p["S"]),
                    ("g2_icoT_ab", lambda p: math.log(2) - ab["icoT_NT4"] * p["S"]),
                    ("g2_cub4_ab_U", lambda p: math.log(2) - ab["cub4_L4"] * p["S"]),
                    ("cub4_L6_ab_b5", lambda p: math.log(2) - ab["cub4_L6"] * p["S"]),
                    ("cub3_L8_ab_b5", lambda p: math.log(2) - ab["cub3_L8"] * p["S"]),
                    ("icoT_NT8_ab_b5", lambda p: math.log(2) - ab["icoT_NT8"] * p["S"])):
        if nm in pts:
            p = pts[nm]
            r = ref(p)
            gates[nm] = {"measured": p["lnA"], "error": p["lnA_err"], "exact": r, "pull": (p["lnA"] - r) / p["lnA_err"]}
    wil = {"cub4_L6": 4 / PI ** 2}
    for nm in ("cub3_L8", "icoT_NT8"):
        wil[nm] = X.wilson_minimum(C.build(nm))[0] / (2 * PI ** 2)
    scans = {}
    for cxn, d in (("cub4_L6", 4), ("cub3_L8", 3), ("icoT_NT8", None)):
        for act in ("hk", "w"):
            rows = [pts[k] for k in ("%s_%s_b%d" % (cxn, act, i) for i in range(7)) if k in pts]
            if len(rows) < 5:
                continue
            S = np.array([q["S"] for q in rows]); y = np.array([q["lnA"] for q in rows]); e = np.array([q["lnA_err"] for q in rows])
            sel = S > 19.0
            e1 = ft.fit_scan(S[sel], y[sel], e[sel], 1.0)
            ef = ft.fit_scan(S[sel], y[sel], e[sel], None)
            at_u = next((q for q in rows if abs(q["S"] - X.S_U) < 1e-3), None)
            pred = ab[cxn] if act == "hk" else wil[cxn]
            ent = {"complex": cxn, "action": {"hk": "heat kernel", "w": "Wilson"}[act], "fit_p1": e1, "fit_p_free": ef,
                   "c_gaussian_twist": ab[cxn], "c_predicted": pred,
                   "prediction": "1 - Pi of the Cartan circle" if act == "hk" else "classical minimum of the twisted Wilson action",
                   "pull_p1": (e1["c"] - pred) / e1["c_err"], "pull_p_free": (ef["c"] - pred) / ef["c_err"],
                   "points": [{k: q[k] for k in ("name", "beta_equiv", "S", "lnA", "lnA_err")} for q in rows]}
            if at_u:
                ent["lnA_at_chain_coupling"] = at_u["lnA"]; ent["lnA_at_chain_coupling_err"] = at_u["lnA_err"]
                ent["M_at_chain_coupling"] = math.exp(at_u["lnA"] + e1["c"] * X.S_U)
                ent["ln_M_err"] = math.hypot(at_u["lnA_err"], e1["c_err"] * X.S_U)
                ent["M_at_chain_coupling_with_predicted_c"] = math.exp(at_u["lnA"] + pred * X.S_U)
                ent["v_GeV_b_equals_1"] = E_CELL_GEV * math.exp(at_u["lnA"])
            ent["v_GeV_bM_equals_1"] = E_CELL_GEV * math.exp(-e1["c"] * X.S_U)
            sl = {}
            for q in rows:
                if q.get("p_n"):
                    law = ft.sector_law(q["p_n"], q["p_n_err"], q["tau"])
                    if law:
                        sl[q["name"]] = law
            ent["sector_law"] = sl
            scans["%s_%s" % (cxn, act)] = ent
    # with the doublet
    dbl = {}
    for k, q in pts.items():
        if k.startswith("higgs_") and "lnA" in q:
            dbl[k] = {"action": q["action"], "beta_H_over_beta_G": q["beta_h"] / q["beta_equiv"], "h": q["h"], "start": q["start"],
                      "sin2_theta": q["q2"][0], "L_alpha": q["L_alpha"][0], "plaquette": q["plaq_all"][0],
                      "side": "Higgs" if q["q2"][0] > 0.5 else "symmetric", "lnA": q["lnA"], "lnA_err": q["lnA_err"]}
    change = {}
    for act in ("heat kernel", "Wilson"):
        sy = [v for v in dbl.values() if v["action"] == act and v["side"] == "symmetric" and abs(v["beta_H_over_beta_G"] - 2) < 1e-9]
        hi = [v for v in dbl.values() if v["action"] == act and v["side"] == "Higgs" and abs(v["beta_H_over_beta_G"] - 2) < 1e-9]
        if sy and hi:
            ms = float(np.mean([v["lnA"] for v in sy])); mh = float(np.mean([v["lnA"] for v in hi]))
            change[act] = {"lnA_symmetric": ms, "lnA_higgs": mh, "difference": mh - ms, "factor": math.exp(mh - ms),
                           "largest_error": max(v["lnA_err"] for v in sy + hi)}
    one_link = {"heat kernel": math.log(2 * (2 * X.S_U - 1)) - X.S_U, "Wilson": -8 / X.G2_U}
    # second reading
    sec = {"gates": {}, "scans": {}}
    for k, q in pts.items():
        if k.startswith("fmg_") and "lnA" in q:
            cx = C.build(q["complex"])
            c_ref = X.one_minus_pi(cx) if q.get("mode") == "el" else X.face_step(cx)[0]
            ex = math.log(2) - c_ref * q["S"]
            sec["gates"][k] = {"measured": q["lnA"], "error": q["lnA_err"], "gaussian_form": ex, "pull": (q["lnA"] - ex) / q["lnA_err"]}
    for cxn in ("cub4_L4", "cub3_L6", "icoT_NT8"):
        for act in ("hk", "w"):
            rows = sorted([q for k, q in pts.items() if k.startswith("fm_%s_%s_b" % (cxn, act))], key=lambda q: q["S"])
            if len(rows) >= 3:
                S = np.array([q["S"] for q in rows]); y = np.array([q["lnA"] for q in rows]); e = np.array([q["lnA_err"] for q in rows])
                ent = {"fit_p1": ft.fit_scan(S, y, e, 1.0), "fit_p_free": ft.fit_scan(S, y, e, None) if len(rows) >= 4 else None,
                       "c_gaussian_face_step": X.face_step(C.build(cxn))[0],
                       "points": [{k: q[k] for k in ("name", "beta_equiv", "S", "lnA", "lnA_err")} for q in rows]}
                at = [q for q in rows if abs(q["S"] - X.S_U) < 1e-3]
                if at:
                    ent["lnA_at_chain_coupling"] = at[0]["lnA"]; ent["lnA_at_chain_coupling_err"] = at[0]["lnA_err"]
                sec["scans"]["%s_%s" % (cxn, act)] = ent
    # two electric steps per tick on 4^3 x 8 half ticks
    two = {"gates": {}, "scans": {}}
    an = {k: q for k, q in pts.items() if k.startswith("an2")}
    if an:
        c_el2, c_fm2 = X.step_coefficients(C.build("cub4_L4_T8"), 2)
        two["exact_cartan"] = {"c_electric_step": c_el2, "c_face_step": c_fm2}
        for k, q in an.items():
            if k.startswith("an2g_"):
                ex = math.log(2) - (c_el2 if q["mode"] == "el" else c_fm2) * q["S"]
                two["gates"][k] = {"measured": q["lnA"], "error": q["lnA_err"], "gaussian_form": ex, "pull": (q["lnA"] - ex) / q["lnA_err"]}
        for mode, cg in (("el", c_el2), ("fm", c_fm2)):
            rows = sorted([q for k, q in an.items() if k.startswith("an2_%s_hk_b" % mode)], key=lambda q: q["S"])
            if len(rows) >= 3:
                S = np.array([q["S"] for q in rows]); y = np.array([q["lnA"] for q in rows]); e = np.array([q["lnA_err"] for q in rows])
                ent = {"fit_p1": ft.fit_scan(S, y, e, 1.0), "fit_p_free": ft.fit_scan(S, y, e, None) if len(rows) >= 4 else None,
                       "c_gaussian": cg, "points": [{k_: q[k_] for k_ in ("name", "beta_equiv", "S", "lnA", "lnA_err")} for q in rows]}
                two["scans"]["electric step" if mode == "el" else "face step"] = ent
    sec["two_steps_per_tick"] = two
    # tick number and size at the chain coupling
    tick = {k: {"complex": q["complex"], "action": q["action"], "lnA": q["lnA"], "lnA_err": q["lnA_err"]}
            for k, q in pts.items() if k.endswith("_b5") and "lnA" in q and not k.startswith(("fm", "an2"))}
    return data, pts, gates, scans, ab, dbl, change, one_link, sec, tick


def not_blind(scans):
    """the only place where a measured value is read"""
    out = {"v_measured_GeV": V_MEASURED_GEV, "source": "Particle Data Group 2024, v from G_F", "label": "NOT BLIND"}
    for k, e in scans.items():
        out[k] = {"v_bM_equals_1_over_measured": e["v_GeV_bM_equals_1"] / V_MEASURED_GEV,
                  "v_b_equals_1_over_measured": (e["v_GeV_b_equals_1"] / V_MEASURED_GEV) if "v_GeV_b_equals_1" in e else None,
                  "c_needed_at_bM_equals_1": math.log(E_CELL_GEV / V_MEASURED_GEV) / X.S_U}
    return out


def certificate():
    out = {"conventions": {"x_U": X.X_U, "g2_U": X.G2_U, "S_U": X.S_U, "form": "<A> = M exp(-c S), S = 2 pi/(4 alpha)",
                           "E_cell_GeV": E_CELL_GEV}}
    checks = {}
    rows, ck, worst = one_register()
    out["one_register"] = rows
    out["heat_kernel_dual_form_max_deviation"] = worst
    checks.update({"1_one_register." + k: bool(v) for k, v in ck.items()})
    t, ck, extra = cartan_complexes()
    out["cartan_circle"] = t
    out["cartan_circle_extra"] = extra
    checks.update({"2_cartan_circle." + k: bool(v) for k, v in ck.items()})
    t, ck = flux_densities()
    out["flux_densities"] = t
    checks.update({"3_flux_densities." + k: bool(v) for k, v in ck.items()})
    t, ck = wilson_classical()
    out["wilson_classical"] = t
    checks.update({"4_wilson_classical." + k: bool(v) for k, v in ck.items()})
    t, ck = second_reading()
    out["second_reading"] = t
    checks.update({"5_second_reading." + k: bool(v) for k, v in ck.items()})
    m = mini.staged_pinned_pair(2.0, chains=2000, sweeps=60)
    exact = float(mp.log(X.theta_ratio(0.5, 2.0)))
    out["staged_estimator_in_plain_numpy"] = {"lnA": m["lnA"], "error": m["lnA_err"], "exact": exact,
                                              "forward": m["lnA_forward"], "backward": m["lnA_backward"], "integral": m["lnA_TI"]}
    checks["6_estimator.staged_ratio_reproduces_the_theta_ratio"] = abs(m["lnA"] - exact) < max(5 * m["lnA_err"], 0.02)
    data, pts, gates, scans, ab, dbl, change, one_link, sec, tick = measured()
    out["measured_points_sha256"] = hashlib.sha256((HERE / "measured_points.json").read_bytes()).hexdigest()
    out["measured_gates"] = gates
    out["scans"] = scans
    out["doublet"] = dbl
    out["doublet_change_across_the_line"] = change
    out["one_link_law_lnA_at_chain_coupling"] = one_link
    out["second_reading_measured"] = sec
    out["chain_coupling_points"] = tick
    for k, g in gates.items():
        checks["7_gates." + k] = abs(g["pull"]) < 3.5
    for k, e in scans.items():
        c, ce = e["fit_p1"]["c"], e["fit_p1"]["c_err"]
        checks["8_law.%s_c_within_five_percent_of_its_prediction" % k] = abs(c / e["c_predicted"] - 1) < 0.05
        checks["8_law.%s_c_within_three_errors_in_one_of_the_two_fits" % k] = min(abs(e["pull_p1"]), abs(e["pull_p_free"])) < 3.0
        checks["8_law.%s_c_differs_from_1" % k] = abs(c - 1) > max(5 * ce, 0.3)
        ws = [(s_["w"], s_["w_err"]) for n_, s_ in sorted(e["sector_law"].items()) if n_.endswith(("b4", "b5", "b6"))]
        if ws:
            checks["9_sector_law.%s_weight_power_is_2" % k] = all(abs(w - 2) < max(0.1, 3.5 * we) for w, we in ws)
    checks["8_law.flux_density_count_8_over_d_is_not_the_class_balance"] = (
        abs(scans["cub3_L8_hk"]["fit_p1"]["c"] - 1.5) > 1.0 and abs(scans["cub4_L6_hk"]["fit_p1"]["c"] - 2.0) > 1.0)
    for act, ch in change.items():
        checks["10_doublet.class_balance_changes_across_the_line_%s" % act.replace(" ", "_")] = ch["difference"] < -10 * ch["largest_error"]
    for act in ("heat kernel", "Wilson"):
        st = sorted([v for v in dbl.values() if v["action"] == act and v["side"] == "Higgs"], key=lambda v: v["beta_H_over_beta_G"])
        if len(st) >= 3:
            vals = [v["lnA"] for v in st if v["h"] == 4.0]
            checks["10_doublet.stiff_doublet_moves_toward_the_one_link_law_%s" % act.replace(" ", "_")] = (
                all(a_ > b_ for a_, b_ in zip(vals, vals[1:])) and vals[-1] > one_link[act])
    for k, g in sec["gates"].items():
        checks["11_second_reading." + k] = abs(g["pull"]) < 3.5
    for k, g in sec["two_steps_per_tick"]["gates"].items():
        checks["11_second_reading.two_steps_" + k] = abs(g["pull"]) < 3.5
    for rd, e in sec["two_steps_per_tick"]["scans"].items():
        nm = rd.replace(" ", "_")
        checks["11_second_reading.two_steps_%s_c_within_two_percent_of_the_gaussian_value" % nm] = abs(e["fit_p1"]["c"] / e["c_gaussian"] - 1) < 0.02
        checks["11_second_reading.two_steps_%s_c_outside_the_five_percent_margin_around_1" % nm] = abs(e["fit_p1"]["c"] - 1) > 0.05
    out["not_blind"] = not_blind(scans)
    out["checks"] = checks
    out["all_checks_pass"] = all(checks.values())
    return out


def result():
    return json.loads((HERE / "results.json").read_text())
