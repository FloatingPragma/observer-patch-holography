"""Derived numbers of the electroweak transition from the production record of the lattice runs.

`production.json` holds the lattice measurements (positions and heights of susceptibility maxima, reweighted
curves, validation runs). This module maps them to temperatures through the dimensional reduction, extrapolates
the peak position in the lattice spacing, and evaluates the named checks.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import lattice as lt
from . import reduction as rd

HERE = Path(__file__).parent

# ---- inputs -----------------------------------------------------------------------------------------------
G2_DECLARED = 4.0 * math.pi / 29.605
ROWS = {
    "declared": dict(
        label="declared chain (slopes 33/5, 1, -3)", status="declared forms; derivation open",
        v=246.60, g2=G2_DECLARED, sin2=0.23074, m_H=125.84, m_t=171.99, alpha_s=0.11834,
        mass_assumption="m_H and m_t are the declared midpoint values"),
    "one_doublet_A": dict(
        label="source as built, one doublet; same quartic and top couplings as the declared row",
        status="source as built; Higgs and top masses assumed",
        v=222.9, g2=4.0 * math.pi / (142.25 * 0.2056), sin2=0.2056,
        m_H=125.84 * 222.9 / 246.60, m_t=171.99 * 222.9 / 246.60, alpha_s=0.1123,
        mass_assumption="m_H / v and m_t / v equal to the declared row"),
    "one_doublet_B": dict(
        label="source as built, one doublet; same Higgs and top masses as the declared row",
        status="source as built; Higgs and top masses assumed",
        v=222.9, g2=4.0 * math.pi / (142.25 * 0.2056), sin2=0.2056, m_H=125.84, m_t=171.99, alpha_s=0.1123,
        mass_assumption="m_H = 125.84 GeV and m_t = 171.99 GeV"),
}
STANDARD_MODEL = dict(rd.SM)
REFERENCE = dict(T_c=159.5, err=1.5, text_value=159.6, text_stat=0.1, cite="D'Onofrio, Rummukainen, Phys. Rev. D93 (2016) 025003",
                 peak_positions={"6": 158.73, "9": 159.01, "16": 159.26}, continuum_peak=159.58,
                 map_curve=dict(T=160.0, x=0.2899, y=-0.0536, g32_over_T=0.39282, T_y0=162.1))
# susceptibility for Standard Model inputs at beta_G = 9, read from the source of figure 6 of arXiv:1508.07161
# (square symbols); chi = V T <(phi^dag phi - <phi^dag phi>)^2>
PUBLISHED_CURVE = ((156.0, 6.22), (157.0, 6.84), (158.0, 8.22), (158.5, 9.24), (158.8, 9.94), (159.0, 10.02),
                   (159.2, 9.96), (159.5, 8.75), (160.0, 6.06), (160.5, 4.30), (161.0, 3.12))
# y(T) of the published Standard Model map, read from the source of figure 1 of arXiv:1508.07161
PUBLISHED_Y = ((158.5, -0.0932), (159.0, -0.0801), (159.5, -0.0666), (160.0, -0.0536), (161.0, -0.0281), (162.0, -0.0025))
X_C = dict(value=0.0983, err=0.0015, m_H_SM=72.0, m_H_SM_err=2.0, cite="Rummukainen et al., Nucl. Phys. B532 (1998) 283")
TABLE4 = ((60, 5, 0.358495, 137.534), (60, 8, 0.3479735, 137.669), (60, 12, 0.3426840, 137.842),
          (60, 20, 0.3387418, 138.019), (35, 8, 0.3450806, 94.181), (70, 8, 0.3491523, 153.620))


def inputs_of(row: dict) -> dict:
    return rd.chain_inputs(row["v"], row["g2"], row["sin2"], row["m_H"], row["m_t"], row["alpha_s"])


def load_production() -> dict:
    return json.loads((HERE / "production.json").read_text())


def _wls(design: np.ndarray, values: np.ndarray, errors: np.ndarray):
    a = design / errors[:, None]
    cov = np.linalg.inv(a.T @ a)
    coef = cov @ a.T @ (values / errors)
    chi2 = float(((a @ coef - values / errors) ** 2).sum())
    return coef, cov, chi2


def continuum_fit(points: list[dict]) -> dict:
    """Global fit  value = c0 + C / beta_G + D(size)  over sets with physical size N * 4 / beta_G.
    The largest size is the reference (D = 0). With one lattice spacing the slope is absent."""
    sizes = sorted({round(p["size"], 6) for p in points})
    ref = sizes[-1]
    others = [s for s in sizes if s != ref]
    with_slope = len({p["beta_G"] for p in points}) >= 2
    ncol = 1 + (1 if with_slope else 0) + len(others)
    first = 2 if with_slope else 1
    design = np.zeros((len(points), ncol))
    for i, p in enumerate(points):
        design[i, 0] = 1.0
        if with_slope:
            design[i, 1] = 1.0 / p["beta_G"]
        s = round(p["size"], 6)
        if s != ref:
            design[i, first + others.index(s)] = 1.0
    vals = np.array([p["value"] for p in points])
    errs = np.array([p["err"] for p in points])
    dof = len(points) - ncol
    if dof < 0:
        raise ValueError("too few sets for the continuum fit")
    coef, cov, chi2 = _wls(design, vals, errs)
    scale = math.sqrt(max(chi2 / dof, 1.0)) if dof > 0 else 1.0
    return dict(continuum=float(coef[0]), err=float(math.sqrt(cov[0, 0]) * scale),
                slope=float(coef[1]) if with_slope else None,
                slope_err=float(math.sqrt(cov[1, 1]) * scale) if with_slope else None, chi2=chi2, dof=dof,
                extrapolated=with_slope, size_reference=ref,
                size_offsets={str(s): [float(coef[first + k]), float(math.sqrt(cov[first + k, first + k]) * scale)]
                              for k, s in enumerate(others)},
                error_scale=scale, points=points)


def width_at_fraction(T: list, chi: list, fraction: float = 0.8):
    """Interval of T on which the curve exceeds the given fraction of its maximum (linear interpolation)."""
    T = np.asarray(T, float); chi = np.asarray(chi, float)
    i = int(np.argmax(chi)); level = fraction * chi[i]
    lo = hi = None
    for j in range(i, 0, -1):
        if chi[j - 1] < level <= chi[j]:
            lo = T[j - 1] + (level - chi[j - 1]) / (chi[j] - chi[j - 1]) * (T[j] - T[j - 1]); break
    for j in range(i, len(T) - 1):
        if chi[j + 1] < level <= chi[j]:
            hi = T[j] + (chi[j] - level) / (chi[j] - chi[j + 1]) * (T[j + 1] - T[j]); break
    return (None if lo is None else float(lo)), (None if hi is None else float(hi))


def certificate(production: dict | None = None) -> dict:
    prod = load_production() if production is None else production
    out: dict = dict(schema="oph.electroweak_crossover.v1")
    out["scope"] = ("Phase structure of the dimensionally reduced SU(2) + Higgs theory and the pseudo-critical "
                    "temperature of the electroweak transition for the couplings of the declared chain. The reduction "
                    "is imported and correct to relative order g^4. The chain's couplings are declared forms. "
                    "The hypercharge gauge field is not simulated.")
    checks: dict[str, bool] = {}

    # ---- 1. constants and lattice relations ----------------------------------------------------------------
    out["constants"] = dict(Sigma=lt.SIGMA, zeta=lt.ZETA, delta=lt.DELTA, c_B=rd.C_B, c_F=rd.C_F, c=rd.C_CONST)
    checks["1_constants.c_B"] = abs(rd.C_B - 1.953808) < 2e-6
    checks["1_constants.c_F"] = abs(rd.C_F - 0.567514) < 2e-6
    checks["1_constants.c"] = abs(rd.C_CONST + 0.348725) < 5e-6
    checks["1_constants.condensate"] = abs(lt.CONDENSATE_CONSTANT - 0.66796) < 2e-5
    rows4 = []
    for m_h, beta_G, beta_H, t_c in TABLE4:
        h2 = (m_h / 80.6) ** 2
        x = -0.00550 + 0.12622 * h2
        _, y = lt.continuum_parameters(beta_G, beta_H, x * beta_H ** 2 / beta_G)
        t = math.sqrt(2.58088 * m_h ** 2 / (0.39818 + 0.15545 * h2 - 0.00190 * h2 ** 2 - y))
        rows4.append(dict(m_H_star=m_h, beta_G=beta_G, beta_H_c=beta_H, x=x, y=y, T_from_relation=t, T_published=t_c))
        checks[f"2_relations.table4_mH{m_h}_betaG{beta_G}"] = abs(t - t_c) < 0.03
    out["lattice_relations_table4"] = rows4

    # ---- 2. the map on Standard Model inputs against the published curves ----------------------------------
    r = rd.reduce(160.0, **STANDARD_MODEL, linear=True)
    t0 = rd.solve_T(lambda x: 0.0, STANDARD_MODEL, linear=True)
    ref = REFERENCE["map_curve"]
    out["map_check_standard_model"] = dict(inputs=STANDARD_MODEL, x=r["x"], y=r["y"], g32_over_T=r["g32"] / 160.0, z=r["z"],
                                           T_y0=t0, published=ref,
                                           tree_level_T_y0=rd.solve_T(lambda x: 0.0, STANDARD_MODEL, tree=True))
    checks["3_map.x_at_160"] = abs(r["x"] - ref["x"]) < 0.003
    checks["3_map.g32_at_160"] = abs(r["g32"] / 160.0 - ref["g32_over_T"]) < 0.001
    checks["3_map.T_y0"] = abs(t0 - ref["T_y0"]) < 0.6
    a = rd.reduce(160.0, **STANDARD_MODEL, linear=True, mu_over_T=math.pi)
    b = rd.reduce(160.0, **STANDARD_MODEL, linear=True, mu_over_T=4 * math.pi)
    checks["3_map.scale_independent"] = abs(a["y"] - b["y"]) < 1e-10 and abs(a["x"] - b["x"]) < 1e-12
    # running of nu^2 from the Higgs self-energy bracket against the printed one-loop running
    sm = STANDARD_MODEL
    g02 = 4 * math.sqrt(2.0) * sm["G_mu"] * sm["mW"] ** 2
    h2 = (sm["mH"] / sm["mW"]) ** 2
    t2 = (sm["mt"] / sm["mW"]) ** 2
    s2 = 4 * math.pi * sm["alpha_s"] / g02
    b1 = rd.vacuum_brackets(100.0, sm["mW"], sm["mH"], sm["mt"], s2)
    b2 = rd.vacuum_brackets(100.0 * math.e, sm["mW"], sm["mH"], sm["mt"], s2)
    lam, gy2 = g02 * h2 / 8.0, g02 * t2 / 2.0
    checks["3_map.running_nu2"] = abs((b2[1] - b1[1]) * g02 - 2.0 * (-2.25 * g02 + 6 * lam + 3 * gy2)) < 1e-9
    checks["3_map.running_g2"] = abs((b2[0] - b1[0]) - 2.0 * (4 * rd.NF / 3.0 - 43 / 6.0)) < 1e-9
    # QCD part of the top relation at mu = m_t: the square of the pole to MS-bar mass relation
    d1 = rd.vacuum_brackets(sm["mt"], sm["mW"], sm["mH"], sm["mt"], s2)[3]
    d0 = rd.vacuum_brackets(sm["mt"], sm["mW"], sm["mH"], sm["mt"], 0.0)[3]
    checks["3_map.top_qcd_term"] = abs(2 * (d1 - d0) * g02 * rd.LOOP + (8.0 / 3.0) * sm["alpha_s"] / math.pi) < 1e-9

    # ---- 3. lattice production: peak positions per set ------------------------------------------------------
    sets = prod["chain_point"]
    lat = {}
    for key, s in sets.items():
        ys = s["yscan"]
        lat[key] = dict(n_side=s["L"], beta_G=s["bG"], x=ys["x"], size=s["L"] * 4.0 / s["bG"], runs=s["nruns"],
                        y_max=ys["y_max"], y_max_err=ys["y_max_err"], chi3_max=ys["chi3_max"],
                        chi3_max_err=ys["chi3_max_err"], half_width_y=ys["half_width_sym"],
                        T_max={row: [t["T_max"], t["T_max_err"]] for row, t in s["tscan"].items()},
                        chi_max={row: [t["chi_max"], t["chi_max_err"]] for row, t in s["tscan"].items()},
                        half_width_T={row: t["half_width_sym"] for row, t in s["tscan"].items()})
    out["lattice_sets"] = lat
    main_x = [v for k, v in lat.items() if abs(v["x"] - 0.296) < 1e-9]
    fit_y = continuum_fit([dict(beta_G=v["beta_G"], size=v["size"], n_side=v["n_side"], value=v["y_max"], err=v["y_max_err"])
                           for v in main_x])
    fit_T = continuum_fit([dict(beta_G=v["beta_G"], size=v["size"], n_side=v["n_side"], value=v["T_max"]["declared"][0],
                                err=v["T_max"]["declared"][1]) for v in main_x if "declared" in v["T_max"]])
    out["continuum_extrapolation"] = dict(y_max=fit_y, T_max_declared=fit_T,
                                          form="value = continuum + C / beta_G + D(physical size); largest size is the reference")
    checks["4_lattice.fit_T_dof_nonnegative"] = fit_T["dof"] >= 0
    checks["4_lattice.fit_T_chi2"] = fit_T["dof"] == 0 or fit_T["chi2"] / fit_T["dof"] < 4.0

    # reference set: 24^3 at beta_G = 8, the set on which every row is reweighted
    key24 = next((k for k, v in lat.items() if v["beta_G"] == 8 and v["n_side"] == 24 and abs(v["x"] - 0.296) < 1e-9), None)
    key36 = next((k for k, v in lat.items() if abs(v["x"] - 0.36) < 1e-9), None)
    ref_set = lat[key24]
    dy_cont = fit_y["continuum"] - ref_set["y_max"]           # shift of the peak from this set to the continuum
    out["reference_set"] = dict(n_side=ref_set["n_side"], beta_G=8, y_max=ref_set["y_max"], y_shift_to_continuum=dy_cont,
                                y_shift_err=fit_y["err"])

    # ---- 4. temperatures per row -----------------------------------------------------------------------------
    rows_out = {}
    for name, row in list(ROWS.items()) + [("standard_model", None)]:
        pars = dict(STANDARD_MODEL) if row is None else inputs_of(row)
        entry = dict(inputs=pars)
        if row is not None:
            entry.update(label=row["label"], status=row["status"], mass_assumption=row["mass_assumption"],
                         v=row["v"], g2=row["g2"], gprime2=row["g2"] * row["sin2"] / (1 - row["sin2"]), sin2=row["sin2"],
                         lambda_tree=row["g2"] * (row["m_H"] / pars["mW"]) ** 2 / 8.0,
                         top_yukawa_tree=math.sqrt(row["g2"] / 2.0) * row["m_t"] / pars["mW"])
        entry["T_y0"] = rd.solve_T(lambda x: 0.0, pars, linear=True)
        entry["curve"] = [dict(T=T, **{k: rd.reduce(T, **pars, linear=True)[k] for k in ("x", "y", "z")},
                               g32_over_T=rd.reduce(T, **pars, linear=True)["g32"] / T)
                          for T in np.arange(round(entry["T_y0"]) - 14.0, round(entry["T_y0"]) + 8.1, 2.0).tolist()]
        src = key36 if name == "one_doublet_B" else key24
        if src is not None and name in sets[src]["tscan"]:
            t = sets[src]["tscan"][name]
            T_lat, T_lat_err = t["T_max"], t["T_max_err"]
            red = rd.reduce(T_lat, **pars, linear=True)
            y_c = red["y"] + dy_cont
            T_c = rd.solve_T(lambda x: y_c, pars, linear=True)
            dydT = (rd.reduce(T_c + 0.5, **pars, linear=True)["y"] - rd.reduce(T_c - 0.5, **pars, linear=True)["y"])
            stat = math.hypot(T_lat_err, fit_y["err"] / dydT)
            entry.update(set_used=src, T_max_lattice=T_lat, T_max_lattice_err=T_lat_err, x_at_T_c=rd.reduce(T_c, **pars, linear=True)["x"],
                         y_at_T_max_lattice=red["y"], y_c_continuum=y_c, T_c=T_c, T_c_err_lattice=stat, dy_dT=dydT,
                         chi_max=t["chi_max"], chi_max_err=t["chi_max_err"], half_width_high_side_GeV=t["half_width_sym"],
                         g32_over_T=rd.reduce(T_c, **pars, linear=True)["g32"] / T_c)
            # variants of the map at fixed y_c
            var = {}
            for label, kw in (("product_form_mu_piT", dict(linear=False, mu_over_T=math.pi)),
                              ("product_form_mu_7T", dict(linear=False)),
                              ("no_hypercharge_in_map", dict(linear=True, gprime=False)),
                              ("tree_level_reduction", dict(tree=True))):
                try:
                    var[label] = rd.solve_T(lambda x: y_c, pars, **kw)
                except Exception:
                    var[label] = None
            p2 = dict(pars)
            p2["alpha_s"] = pars["alpha_s"] / (1.0 + pars["alpha_s"] * 23.0 / (12.0 * math.pi) * math.log(pars["mW"] ** 2 / pars["mZ"] ** 2))
            var["alpha_s_run_to_mW"] = rd.solve_T(lambda x: y_c, p2, linear=True)
            entry["variants_T_c"] = var
            lo80, hi80 = width_at_fraction(t["T"], t["chi"], 0.8)
            entry["interval_above_80_percent_of_maximum_GeV"] = [lo80, hi80]
            entry["width_at_80_percent_GeV"] = (hi80 - lo80) if lo80 is not None and hi80 is not None else None
            entry["truncation_error"] = 0.01 * T_c
            entry["T_c_err_total"] = math.hypot(stat, 0.01 * T_c)
            sec = rd.cosmic_time_seconds(T_c)
            entry["cosmic_time_seconds"] = sec
            entry["cosmic_time_ticks"] = rd.ticks(sec)
            entry["kind"] = "crossover"
        rows_out[name] = entry
    out["rows"] = rows_out
    dec = rows_out["declared"]
    if "T_c" in dec:
        base = ROWS["declared"]

        def shifted(**kw):
            r2 = dict(base); r2.update(kw)
            return rd.solve_T(lambda x: dec["y_c_continuum"], inputs_of(r2), linear=True) - dec["T_c"]
        dec["sensitivity_GeV"] = {
            "v_plus_1_percent": shifted(v=base["v"] * 1.01), "g2_plus_1_percent": shifted(g2=base["g2"] * 1.01),
            "sin2_plus_1_percent": shifted(sin2=base["sin2"] * 1.01), "m_H_plus_1_GeV": shifted(m_H=base["m_H"] + 1.0),
            "m_t_plus_1_GeV": shifted(m_t=base["m_t"] + 1.0), "alpha_s_plus_0.001": shifted(alpha_s=base["alpha_s"] + 0.001)}
    if "T_c" in dec:
        out["comparison_not_blind"] = dict(
            T_c_declared=dec["T_c"], err=dec["T_c_err_total"], reference=REFERENCE,
            difference=dec["T_c"] - REFERENCE["T_c"],
            pull=(dec["T_c"] - REFERENCE["T_c"]) / math.hypot(dec["T_c_err_total"], REFERENCE["err"]),
            pipeline_on_standard_model_inputs=rows_out["standard_model"].get("T_c"),
            pipeline_minus_published_continuum_peak=(rows_out["standard_model"].get("T_c", float("nan")) - REFERENCE["continuum_peak"]))
        checks["5_rows.declared_is_crossover"] = dec["x_at_T_c"] > X_C["value"] + 5 * X_C["err"]
        checks["5_rows.T_max_consistent_with_y_max"] = abs(
            rd.solve_T(lambda x: lat[key24]["y_max"], inputs_of(ROWS["declared"]), linear=True) - dec["T_max_lattice"]) < 0.15
        sm_row = rows_out["standard_model"]
        if "T_max_lattice" in sm_row:
            interp = REFERENCE["peak_positions"]["6"] + (REFERENCE["peak_positions"]["9"] - REFERENCE["peak_positions"]["6"]) * (1 / 6 - 1 / 8) / (1 / 6 - 1 / 9)
            out["comparison_not_blind"]["standard_model_peak_beta_G_8"] = dict(
                this_pipeline=sm_row["T_max_lattice"], err=sm_row["T_max_lattice_err"], published_interpolated=interp)
            checks["5_rows.standard_model_peak_matches_published"] = abs(sm_row["T_max_lattice"] - interp) < 0.5
            ts = sets[key24]["tscan"]["standard_model"]
            curve = []
            for T, pub in PUBLISHED_CURVE:
                if ts["T"][0] <= T <= ts["T"][-1]:
                    curve.append(dict(T=T, published=pub, this=float(np.interp(T, ts["T"], ts["chi"])),
                                      this_err=float(np.interp(T, ts["T"], ts["chi_err"]))))
            out["comparison_not_blind"]["standard_model_curve"] = curve
            # position of the maximum in y: published SU(2) x U(1) simulation against this SU(2) simulation
            y_pub = float(np.interp(interp, [t for t, _ in PUBLISHED_Y], [y for _, y in PUBLISHED_Y]))
            dy = y_pub - sm_row["y_at_T_max_lattice"]
            out["hypercharge_field_estimate"] = dict(
                y_max_published_with_U1=y_pub, y_max_this_without_U1=sm_row["y_at_T_max_lattice"], dy=dy,
                dT_GeV_declared=dy / dec["dy_dT"],
                note="difference of the peak position in y at beta_G = 8 between the published simulation with the "
                     "U(1) field (z = 0.31) and this simulation (z = 0); read from figures, accuracy about 0.003 in y")
            checks["5_rows.standard_model_curve_matches_published"] = bool(curve) and max(
                abs(c["this"] / c["published"] - 1.0) for c in curve) < 0.15

    # ---- 5. endpoint and the Higgs mass scan of the declared chain --------------------------------------------
    val = prod.get("validation", {})
    ycx = []                                                  # (x, y_max) at beta_G = 8 on the largest volume measured
    fo = val.get("first_order", {})
    for L in ("16", "12", "8"):
        if L in fo.get("sets", {}):
            ycx.append((fo["x"], fo["sets"][L]["y_max"]))
            break
    for p in val.get("endpoint", {}).get("points", []):
        for L in ("24", "16"):
            if L in p["sets"]:
                ycx.append((p["x"], p["sets"][L]["y_max"]))
                break
    co = val.get("crossover", {})
    for L in ("32", "24", "16"):
        if L in co.get("sets", {}):
            ycx.append((co["x"], co["sets"][L]["y_max"]))
            break
    if key24:
        ycx.append((0.296, lat[key24]["y_max"]))
    if key36:
        ycx.append((0.36, lat[key36]["y_max"]))
    ycx.sort()
    out["peak_position_vs_x_beta_G_8"] = [dict(x=a, y_max=b) for a, b in ycx]
    if len(ycx) >= 3 and "T_c" in dec:
        xs = np.array([a for a, _ in ycx]); yv = np.array([b for _, b in ycx])

        def y_of_x(x):
            return float(np.interp(x, xs, yv))
        base = ROWS["declared"]
        scan = []
        for m_h in np.arange(50.0, 156.0, 5.0).tolist() + [base["m_H"]]:
            pars = rd.chain_inputs(base["v"], base["g2"], base["sin2"], m_h, base["m_t"], base["alpha_s"])
            try:
                T = rd.solve_T(y_of_x, pars, linear=True, lo=30.0)
            except Exception:
                continue
            red = rd.reduce(T, **pars, linear=True)
            scan.append(dict(m_H=m_h, T_peak_beta_G_8=T, x=red["x"], y=red["y"],
                             kind="first-order" if red["x"] < X_C["value"] else "crossover"))
        scan.sort(key=lambda r_: r_["m_H"])
        out["higgs_mass_scan_declared_chain"] = dict(
            note="other inputs of the declared chain fixed; peak position at beta_G = 8 without continuum extrapolation",
            rows=scan)

        def x_at_transition(m_h):
            pars = rd.chain_inputs(base["v"], base["g2"], base["sin2"], m_h, base["m_t"], base["alpha_s"])
            try:
                T = rd.solve_T(y_of_x, pars, linear=True, lo=30.0)
            except Exception:
                return -1.0
            return rd.reduce(T, **pars, linear=True)["x"]

        def solve_mh(xc):
            lo, hi = 15.0, 150.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if x_at_transition(mid) < xc:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)
        mhc = solve_mh(X_C["value"])
        out["endpoint"] = dict(x_c=X_C, m_H_c_declared_chain=mhc,
                               m_H_c_err_from_x_c=0.5 * abs(solve_mh(X_C["value"] + X_C["err"]) - solve_mh(X_C["value"] - X_C["err"])),
                               m_H_at_x={f"{xv:.5f}": (solve_mh(xv) if 16.0 < solve_mh(xv) < 149.0 else None)
                                         for xv in (0.036, 0.06445, 0.0897, 0.115, 0.13)},
                               m_H_tree_at_x={f"{xv:.5f}": inputs_of(base)["mW"] * math.sqrt(8.0 * xv)
                                              for xv in (0.036, 0.06445, 0.0897, 0.115, 0.13)},
                               note="m_H_at_x: Higgs mass at which the declared chain reaches x on its transition line; "
                                    "None where the map has no such mass. m_H_tree_at_x: tree level, x = m_H^2 / (8 m_W^2).")
        sm_pars = dict(STANDARD_MODEL)

        def x_sm(m_h):
            p = dict(sm_pars); p["mH"] = m_h
            try:
                T = rd.solve_T(y_of_x, p, linear=True, lo=30.0)
            except Exception:
                return -1.0
            return rd.reduce(T, **p, linear=True)["x"]
        lo, hi = 45.0, 150.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if x_sm(mid) < X_C["value"]:
                lo = mid
            else:
                hi = mid
        out["endpoint"]["m_H_c_standard_model_inputs_this_map"] = 0.5 * (lo + hi)
        checks["6_endpoint.standard_model_value_reproduced"] = abs(0.5 * (lo + hi) - X_C["m_H_SM"]) < 2 * X_C["m_H_SM_err"] + 1.0
        checks["6_endpoint.declared_higgs_above_endpoint"] = base["m_H"] > mhc

    # ---- 6. validation summaries --------------------------------------------------------------------------------
    bp = val.get("broken_phase")
    if bp:
        fit = bp["fits"].get("linear_bG>=8") or next(iter(bp["fits"].values()))
        out["validation_broken_phase"] = dict(fit=fit, literature=bp["literature"], points=bp["points"])
        checks["7_validation.broken_phase_condensate"] = abs(fit["continuum"] - bp["literature"]["value"]) < 3 * math.hypot(
            fit["err"], bp["literature"]["err"])
    out["validation"] = {k: val[k] for k in ("first_order", "crossover", "endpoint") if k in val}
    for k, g in prod.get("gates", {}).items():
        checks["0_gates." + k] = bool(g.get("ok"))
    out["gates"] = prod.get("gates", {})
    out["checks"] = checks
    out["all_checks_pass"] = all(checks.values())
    out["input_sha256"] = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                           for name in ("production.json", "lattice.py", "reduction.py", "reweight.py", "kernels.py")
                           if (HERE / name).exists()}
    return json.loads(json.dumps(out))
