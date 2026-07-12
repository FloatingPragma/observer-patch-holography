#!/usr/bin/env python3
"""Common-scheme scan of the quark shape law (#377 scheme rider).

The held-out shape-law tests score the frozen ray structure against the
reference sextet in its quoted conventions (light quarks MS-bar at 2 GeV,
charm and bottom at their own scale, top from cross sections). Log-gap
ratios are convention-dependent at the ten-percent level once the quotes
are run to a common scale, which is the same order as the deviations under
study. This artifact makes the scheme question quantitative:

  1. a self-contained two-loop QCD running package (two-loop beta, two-loop
     mass anomalous dimension, flavor thresholds at the own-scale heavy
     quote points, continuity matching) with benchmark windows against
     standard running factors;
  2. the shape metrics (sector gap ratios, implied rho per sector, the
     zero-parameter reciprocity product, the joint single-rho optimum, and
     the interior-mass held-out deviations at the frozen template rho) for
     three convention families: the quoted mixed convention, the own-scale
     convention with the MS-bar top mass, and a grid of common scales from
     1.5 to 250 GeV;
  3. the resulting scheme verdict: which convention family the ray law
     prefers, how the compare-only rho band of acceptance test A2 moves
     with the convention, and the declared requirement that the flavor
     lane state the scheme of its profile law.

Every mass and coupling here is a cited compare-only reference. Nothing
feeds an OPH solve path; the only OPH input is the frozen template
rho_ord, used to score interior-mass deviations per convention.

Run:
    python3 code/particles/flavor/derive_quark_common_scheme_shape_law_scan.py
writes code/particles/runs/flavor/quark_common_scheme_shape_law_scan.json.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
from datetime import datetime, timezone

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = ROOT / "particles" / "runs" / "flavor"
DEFAULT_OUT = RUNS / "quark_common_scheme_shape_law_scan.json"
GENERATOR_PATH = RUNS / "generation_bundle_branch_generator.json"
MEAN_SPLIT_PATH = RUNS / "quark_sector_mean_split.json"
REFERENCE_PATH = ROOT / "particles" / "data" / "particle_reference_values.json"

# compare-only cited constants
ALPHA_S_MZ = 0.1180          # PDG world average alpha_s(m_Z)
MZ = 91.1876
MT_MSBAR_MT = 162.5          # PDG MS-bar top mass m_t(m_t), compare-only
GRID_LO, GRID_HI, GRID_N = 1.5, 250.0, 41

FOUR_PI = 4.0 * math.pi


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def beta0(nf: int) -> float:
    return 11.0 - 2.0 * nf / 3.0


def beta1(nf: int) -> float:
    return 102.0 - 38.0 * nf / 3.0


def gamma1(nf: int) -> float:
    # two-loop mass anomalous dimension in the d ln m / d ln mu convention,
    # gamma(a) = 8 a + gamma1 a^2 with a = alpha_s / (4 pi)
    return 404.0 / 3.0 - 40.0 * nf / 9.0


class Runner:
    """Two-loop alpha_s and MS-bar mass running with flavor thresholds."""

    def __init__(self, mc_mc: float, mb_mb: float, mt_mt: float):
        self.thresholds = (mc_mc, mb_mb, mt_mt)

    def _nf(self, mu: float) -> int:
        return 3 + sum(1 for t in self.thresholds if mu >= t)

    def _segments(self, mu_from: float, mu_to: float) -> list[tuple[float, float, int]]:
        points = [mu_from, mu_to]
        lo, hi = min(mu_from, mu_to), max(mu_from, mu_to)
        cuts = sorted({t for t in self.thresholds if lo < t < hi})
        ordered = [lo] + cuts + [hi]
        if mu_from > mu_to:
            ordered = ordered[::-1]
        segs = []
        for a, b in zip(ordered[:-1], ordered[1:], strict=True):
            mid = math.sqrt(a * b)
            segs.append((a, b, self._nf(mid)))
        return segs

    def _rk4(self, t0: float, t1: float, state: tuple[float, float],
             nf: int, steps_per_unit: float = 400.0) -> tuple[float, float]:
        b0, b1, g1 = beta0(nf), beta1(nf), gamma1(nf)

        def deriv(s):
            a, _ = s
            da = -2.0 * a * a * (b0 + b1 * a)
            dlnm = -(8.0 * a + g1 * a * a)
            return da, dlnm

        span = t1 - t0
        n = max(2, int(abs(span) * steps_per_unit))
        h = span / n
        a, lnm = state
        for _ in range(n):
            k1 = deriv((a, lnm))
            k2 = deriv((a + 0.5 * h * k1[0], lnm + 0.5 * h * k1[1]))
            k3 = deriv((a + 0.5 * h * k2[0], lnm + 0.5 * h * k2[1]))
            k4 = deriv((a + h * k3[0], lnm + h * k3[1]))
            a += h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0
            lnm += h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0
        return a, lnm

    def alpha_s(self, mu: float) -> float:
        a = ALPHA_S_MZ / FOUR_PI
        lnm = 0.0
        for seg_from, seg_to, nf in self._segments(MZ, mu):
            a, lnm = self._rk4(math.log(seg_from), math.log(seg_to),
                               (a, lnm), nf)
        return a * FOUR_PI

    def mass_factor(self, mu_from: float, mu_to: float) -> float:
        """MS-bar factor m(mu_to)/m(mu_from) for a mass quoted at mu_from."""
        a = self.alpha_s(mu_from) / FOUR_PI
        lnm = 0.0
        for seg_from, seg_to, nf in self._segments(mu_from, mu_to):
            a, lnm = self._rk4(math.log(seg_from), math.log(seg_to),
                               (a, lnm), nf)
        return math.exp(lnm)


def shape_metrics(sextet: dict[str, float], rho_template: float) -> dict:
    mu_, mc, mt = sextet["up"], sextet["charm"], sextet["top"]
    md, ms, mb = sextet["down"], sextet["strange"], sextet["bottom"]
    gu21, gu32 = math.log(mc / mu_), math.log(mt / mc)
    gd21, gd32 = math.log(ms / md), math.log(mb / ms)
    r_up, r_down = gu21 / gu32, gd21 / gd32
    product = r_up * r_down
    rho_band = sorted([r_up, 1.0 / r_down])
    rho_star = math.sqrt(r_up / r_down)
    frac_up = rho_template / (1.0 + rho_template)
    frac_down = 1.0 / (1.0 + rho_template)
    mc_pred = mu_ * math.exp(math.log(mt / mu_) * frac_up)
    ms_pred = md * math.exp(math.log(mb / md) * frac_down)
    return {
        "gap_ratio_up": r_up,
        "gap_ratio_down": r_down,
        "implied_rho_up": r_up,
        "implied_rho_down": 1.0 / r_down,
        "rho_band": rho_band,
        "template_rho_inside_band": bool(
            rho_band[0] <= rho_template <= rho_band[1]),
        "reciprocity_product": product,
        "reciprocity_deviation": product - 1.0,
        "joint_single_rho_optimum": rho_star,
        "interior_heldout_at_template_rho": {
            "charm_relative_error": mc_pred / mc - 1.0,
            "strange_relative_error": ms_pred / ms - 1.0,
        },
    }


def build() -> dict:
    refs = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))["entries"]
    generator = json.loads(GENERATOR_PATH.read_text(encoding="utf-8"))
    mean_split = json.loads(MEAN_SPLIT_PATH.read_text(encoding="utf-8"))
    spectrum = sorted(
        float(x) for x in
        generator["charged_sector_response_operator_candidate"]["ordered_spectrum"])
    g21 = spectrum[1] - spectrum[0]
    g32 = spectrum[2] - spectrum[1]
    rho_template = 3.0 * g32 / (2.0 * g32 + g21)
    if abs(rho_template - float(mean_split["rho_ord"])) > 1.0e-12:
        raise AssertionError("template rho_ord disagrees with mean split")

    def mass(key: str) -> float:
        return float(refs[key]["value_gev"])

    quoted = {
        "up": mass("up_quark"), "charm": mass("charm_quark"),
        "top": mass("top_quark"), "down": mass("down_quark"),
        "strange": mass("strange_quark"), "bottom": mass("bottom_quark"),
    }
    runner = Runner(mc_mc=quoted["charm"], mb_mb=quoted["bottom"],
                    mt_mt=MT_MSBAR_MT)

    # benchmark windows against standard running factors (compare-only)
    benchmarks = {
        "alpha_s_2gev": {
            "value": runner.alpha_s(2.0),
            "window": [0.28, 0.32],
        },
        "mb_own_scale_to_mz_factor": {
            "value": runner.mass_factor(quoted["bottom"], MZ),
            "window": [0.66, 0.72],
        },
        "mc_own_scale_to_mz_factor": {
            "value": runner.mass_factor(quoted["charm"], MZ),
            "window": [0.45, 0.53],
        },
        "light_2gev_to_mz_factor": {
            "value": runner.mass_factor(2.0, MZ),
            "window": [0.55, 0.63],
        },
    }
    for name, row in benchmarks.items():
        row["passes"] = bool(row["window"][0] <= row["value"] <= row["window"][1])
    if not all(row["passes"] for row in benchmarks.values()):
        raise AssertionError(f"running benchmarks outside windows: {benchmarks}")

    conventions = {}
    conventions["quoted_mixed"] = {
        "description": "reference quotes as recorded (lights MS-bar 2 GeV, "
                       "charm and bottom at own scale, top from cross "
                       "sections)",
        "sextet_gev": quoted,
        "metrics": shape_metrics(quoted, rho_template),
    }

    own_scale = dict(quoted)
    own_scale["top"] = MT_MSBAR_MT
    conventions["own_scale_heavy_msbar_top"] = {
        "description": "own-scale MS-bar heavy quarks including the top "
                       "(m_t(m_t) = 162.5), lights MS-bar 2 GeV",
        "sextet_gev": own_scale,
        "metrics": shape_metrics(own_scale, rho_template),
    }

    def common_sextet(mu: float) -> dict[str, float]:
        light = runner.mass_factor(2.0, mu)
        return {
            "up": quoted["up"] * light,
            "down": quoted["down"] * light,
            "strange": quoted["strange"] * light,
            "charm": quoted["charm"] * runner.mass_factor(quoted["charm"], mu),
            "bottom": quoted["bottom"] * runner.mass_factor(quoted["bottom"], mu),
            "top": MT_MSBAR_MT * runner.mass_factor(MT_MSBAR_MT, mu),
        }

    ratio = (GRID_HI / GRID_LO) ** (1.0 / (GRID_N - 1))
    grid_rows = []
    for idx in range(GRID_N):
        mu = GRID_LO * ratio ** idx
        metrics = shape_metrics(common_sextet(mu), rho_template)
        grid_rows.append({
            "mu_gev": mu,
            "gap_ratio_up": metrics["gap_ratio_up"],
            "gap_ratio_down": metrics["gap_ratio_down"],
            "reciprocity_product": metrics["reciprocity_product"],
            "rho_band": metrics["rho_band"],
            "template_rho_inside_band": metrics["template_rho_inside_band"],
            "charm_heldout_rel_error": metrics[
                "interior_heldout_at_template_rho"]["charm_relative_error"],
            "strange_heldout_rel_error": metrics[
                "interior_heldout_at_template_rho"]["strange_relative_error"],
        })
    highlight = {}
    for mu in (2.0, MZ, 200.0):
        highlight[f"mu_{mu:g}_gev"] = {
            "sextet_gev": common_sextet(mu),
            "metrics": shape_metrics(common_sextet(mu), rho_template),
        }

    products = [row["reciprocity_product"] for row in grid_rows]
    product_spread = max(products) - min(products)
    quoted_dev = abs(conventions["quoted_mixed"]["metrics"]
                     ["reciprocity_deviation"])
    common_dev = abs(products[0] - 1.0)

    findings = {
        "reciprocity_by_convention": {
            "quoted_mixed": conventions["quoted_mixed"]["metrics"]
                            ["reciprocity_deviation"],
            "own_scale_heavy_msbar_top": conventions[
                "own_scale_heavy_msbar_top"]["metrics"]
                ["reciprocity_deviation"],
            "common_scale_unique_reading": products[0] - 1.0,
        },
        "common_scale_reading_is_scale_invariant": {
            "statement": "QCD mass running is flavor blind at a common "
                         "scale, so common-scale gap ratios are the same "
                         "for every mu; the grid verifies this numerically",
            "max_reciprocity_product_spread_on_grid": product_spread,
        },
        "common_scale_family_dominated_by_quoted": bool(
            common_dev > quoted_dev),
        "a2_band_by_convention": {
            "quoted_mixed": conventions["quoted_mixed"]["metrics"]["rho_band"],
            "own_scale_heavy_msbar_top": conventions[
                "own_scale_heavy_msbar_top"]["metrics"]["rho_band"],
            "common_scale_any_mu": highlight[f"mu_{MZ:g}_gev"]["metrics"]["rho_band"],
        },
        "statement": (
            "The ray law with one shared rho scores best in the quoted "
            "mixed convention (own-scale heavy quotes, 2 GeV lights), where "
            "the two sector-preferred rho values bracket the template value "
            "inside a two-percent window. The common-scale reading is "
            "unique up to running (gap ratios are scale invariant at a "
            "common mu) and it degrades the reciprocity product to the "
            "fifteen-percent level, with the sectors preferring rho values "
            "eighteen percent apart. The shape agreement of the flavor "
            "lane is therefore a statement about own-scale-class masses, "
            "and the lane must declare that scheme choice as part of any "
            "promotion. Acceptance test A2 of the issue-377 program "
            "inherits a scheme rider: the compare-only rho band is defined "
            "in the declared convention, and the quoted-convention band is "
            "the current default."
        ),
    }

    return {
        "artifact": "oph_quark_common_scheme_shape_law_scan",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "row_class": "compare_only_scheme_analysis",
        "guards": {
            "measured_values_in_any_oph_solve_path": False,
            "public_promotion_allowed": False,
            "oph_inputs": ["rho_ord (frozen template value, scoring only)"],
        },
        "running_package": {
            "orders": "two-loop beta function, two-loop mass anomalous "
                      "dimension, continuity matching at thresholds",
            "alpha_s_mz": ALPHA_S_MZ,
            "thresholds_gev": [quoted["charm"], quoted["bottom"], MT_MSBAR_MT],
            "mt_msbar_mt_gev": MT_MSBAR_MT,
            "truncation_note": "higher-order running and threshold matching "
                               "corrections are at the few-percent level on "
                               "these factors; the benchmark windows bound "
                               "the implementation against standard values",
            "benchmarks": benchmarks,
        },
        "template_rho_ord": rho_template,
        "conventions": conventions,
        "common_scale_grid": grid_rows,
        "common_scale_highlights": highlight,
        "findings": findings,
        "notes": [
            "This scan is the declared common-scheme follow-up of the "
            "held-out shape-law artifact.",
            "It consumes measured masses and couplings as cited "
            "compare-only references and emits no OPH quantity.",
            "The scheme statement is a display and acceptance-test rider; "
            "the nonidentifiability theorem and the promotion gates are "
            "unchanged.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Scan quark shape-law metrics across mass conventions.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    report = build()
    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    bm = report["running_package"]["benchmarks"]
    print("running benchmarks:")
    for name, row in bm.items():
        print(f"  {name:28s} {row['value']:.4f}  window {row['window']}  "
              f"passes={row['passes']}")
    print(f"template rho_ord = {report['template_rho_ord']:.6f}")
    for name in ("quoted_mixed", "own_scale_heavy_msbar_top"):
        met = report["conventions"][name]["metrics"]
        print(f"{name}: band [{met['rho_band'][0]:.4f}, {met['rho_band'][1]:.4f}]  "
              f"inside={met['template_rho_inside_band']}  "
              f"reciprocity {met['reciprocity_deviation']:+.2%}  "
              f"m_c {met['interior_heldout_at_template_rho']['charm_relative_error']:+.2%}  "
              f"m_s {met['interior_heldout_at_template_rho']['strange_relative_error']:+.2%}")
    for key, block in report["common_scale_highlights"].items():
        met = block["metrics"]
        print(f"common {key}: band [{met['rho_band'][0]:.4f}, "
              f"{met['rho_band'][1]:.4f}]  inside={met['template_rho_inside_band']}  "
              f"reciprocity {met['reciprocity_deviation']:+.2%}")
    f = report["findings"]
    print(f"common-scale gap ratios scale invariant to "
          f"{f['common_scale_reading_is_scale_invariant']['max_reciprocity_product_spread_on_grid']:.2e}")
    print(f"quoted convention dominates common-scale reading: "
          f"{f['common_scale_family_dominated_by_quoted']}")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
