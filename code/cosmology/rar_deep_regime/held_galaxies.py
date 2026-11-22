"""Cross-predict the total radial acceleration of held SPARC galaxies in the deep cut."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

import exterior_profile
import matched_observable_diagnostic as matched
import rar_deep_regime_diagnostic as base
import sparc_full_rar_calibration as empirical

RECEIPT = Path(__file__).with_name("receipts") / "held_galaxies.json"
FOLDS = 5


def fold_for(name):
    digest = hashlib.sha256(b"OPH-deep-held-v1\0" + name.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % FOLDS


def conditional_total(gbar, a0):
    return gbar + np.sqrt(gbar * a0)


def fit_scale(model, train_bar, train_truth, train_weight):
    def loss(log_a0):
        residual = train_truth - np.log10(model(train_bar, 10.0**log_a0))
        return np.average(residual**2, weights=train_weight)

    result = minimize_scalar(
        loss,
        bounds=base.LOG10_A0_BOUNDS,
        method="bounded",
        options={"xatol": 1e-12},
    )
    if (
        not result.success
        or min(abs(result.x - bound) for bound in base.LOG10_A0_BOUNDS) < 1e-4
    ):
        raise ValueError("Training optimizer failed or reached boundary")
    return 10.0**result.x


def cross_predict(acc):
    keep = acc["g_bar"] < matched.PRIMARY_FRACTION * base.DEEP_CUT_REFERENCE_A0_M_S2
    gbar, gobs, names = (acc[key][keep] for key in ("g_bar", "g_obs", "names"))
    galaxies, groups, counts = np.unique(names, return_inverse=True, return_counts=True)
    weights = 1.0 / counts[groups]
    folds = np.array([fold_for(name) for name in galaxies])[groups]
    truth = np.log10(gobs)
    models = {"conditional_deep": conditional_total, "empirical_rar": empirical.rar}
    predictions = {name: np.empty_like(truth) for name in models}
    predictions["baryons_only"] = np.log10(gbar)
    predictions["rescaled_baryons"] = np.empty_like(truth)
    fitted_scales = {name: [] for name in models}
    multipliers = []
    for fold in range(FOLDS):
        hold = folds == fold
        train = ~hold
        if not hold.any() or not train.any():
            raise ValueError("Empty fold")
        train_bar, train_truth, train_weight = gbar[train], truth[train], weights[train]
        log_multiplier = np.average(
            train_truth - np.log10(train_bar), weights=train_weight
        )
        multipliers.append(float(10.0**log_multiplier))
        predictions["rescaled_baryons"][hold] = np.log10(gbar[hold]) + log_multiplier
        for name, model in models.items():
            a0 = fit_scale(model, train_bar, train_truth, train_weight)
            fitted_scales[name].append(a0)
            predictions[name][hold] = np.log10(model(gbar[hold], a0))

    metrics, square_error_by_galaxy = {}, {}
    for name, pred in predictions.items():
        residual = truth - pred
        square_error_by_galaxy[name] = np.bincount(groups, weights=residual**2) / counts
        metrics[name] = {
            "rmse_dex": float(np.sqrt(np.mean(square_error_by_galaxy[name]))),
            "mean_observed_minus_predicted_dex": float(
                np.average(residual, weights=weights)
            ),
            "mae_dex": float(np.average(np.abs(residual), weights=weights)),
        }
    rng = np.random.default_rng(20261122)
    boot = rng.integers(len(galaxies), size=(4000, len(galaxies)))
    differences = {}
    for other in ("empirical_rar", "baryons_only", "rescaled_baryons"):
        diff = np.sqrt(
            square_error_by_galaxy["conditional_deep"][boot].mean(axis=1)
        ) - np.sqrt(square_error_by_galaxy[other][boot].mean(axis=1))
        differences[other] = {
            "conditional_minus_control_rmse_dex": metrics["conditional_deep"][
                "rmse_dex"
            ]
            - metrics[other]["rmse_dex"],
            "fixed_predictions_paired_galaxy_bootstrap_95pct": np.percentile(
                diff, [2.5, 97.5]
            ).tolist(),
        }
    return {
        "galaxies": len(galaxies),
        "radii": len(names),
        "nonpositive_observed_excess_retained": int(np.sum(gobs <= gbar)),
        "held_galaxies_per_fold": [
            len(np.unique(names[folds == f])) for f in range(FOLDS)
        ],
        "held_radii_per_fold": [int(np.sum(folds == f)) for f in range(FOLDS)],
        "training_scales_m_s2_in_fold_order": fitted_scales,
        "training_baryon_multipliers_in_fold_order": multipliers,
        "held_total_acceleration": metrics,
        "paired_galaxy_differences": differences,
    }


def same_ring_control(t1, t2):
    outer = matched.outermost_retained_points(t1, t2)
    channel = matched.channel_b(t1, t2, matched.PRIMARY_FRACTION)
    meta = {r["name"]: r for r in t1}
    names = channel["names"]
    gm = np.array(
        [
            base.G_SI
            * base.M_SUN
            * 1e9
            * (
                base.UPSILON_DISK * meta[n]["L36_gsun"]
                + base.GAS_HELIUM * meta[n]["MHI_gsun"]
            )
            for n in names
        ]
    )
    radius, vbar2, vobs = (
        np.array([outer[n][key] for n in names])
        for key in ("radius_m", "vbar2", "vobs_kms")
    )

    def estimates(speed):
        excess = (speed * base.KM) ** 2 - vbar2
        radial = (excess / radius) ** 2 / (vbar2 / radius)
        return np.log10(excess**2 / gm), np.log10(radial)

    btfr, radial = estimates(vobs)
    btfr_changed, radial_changed = estimates(1.2 * vobs)
    np.testing.assert_allclose(btfr, channel["log_a0"], atol=3e-14, rtol=0)
    geometry = np.log10(radius * vbar2 / gm)
    np.testing.assert_allclose(btfr - radial, geometry, atol=3e-14, rtol=0)
    np.testing.assert_allclose(
        btfr_changed - radial_changed, geometry, atol=3e-14, rtol=0
    )
    return {
        "outer_channel_b_galaxies": len(names),
        "control": "Multiply only outer Vobs by 1.2 after selection, holding baryons, radii and membership fixed.",
        "max_original_ratio_minus_baryon_geometry_abs_dex": float(
            np.max(np.abs(btfr - radial - geometry))
        ),
        "max_perturbed_minus_original_ratio_abs_dex": float(
            np.max(np.abs((btfr_changed - radial_changed) - (btfr - radial)))
        ),
        "median_individual_btfr_log_scale_increase_dex": float(
            np.median(btfr_changed - btfr)
        ),
        "mean_log10_geometry_ratio": float(np.mean(geometry)),
        "scope": "Single-ring ratio only. The published channel A fit also uses interior deep radii.",
    }


def training_sample(acc, excluded):
    keep = ~np.isin(acc["names"], excluded) & (
        acc["g_bar"] < matched.PRIMARY_FRACTION * base.DEEP_CUT_REFERENCE_A0_M_S2
    )
    galaxies, membership, counts = np.unique(
        acc["names"][keep], return_inverse=True, return_counts=True
    )
    train_bar, train_truth = acc["g_bar"][keep], np.log10(acc["g_obs"][keep])
    weights = 1.0 / counts[membership]
    return galaxies, membership, train_bar, train_truth, weights


def run():
    source = empirical.verify_source_snapshot()
    t1, t2 = base.read_table1(), base.read_table2()
    acc = base.accelerations(t1, t2)
    groups = exterior_profile.exterior_rings(t1, t2)
    held = [galaxy["name"] for galaxy, _ in groups]
    galaxies, membership, train_bar, train_truth, weights = training_sample(acc, held)
    scale = fit_scale(conditional_total, train_bar, train_truth, weights)
    rng = np.random.default_rng(20261123)
    galaxy_draws = rng.multinomial(
        len(galaxies), np.full(len(galaxies), 1.0 / len(galaxies)), size=2000
    )
    scales = [
        fit_scale(conditional_total, train_bar, train_truth, draw[membership] * weights)
        for draw in galaxy_draws
    ]
    exterior = exterior_profile.measure(groups, t1, t2)
    exterior["amplitude_transfer_from_disjoint_galaxies"] = {
        "training": "Same primary deep radial cut and equal-galaxy log-acceleration loss as the five-fold analysis. Exclude every exterior test galaxy. Transfer its single fitted a0 from gobs = gbar + sqrt(gbar a0) to the compact photometric exterior profile.",
        "training_galaxies": len(galaxies),
        "training_radii": len(train_bar),
        "fitted_a0_m_s2": scale,
        "training_galaxy_bootstrap": "2000 galaxy-cluster resamples with replacement, refitting each calibration, seed 20261123. Conditions on the catalogue inputs.",
        "training_galaxy_bootstrap_a0_95pct_m_s2": np.percentile(
            scales, [2.5, 97.5]
        ).tolist(),
        "minimum_maximum_velocity_error_lower_bound_by_geometry_box": {
            str(width): exterior_profile.minimum_error(
                lambda k: exterior_profile.fixed_scale_gap(groups, width, k, scale)
            )
            for width in (0, 1, 2)
        },
    }
    exterior["paired_baryon_source_transfer"] = exterior_profile.compare_baryon_sources(
        groups,
        scale,
        np.percentile(scales, [2.5, 97.5]),
        exterior["quoted_velocity_error_threshold"],
    )
    exterior["paired_baryon_source_transfer"]["limiting_ring_primary_source_check"] = (
        exterior_profile.ngc0024_projection_control(
            groups,
            scale,
            np.percentile(scales, 97.5),
            exterior["quoted_velocity_error_threshold"],
        )
    )
    exterior["paired_baryon_source_transfer"]["external_geometry_sensitivity"] = (
        exterior_profile.alternate_geometry(
            groups, scale, exterior["quoted_velocity_error_threshold"]
        )
    )
    # Probe two training-selection differences without fitting any held speed.
    bulges = sorted({row["name"] for row in t2 if row["vbul_kms"] != 0})
    training_variants = []
    for error_cut in (True, False):
        axes = acc if error_cut else base.accelerations(t1, t2, max_rel_error=np.inf)
        for bulge_free in (False, True):
            members, _, bars, truth, weight = training_sample(
                axes, held + (bulges if bulge_free else [])
            )
            value = fit_scale(conditional_total, bars, truth, weight)
            training_variants.append(
                {
                    "relative_velocity_error_cut": error_cut,
                    "bulge_free_only": bulge_free,
                    "training_galaxies": len(members),
                    "training_rings": len(bars),
                    "fitted_scale_m_s2": value,
                }
            )
    exterior["paired_baryon_source_transfer"]["training_selection_sensitivity"] = (
        training_variants
    )
    return {
        "source_file_receipts_sha256": source["file_receipts_sha256"],
        "protocol": {
            "split": "First eight bytes of SHA256(UTF8('OPH-deep-held-v1') + NUL + UTF8(galaxy)), big-endian, modulo five. Each galaxy is held once.",
            "weighting": "Equal galaxy weights in training and scoring; equal weights per retained radius within a galaxy.",
            "selection": "Quality <=2, inclination >=30 deg, positive speeds and baryon acceleration, gbar < 1.2e-11 m/s^2; primary relative velocity error <=0.10. Disk/bulge M/L 0.5/0.7.",
            "interval": "4000 paired galaxy bootstrap draws, seed 20261122, conditional on the five training fits. Excludes refitting and shared catalogue systematics.",
            "shared_inputs": "Observed and baryon accelerations use the same radius; the galaxy distance, orientation and baryonic mass model are catalogue inputs.",
            "source_identification": "SPARC records rotation speeds and baryonic components. No independent measurement of the modular collar charge or conversion from that charge to acceleration is present.",
            "scope": "Retrospective within-catalogue prediction. The empirical RAR was originally inferred using SPARC; holding out galaxies here tests transfer of fitted scale, not discovery on new observations.",
        },
        "primary": cross_predict(acc),
        "without_relative_velocity_error_cut": cross_predict(
            base.accelerations(t1, t2, max_rel_error=np.inf)
        ),
        "shared_ring_inputs": same_ring_control(t1, t2),
        "compact_exterior_profile": exterior,
    }


def verify(expected, measured, key=""):
    if isinstance(expected, dict):
        if measured.keys() != expected.keys():
            raise ValueError(f"Receipt keys changed at {key}")
        for field in expected:
            verify(expected[field], measured[field], field)
    elif isinstance(expected, list):
        for a, b in zip(expected, measured, strict=True):
            verify(a, b, key)
    elif isinstance(expected, float):
        # Near-zero cancellation residuals in dex depend on the arithmetic library.
        absolute_tolerance = 3e-14 if key.startswith("max_") else 1e-20
        np.testing.assert_allclose(
            measured, expected, rtol=1e-6, atol=absolute_tolerance
        )
    elif measured != expected:
        raise ValueError(f"Receipt changed at {key}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="write the archived receipt"
    )
    args = parser.parse_args()
    observed = run()

    if args.write:
        RECEIPT.write_text(json.dumps(observed, indent=2, allow_nan=False) + "\n")
    else:
        verify(json.loads(RECEIPT.read_text()), observed)
    print(json.dumps(observed["primary"]["held_total_acceleration"], indent=2))
