"""Pass a cosmic-age history through the retained DESI DR1 Figure 5 slope windows."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import erfc, sqrt
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from level202_chronometers import (
    CORR,
    H,
    MPC_KM,
    STAT,
    SYST,
    Z,
    background,
    screen_report,
)

PDF_SHA = "3ac2b5c5a5b2051108e0a925a45e99bcd5dbe8f0259ac80e2aadd7252de61a60"
BREAKS = (0.55, 0.78)
# Cyan retained points in Loubser (2025), Fig. 5. Bin i has z = 0.305 + i/100.
CLIPPED = [2, 3, 8, 14, 15, 16, 18, 21, 22, 24, 31, 32, 33, 36, 37, 42, 55]
BINS = np.delete(np.arange(70), CLIPPED)
SUB_Z = 0.305 + BINS / 100
CODE_ROOT = Path(__file__).resolve().parents[1]


def slopes(values):
    design = np.column_stack(
        [np.ones_like(SUB_Z), SUB_Z] + [np.maximum(SUB_Z - k, 0) for k in BREAKS]
    )
    coeff = np.linalg.lstsq(design, values, rcond=None)[0]
    return np.cumsum(coeff[1:])


def check_figure(path):
    """Read the vertical cyan error-bar centres and check the printed slopes."""
    import pymupdf

    if sha256(path.read_bytes()).hexdigest() != PDF_SHA:
        raise ValueError("Expected arXiv:2511.02730v1")
    with pymupdf.open(path) as doc:
        page = doc[5]
        # In PDF points, these are the z=.30 and z=.55 vertical guide positions.
        left, first_break = 139.0021209716797, 258.8000793457031
        # The 1.4–2.6 y tick labels differ by 180.475 PDF points.
        scale_y = 1.2 / 180.475
        centres = {}
        for shape in page.get_drawings():
            box = shape["rect"]
            if shape["color"] != (0, 1, 1) or shape["fill"] is not None or box.width:
                continue
            z = 0.30 + 0.25 * (box.x0 - left) / (first_break - left)
            i = round(100 * (z - 0.305))
            if 0 <= i < 70 and abs(z - (0.305 + i / 100)) < 1e-5:
                if i in centres:
                    raise ValueError("Duplicate plotted bin")
                centres[i] = -scale_y * (box.y0 + box.y1) / 2
    if sorted(centres) != BINS.tolist():
        raise ValueError("Retained points disagree with Figure 5")
    measured = slopes([centres[i] for i in BINS])
    if not np.array_equal(np.round(measured, 3), [-0.462, -0.299, -0.301]):
        raise ValueError("Digitized slopes disagree with the printed Figure 5 fit")
    return measured


def reported_rates(vacuum, matter):
    # Cosmic age up to an additive constant. This form admits the pure-vacuum
    # limit and avoids subtracting large ages when the matter density is small.
    u = 1 + SUB_Z
    ages = (
        2 * np.log(vacuum + np.sqrt(vacuum**2 + matter * u**3)) / 3 - np.log(u)
    ) / vacuum
    return -1 / ((1 + Z) * slopes(ages))


def observer_zero_pulse_readout():
    """Count measured periods, including sources without a paired radar echo."""
    source = Path("cosmology/screen_coherence/level202_distance_redshift.json")
    raw = (CODE_ROOT / source).read_bytes()
    survey = json.loads(raw)
    panels = []
    for tick in survey["query_epochs"]:
        rows = survey["panels"][str(tick)]
        shifts = Counter()
        with_echo = 0
        for row in rows:
            arrivals = row["beacon_received"]
            if arrivals is None or None in arrivals:
                continue
            emit_a, emit_b = row["beacon_sent"]
            shift = Fraction(arrivals[1] - arrivals[0], emit_b - emit_a) - 1
            shifts[shift] += 1
            with_echo += row["echo_at"] is not None
        panels.append(
            dict(
                query_tick=tick,
                pulse_shift_counts={str(z): n for z, n in sorted(shifts.items())},
                missing_pulse_pairs=len(rows) - sum(shifts.values()),
                complete_pulse_and_echo=with_echo,
                pulse_pairs_in_galaxy_fit_range=sum(
                    n for z, n in shifts.items() if Fraction(3, 10) < z < 1
                ),
            )
        )
    return dict(
        source=f"code/{source.as_posix()}",
        source_sha256=sha256(raw).hexdigest(),
        receiver=survey["receiver"],
        enrolled_beacons_per_epoch=survey["enrolled"],
        panels=panels,
        chronometer_H_km_s_Mpc=None,
        instrument="The material beacons supply pulse periods and radar echoes. A galaxy-clock rate requires spectroscopic redshifts and independently calibrated ages of comoving tracers.",
    )


def measure():
    screen = screen_report()["level202"]
    h = float(screen["H_dS_s_minus_1"]) * MPC_KM
    sigma = np.hypot(STAT, SYST)
    cov = np.outer(sigma, sigma) * CORR
    chol = np.linalg.cholesky(cov)

    def calibrate(anchor):
        return brentq(lambda b: reported_rates(h, b)[0] - anchor, 0, 1e5)

    mat = calibrate(H[0])
    predictions = reported_rates(h, mat)
    # Implicit derivatives dH_i/dH_anchor from the smooth slope estimator.
    delta = mat * 1e-5
    grad = (reported_rates(h, mat + delta) - reported_rates(h, mat - delta)) / (
        2 * delta
    )
    jac = np.column_stack([-grad[1:] / grad[0], np.eye(2)])
    uncertainty = np.sqrt(np.diag(jac @ cov @ jac.T))
    diff = H[1:] - predictions[1:]

    def joint_residual(b):
        return np.linalg.solve(chol, H - reported_rates(h, b[0]))

    fit = least_squares(joint_residual, [mat], bounds=(0, np.inf))
    if not fit.success:
        raise ValueError(fit.message)
    archival = (
        CODE_ROOT
        / "cosmology/fixed_capacity_wlaw/runtime/official_desi_dr2_fz13_retrospective.json"
    )
    chain = json.loads(archival.read_text())["base_lcdm_capacity_display"]
    reference = chain["combined"]["Lambda_lP2"]
    nominal = float(screen["Lambda_lP2"])
    return dict(
        source="https://arxiv.org/pdf/2511.02730v1",
        source_sha256=PDF_SHA,
        slopes_printed_in_figure_5=[-0.462, -0.299, -0.301],
        excluded_bin_indices=CLIPPED,
        assumptions="Same clipped bins and unweighted continuous piecewise fit as reconstructed from the published figure. Linear spectral index versus stellar age with a common formation-time offset. The published covariance, median redshifts and breakpoints are held fixed. Galaxy selection and clipping are not resimulated.",
        forecast_inputs="Assigned de Sitter screen area, flat Friedmann evolution with dust and constant Lambda, and the first astronomical bin for the matter calibration. The separate pulse readout uses the simulated emitter and receiver records.",
        background=background(h, mat),
        reported_H_for_three_bins=predictions.tolist(),
        clock_anchor_variation=[
            background(h, calibrate(H[0] + s * sigma[0])) for s in (-1, 1)
        ],
        measured_minus_predicted=diff.tolist(),
        linearized_residual_sigma=uncertainty.tolist(),
        marginal_standardized_residuals=(diff / uncertainty).tolist(),
        gaussian_two_sided_marginal_tails=[
            erfc(abs(x) / sqrt(2)) for x in diff / uncertainty
        ],
        all_bin_chi2=2 * fit.cost,
        approximate_tail_two_dof=np.exp(-fit.cost),
        official_flat_BAO_CMB=dict(
            source=chain["source_directory"],
            Lambda_lP2_mean=reference["weighted_mean"],
            Lambda_lP2_sd=reference["weighted_std"],
            level202=nominal,
            level202_minus_mean_in_posterior_sd=(nominal - reference["weighted_mean"])
            / reference["weighted_std"],
            interpretation="The fixed area assigned to the screen is identified with a de Sitter horizon. The distance to the flat-LambdaCDM posterior uses CMB data shared with level selection. Determining a vacuum stress from observer records requires an independent gravitational measurement.",
        ),
        observer_zero_pulse_readout=observer_zero_pulse_readout(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-figure", type=Path)
    args = parser.parse_args()
    if args.verify_figure:
        check_figure(args.verify_figure)
    print(json.dumps(measure(), indent=2, allow_nan=False))
