"""Test a universal finite edge against a resolved gas-rich dwarf and weak lensing."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import measure
from .. import lensing_profile
from ..udg_transfer.measure import GM_SUN, KPC


HERE = Path(__file__).resolve().parent


def edge_bound(mass_solar, radius_kpc, vmax_kms, low_scale, nominal_scale):
    gm = GM_SUN * mass_solar
    radius = radius_kpc * KPC
    # From Vdark² <= Vobs² and Vdark²=√(GM a0) √(gb/gt) outside the edge.
    ratio = (vmax_kms * 1000) ** 2 * radius / gm
    if (vmax_kms * 1000) ** 4 >= gm * low_scale:
        raise ValueError(
            "Dwarf observation does not require an outer edge at this calibration"
        )
    return dict(
        mass_solar=mass_solar,
        radius_kpc=radius_kpc,
        allowed_speed_kms=vmax_kms,
        max_halo_to_baryonic_total_mass=ratio,
        min_gt_m_s2=nominal_scale / ratio**2,
        min_gt_at_low_scale=low_scale / ratio**2,
    )


def profile_surface_bound(gbar, source_gbar, source_mass_ratio):
    """Maximum weak differential surface density for the common sublinear profile."""
    return (
        (4 / math.pi)
        * gbar
        * (1 + source_mass_ratio * np.maximum(1.0, np.sqrt(source_gbar / gbar)))
        / lensing_profile.FOUR_G
    )


def run():
    primary = json.loads(lensing_profile.RECEIPT.read_text())
    weights = json.loads((HERE / "results.json").read_text())["surveys"]
    rings = HERE.parent / "udg_transfer" / "rings.json"
    ring = json.loads(rings.read_text())["panels"]["figure8"]["rings"][-1]
    v = ring["circular_kms"]
    upper_v = v + 2 * ring["plotted_error_kms"]
    radius = ring["radius_kpc"]
    low = primary["sparc_cluster_bootstrap_95pct_m_s2"][0]
    a0 = primary["fitted_sparc_a0_m_s2"]
    cases = {}
    for label, m, r, vmax, lens_mass in (
        ("central", 1.47e9, radius, v, 1.0),
        ("dwarf_two_errors", 1.19e9, radius * (78.7 + 3) / 78.7, upper_v, 1.0),
        (
            "tilt_and_lens_mass",
            1.19e9,
            radius * (78.7 + 3) / 78.7,
            upper_v * math.sin(math.radians(31)) / math.sin(math.radians(27)),
            10**0.2,
        ),
    ):
        bound = edge_bound(m, r, vmax, low, a0)
        comparisons = {}
        for name, survey in weights.items():
            p = primary["profiles"][name]
            lower = np.asarray(p["lower_bin_edge_m_s2"])
            high = lens_mass * (2 * np.asarray(p["g_bar_m_s2"]) - lower)[lower >= 1e-13]
            if high.max() >= bound["min_gt_at_low_scale"]:
                raise ValueError(
                    "One or more source bins can see the interior form and require explicit projection"
                )
            max_esd = (
                (4 / math.pi)
                * high
                * (1 + bound["max_halo_to_baryonic_total_mass"])
                / lensing_profile.FOUR_G
            )
            upper = float(np.array(survey["positive_weights"]) @ max_esd)
            dwarf_gbar = GM_SUN * m / (r * KPC) ** 2
            broad_by_bin = profile_surface_bound(
                high, dwarf_gbar, bound["max_halo_to_baryonic_total_mass"]
            )
            broad_upper = float(np.asarray(survey["positive_weights"]) @ broad_by_bin)
            comparisons[name] = dict(
                sublinear_shared_profile=dict(
                    upper_by_bin=broad_by_bin.tolist(),
                    weighted_upper=broad_upper,
                    gap_over_statistical_se=(survey["weighted_signal"] - broad_upper)
                    / survey["statistical_se"],
                ),
                largest_allowed_gbar_m_s2=float(high.max()),
                upper_by_bin=max_esd.tolist(),
                weighted_upper=upper,
                observed=survey["weighted_signal"],
                statistical_se=survey["statistical_se"],
                cross_over_statistical_se=survey["cross_over_statistical_se"],
                gap_over_statistical_se=(survey["weighted_signal"] - upper)
                / survey["statistical_se"],
            )
        cases[label] = dict(**bound, lens_mass_factor=lens_mass, surveys=comparisons)
    return dict(
        dwarf_figure_data_sha256=hashlib.sha256(rings.read_bytes()).hexdigest(),
        lens_archive_sha256=primary["archive_sha256"],
        cases=cases,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    for label, row in result["cases"].items():
        print(
            label,
            row["min_gt_m_s2"],
            [
                (k, v["weighted_upper"], v["gap_over_statistical_se"])
                for k, v in row["surveys"].items()
            ],
        )
