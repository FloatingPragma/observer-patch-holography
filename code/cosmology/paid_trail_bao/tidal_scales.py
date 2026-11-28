"""Test the exterior tide on galaxy disks at every positive physical hop length."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np

from . import moving_budget as budget
from . import tidal_scale_kernel as kernel

HERE = budget.HERE


def source_data():
    lab = json.loads((budget.LAB / "exterior_tide.json").read_bytes())
    energy = json.loads((budget.LAB / "moving_records.json").read_bytes())["currents"][
        "sharp_ordinary_standards"
    ]
    mu = Q(energy["1"]["energy"]) / Q(energy["2"]["energy"])
    profiles = {
        k: kernel.shape_data(lab, k, ref)
        for k, ref in [("source", "single"), ("positive", "reference_plus_single")]
    }
    return profiles, mu


def pixel_maps():
    maps = json.loads((HERE / "closed_inventory_maps.json").read_bytes())
    return {
        f"{tag}/{cap}": panel["closed"]
        for tag, menu in maps.items()
        for cap, panel in menu.items()
    }


def physical_bounds(raw, factor, unit):
    # Bound pi and round both the integer and the positive conversion factor
    # outward before multiplying. A negative entry reverses the factor bounds.
    pi_lo = Q("3.14159265358979323846264338327950288")
    pi_hi = pi_lo + Q(1, 10**35)
    gain = factor * unit / kernel.WEIGHT
    assert gain > 0
    low_gain = np.nextafter(float(gain / pi_hi), -np.inf)
    high_gain = np.nextafter(float(gain / pi_lo), np.inf)
    xlo = np.nextafter(np.asarray(raw, float), -np.inf)
    xhi = np.nextafter(np.asarray(raw, float), np.inf)
    return (
        np.nextafter(xlo * np.where(xlo < 0, high_gain, low_gain), -np.inf),
        np.nextafter(xhi * np.where(xhi < 0, low_gain, high_gain), np.inf),
    )


def target(tab):
    radii, v, e, b = [
        [Q(str(x)) for x in tab[k]]
        for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
    ]
    low = [
        (max(0, u - 2 * s) ** 2 - z * z) / r
        for r, u, s, z in zip(radii, v, e, b, strict=True)
    ]
    high = [
        ((u + 2 * s) ** 2 - z * z) / r
        for r, u, s, z in zip(radii, v, e, b, strict=True)
    ]
    return radii, low, high


def score_disk(tab, stocks, profiles, mu, G):
    rational, low, high = target(tab)
    radius = np.array(tab["radius_kpc"])
    target_lo = np.nextafter(np.asarray(low, float), -np.inf)
    target_hi = np.nextafter(np.asarray(high, float), np.inf)
    ticks = math.lcm(*(q.denominator for q in rational))
    edges = np.array([int(q * ticks) for q in rational])
    dlo, dhi, unit = kernel.density_ticks(radius, stocks)
    scales = kernel.length_edges(radius[-1])
    answer = {}
    budgets = {}
    for source, (knots, values, ordinary) in profiles.items():
        slices = kernel.kernel_bins(knots, values)
        moment = slices[-1]
        lo = []
        hi = []
        free_lo = []
        free_hi = []
        for x in edges:
            lower, upper, ann_lower, ann_upper = kernel.disk_boxes(
                int(x), edges, dlo, dhi, scales, ticks, slices
            )
            if source == "positive":
                lower, upper = -upper, -lower
                ann_lower, ann_upper = -ann_upper, -ann_lower
            lo.append(lower)
            hi.append(upper)
            free_lo.append(ann_lower)
            free_hi.append(ann_upper)
        # axes: observer ring, imaged source map, physical length interval.
        lo, hi, free_lo, free_hi = map(np.asarray, (lo, hi, free_lo, free_hi))
        # All longer scales give zero, included explicitly in every hull.
        lo = np.pad(lo, ((0, 0), (0, 0), (0, 1)))
        hi = np.pad(hi, ((0, 0), (0, 0), (0, 1)))
        panels = {}
        required = {}
        for j, ref in ordinary.items():
            factor = G * mu * moment / ref
            lower = physical_bounds(lo, factor, unit)[0]
            upper = physical_bounds(hi, factor, unit)[1]
            flo = physical_bounds(free_lo, factor, unit)[0]
            fhi = physical_bounds(free_hi, factor, unit)[1]
            # Both enclosures apply to a common scale. Intersect them because
            # their geometry and mass intervals are evaluated in different orders.
            lower = np.maximum(lower, flo[:, :, None])
            upper = np.minimum(upper, fhi[:, :, None])
            allowed = (lower <= target_hi[:, None, None] + 1e-9) & (
                upper >= target_lo[:, None, None] - 1e-9
            )
            mixed = (lower.min(2) <= target_hi[:, None] + 1e-9) & (
                upper.max(2) >= target_lo[:, None] - 1e-9
            )
            independent = (flo <= target_hi[:, None] + 1e-9) & (
                fhi >= target_lo[:, None] - 1e-9
            )
            assert np.all(~allowed.any(2) | mixed) and np.all(~mixed | independent)
            panels[str(j)] = [
                dict(
                    best_one_length_upper=int(cnt.max()),
                    common_mixture_ring_upper=int(mix.sum()),
                    separate_annulus_ring_upper=int(free.sum()),
                    rings=len(radius),
                )
                for cnt, mix, free in zip(
                    allowed.sum(0), mixed.T, independent.T, strict=True
                )
            ]
            # Grant the most favorable map separately at every radius.
            envelope = fhi.max(1)
            gain = target_lo / np.maximum(envelope, 1e-16)
            w = int(np.argmax(gain))
            required[str(j)] = dict(
                minimum_peak_exposure_ratio=float(gain[w]),
                radius_kpc=float(radius[w]),
                required_inward=float(target_lo[w]),
                all_annuli_upper=float(envelope[w]),
            )
        answer[source] = panels
        budgets[source] = required
    return answer, budgets


def hierarchy(profiles, mu):
    out = {}
    for name, (knots, K, ordinary) in profiles.items():
        moments = {
            p: kernel.integral(knots, K, p, Q(knots[0]), Q(knots[-1])) for p in (1, 2)
        }
        nz = [k for k, v in zip(knots, K, strict=True) if v]
        support = [min(nz) - 1, max(nz) + 1]
        out[name] = dict(
            support=support,
            moment_u=str(moments[1]),
            moment_u2=str(moments[2]),
            log_density_enclosed_per_unit_log_mass={
                str(k): str(mu * moments[2] / q) for k, q in ordinary.items()
            },
            linear_density_enclosed_per_radius_per_unit_mass_per_length={
                str(k): str(mu * moments[1] / q) for k, q in ordinary.items()
            },
        )
    return out


def outer_hierarchy(inputs, maps, profiles, mu, G):
    """Minimum mass of the central uniform-length hierarchy spanning the
    outer half by physical radius; every radius must see the full support.
    """
    shape = hierarchy(profiles, mu)["source"]
    u_min, u_max = shape["support"]
    result = {}
    for name, tab in inputs.items():
        r, v, e, b = [
            np.array(tab[k])
            for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
        ]
        outer = r >= r[-1] / 2
        lower = float(
            max(
                0, ((np.maximum(v[outer] - 2 * e[outer], 0)) ** 2 - b[outer] ** 2).max()
            )
        )
        upper = float((((v[outer] + 2 * e[outer]) ** 2) - b[outer] ** 2).min())
        delta_L = float(r[-1] / u_min - r[outer][0] / u_max)
        supply = max((row[name][-1], label) for label, row in maps.items())
        estimates = {}
        for ref, factor in shape[
            "linear_density_enclosed_per_radius_per_unit_mass_per_length"
        ].items():
            required = lower * delta_L / float(G * Q(factor))
            estimates[ref] = dict(
                minimum_cluster_mass_msun=required,
                ratio_to_largest_map=required / supply[0],
            )
        result[name] = dict(
            r_inner_kpc=float(r[outer][0]),
            r_outer_kpc=float(r[-1]),
            outer_rings=int(outer.sum()),
            minimum_band_width_kpc=delta_L,
            required_excess_v2_band=[lower, upper],
            flat_excess_compatible=lower <= upper,
            largest_patch_supply_msun=supply[0],
            supply_map=supply[1],
            ordinary=estimates,
        )
    return result


def run():
    maps = pixel_maps()
    inputs, *_ = budget.inventory.rp.galaxy.transfer.budget.load_inputs()
    profiles, mu = source_data()
    G = Q(str(budget.inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G))
    disks = {}
    budgets = {}
    for name, tab in inputs.items():
        raw, budgets[name] = score_disk(
            tab, [x[name] for x in maps.values()], profiles, mu, G
        )
        disks[name] = {
            s: {k: dict(zip(maps, v, strict=True)) for k, v in panel.items()}
            for s, panel in raw.items()
        }
    summary = {}
    for source, (_, _, cal) in profiles.items():
        rows = {}
        for radius in cal:
            for tag in maps:
                vals = [panel[source][str(radius)][tag] for panel in disks.values()]
                rows[f"{radius}/{tag}"] = {
                    k: dict(
                        rings=sum(v[k] for v in vals),
                        whole_curve_upper=sum(v[k] == v["rings"] for v in vals),
                    )
                    for k in (
                        "best_one_length_upper",
                        "common_mixture_ring_upper",
                        "separate_annulus_ring_upper",
                    )
                }
        summary[source] = rows
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "TIDAL_SCALES.md",
                budget.LAB / "exterior_tide.json",
                HERE / "closed_inventory_maps.json",
                HERE / "source_profile.json",
                budget.LAB / "moving_records.json",
            )
        },
        hierarchy=hierarchy(profiles, mu),
        compact_outer_half=outer_hierarchy(inputs, maps, profiles, mu, G),
        summary=summary,
        disks=disks,
        exposure=budgets,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    raw = json.dumps(result, indent=2, allow_nan=False) + "\n"
    path = HERE / "tidal_scales.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    for s, rows in result["summary"].items():
        print(
            s,
            {
                k: max(v[k]["rings"] for v in rows.values())
                for k in next(iter(rows.values()))
            },
        )
