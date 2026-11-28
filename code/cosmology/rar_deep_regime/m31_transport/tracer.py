"""Independent M31 halo-star dispersion from the frozen baryonic and record fields."""

import argparse
from functools import lru_cache
from hashlib import sha256
from io import BytesIO
import json
import re
import tarfile
from urllib.request import urlopen

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import roots_legendre

from cosmology.rar_deep_regime.m31_transport import measure as m

# SPLASH III Table 3, Figure 7 and SPLASH I surface brightness fit for the
# dynamically hot sample. The latter's stated uncertainty is 0.3, not the
# 0.2 uncertainty for the sample retaining substructure.
TABLE = "primary_params_results.tex"
FIG = "fitparams_vs_radius_M31meanprior_normal_substprior_5sigma_normal_M31_halo_sigma.pdf"
HASH = {
    TABLE: "c19386e727014874345675a74bafd6d081a66601398145253e8b1820c5f9d3f3",
    FIG: "ffe5718bdfe78f456e6eddcaff78b1d6ce89bdbf99925f759b0b44b6edff5aaa",
    "ms.tex": "78bd2e66fa0f58d2e4895ef1bc7c43ad634228b146b92d7dcc86021168460a7f",
}
R_SPLASH = np.array([27.0590749643, 33.2049070426, 38.1767001523])
R_M31 = R_SPLASH * 785 / 783
DISPERSION = [98.0, 6.6, 7.2]


def verify_inputs():
    from pypdf import PdfReader

    raw = {}
    for ref, names in (("1711.02700v1", [TABLE, FIG]), ("1210.3362v1", ["ms.tex"])):
        with urlopen("https://arxiv.org/src/" + ref, timeout=30) as download:
            with tarfile.open(fileobj=BytesIO(download.read())) as archive:
                for name in names:
                    raw[name] = archive.extractfile(name).read()
                    assert sha256(raw[name]).hexdigest() == HASH[name]
    tex = raw[TABLE].decode("latin1")
    rows = [line.split("&") for line in tex.splitlines() if re.match(r"^\d", line)]
    overlap = [row for row in rows if float(row[0]) < 38.09 and float(row[1]) > 28]
    assert len(overlap) == 1
    cells = overlap[0]
    assert [float(s) for s in cells[:4]] == [24, 40, 382, 240]
    value, up, down = map(float, re.findall(r"\d+\.\d+", cells[7]))
    np.testing.assert_array_equal([value, down, up], DISPERSION)
    assert r"$-2.2\pm0.3$" in raw["ms.tex"].decode("latin1")
    # The first block holds seven horizontal error bars, the second the seven
    # vertical ones on the nominal-bin plot. Axis ticks fix plot coordinates.
    stream = (
        PdfReader(BytesIO(raw[FIG])).pages[0].get_contents().get_data().decode("latin1")
    )
    blocks = stream.split(" re W n /A2 gs 1 j 1.15 w 0 G")[1:3]
    pattern = r"([\d.]+) ([\d.]+) m\n([\d.]+) ([\d.]+) l\n\nS"
    horiz, vert = [np.array(re.findall(pattern, b), float) for b in blocks]
    assert len(vert) == len(horiz) == 7
    xticks = re.findall(r"([\d.]+) 52.56 m\n[\d.]+ 56.56 l", stream)
    assert "(0) Tj ET" in stream and "(200) Tj ET" in stream
    tick = np.array(xticks, float)
    np.testing.assert_allclose(tick, np.linspace(tick[0], tick[4], 5))
    xs = np.array([horiz[2, 0], vert[2, 0], horiz[2, 2]])
    np.testing.assert_allclose(
        (xs - tick[0]) * 200 / (tick[4] - tick[0]), R_SPLASH, atol=1e-8
    )


@lru_cache
def unit_nodes(n):
    x, wt = roots_legendre(n)
    return (x + 1) / 2, wt / 2


def fields():
    source = m.source(
        json.loads(m.INPUT.read_text())["HI_radius_kpc_and_surface_msun_pc2"], 0.05
    )
    setting = json.loads(m.RESULT.read_text())
    mass = 7.1e10 * source.disk + 6.4e10 * source.bulge + source.gas
    cumulative = np.r_[0, np.cumsum(mass)]

    def baryon(r):
        return m.G * (1.4e8 + np.interp(r, source.edge, cumulative)) / r

    result = {"baryons": baryon}
    for name, mu in (("finite", setting["mu_per_kpc"]), ("lossless", 0)):
        w = m.ru(source, mass / mass.sum(), mu)
        inside = cumulative_trapezoid(w**2, source.edge, initial=0)

        def v2(r, mu=mu, w=w, inside=inside):
            # Exact solution outside the last registered injection shell.
            d = np.maximum(0, r - source.edge[-1])
            tail = d if mu == 0 else -np.expm1(-2 * mu * d) / (2 * mu)
            length = np.interp(r, source.edge, inside) + w[-1] ** 2 * tail
            return (
                baryon(r)
                + m.G * setting["lossless_outer_mass_slope_msun_kpc"] * length / r
            )

        result[name] = v2
    return result


def line_of_sight(radial, R, gamma, beta, core=0, breaks=(), n=72):
    # r=R/secant denominator eliminates the singularity at the tangent point.
    u, wt = unit_nodes(n)
    bounds = np.r_[
        0, np.arccos(R / np.asarray([t for t in breaks if t > R])), np.pi / 2
    ]
    bounds = np.unique(np.sort(bounds))
    theta = (bounds[:-1, None] + np.diff(bounds)[:, None] * u).ravel()
    coeff = (np.diff(bounds)[:, None] * wt).ravel()
    c = np.cos(theta)
    r = R / c
    weight = coeff * c ** (gamma - 2) * (1 + (core / r) ** 2) ** (-gamma / 2)
    return weight @ ((1 - beta * c * c) * radial(r)) / weight.sum()


def jeans(v2, R, gamma=3.2, beta=0, core=0, force_beyond=0, n=72):
    # Substitute t=r/s in the right-hand Jeans integral. For the optional
    # support diagnostic integrate only force samples at s >= force_beyond.
    u, wt = unit_nodes(n)

    def radial(r):
        top = np.minimum(1, r / force_beyond) if force_beyond else np.ones_like(r)
        t = top[:, None] * u
        core_weight = (
            (1 + (core / r[:, None]) ** 2) / (1 + (core * t / r[:, None]) ** 2)
        ) ** (gamma / 2)
        return (
            (v2(r[:, None] / t) * t ** (gamma - 2 * beta - 1) * core_weight) @ wt * top
        )

    return line_of_sight(radial, R, gamma, beta, core=core, breaks=[force_beyond], n=n)


def mass_floor(rows, shift):
    r, v, e, vb = (
        np.array([row[k] for row in rows])
        for k in ("r_kpc", "observed_kms", "quoted_error_kms", "ordinary_speed")
    )
    gm = np.maximum.accumulate(np.maximum(0, r * ((v - shift * e) ** 2 - vb**2)))
    return r, np.diff(gm, prepend=0)


def floor_extra(rings, increments, R, gamma, beta):
    # Integrating a nondecreasing spherical enclosed mass gives an explicit
    # sum of exterior Kepler responses. Every increment has positive weight.
    h = gamma + 1 - 2 * beta

    def radial(r):
        return (np.minimum(1, r[:, None] / rings) ** h) @ increments / (h * r)

    return line_of_sight(radial, R, gamma, beta, breaks=rings)


def run():
    compiled = fields()
    gas = json.loads(m.RESULT.read_text())["surveys"]
    cases = {}
    edges = np.unique(np.r_[np.linspace(R_M31[0], R_M31[-1], 25), R_M31])
    at = np.searchsorted(edges, R_M31)
    for gamma in (2.9, 3.2, 3.5):
        for beta in (-0.5, 0, 0.5):
            case = {}

            # Positive spherical enclosed mass makes R*sigma_los^2(R)
            # nondecreasing for a power-law tracer with constant anisotropy.
            # The two adjacent endpoint rectangles bound every radius.
            def support_bounds(variance):
                lo = np.min(edges[:-1] * variance[:-1] / edges[1:])
                hi = np.max(edges[1:] * variance[1:] / edges[:-1])
                return [np.sqrt(lo), np.sqrt(hi)]

            for name, v2 in compiled.items():
                variance = np.array([jeans(v2, R, gamma, beta) for R in edges])
                case[name] = np.sqrt(variance[at]).tolist()
                case[name + "_all_sampled_radii_bounds"] = support_bounds(variance)
                if name == "baryons":
                    ordinary = variance
            for survey, profile in gas.items():
                # Conditional continuation of the outermost measured speed
                # throughout the tracer's force kernel, not a measured halo.
                flat_speed = profile["rows"][-1]["observed_kms"]
                case[survey + "_constant_endpoint_speed_kms"] = flat_speed * np.sqrt(
                    (1 - beta * (gamma - 1) / gamma) / (gamma - 2 * beta)
                )
                cases_floor = []
                lower_bounds = []
                for shift in (0, 1, 2):
                    rings, inc = mass_floor(profile["rows"], shift)
                    variance = ordinary + np.array(
                        [floor_extra(rings, inc, R, gamma, beta) for R in edges]
                    )
                    cases_floor.append(np.sqrt(variance[at]).tolist())
                    lower_bounds.append(support_bounds(variance)[0])
                case[survey + "_minimum_by_zero_one_two_error"] = cases_floor
                case[survey + "_all_sampled_radii_minimum_by_zero_one_two_error"] = (
                    lower_bounds
                )
            cases[f"gamma{gamma}_beta{beta:g}"] = case
    R = R_M31[1]
    support = {
        str(limit): jeans(compiled["finite"], R, force_beyond=limit)
        / jeans(compiled["finite"], R)
        for limit in (38.09, 90, 175)
    }
    cored = {name: np.sqrt(jeans(v2, R, core=5.2)) for name, v2 in compiled.items()}
    return {
        "source_hashes": HASH,
        "protocol_sha256": sha256(
            (m.HERE / "TRACER_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "observed_sigma_minus_plus_kms": DISPERSION,
        "catalog_bin_kpc": [24, 40],
        "counts_all_likely_M31": [382, 240],
        "plotted_support_median_kpc": R_SPLASH.tolist(),
        "common_distance_support_median_kpc": R_M31.tolist(),
        "cases_in_order_min_median_max_radius": cases,
        "isotropic_finite_variance_fraction_from_force_beyond_kpc": support,
        "isotropic_cored_tracer_at_median_kms": cored,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    if args.write and not args.fetch:
        p.error("writing requires verification of the primary stellar records")
    if args.fetch:
        verify_inputs()
    out = run()
    path = m.HERE / "tracer.json"
    if args.write:
        with path.open("x") as handle:
            json.dump(out, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        m.stationary.check(json.loads(path.read_text()), out)
    print(out["cases_in_order_min_median_max_radius"]["gamma3.2_beta0"])
