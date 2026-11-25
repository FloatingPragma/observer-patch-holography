"""Compare stellar spectra, photometry and an Einstein ring on each lens."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, simpson
from scipy.special import chndtr, gammainc, gammaincinv, roots_legendre

from .. import stellar_lensing as slacs

HERE = Path(__file__).parent
G = 4.300917270e-6  # kpc (km/s)^2 / M_sun
KPC_KM = 3.08567758149e16
MASS = 1e11  # M_sun; the positive root for stellar mass uses this unit


def extract(catalogue):
    text = catalogue.read_text()
    spectra = {r["Name"]: r for r in slacs.catalogue_table(text, 4)}
    rows = []
    for ring in slacs.catalogue_table(text, 5):
        r = spectra[ring["Name"]]
        if (r["Mph"], r["Mul"], r["Lens"], ring["Good?"]) != ("E", "S", "A", "Yes"):
            continue
        if not all(r[key] for key in ("Re", "sigma", "e_sigma")):
            continue
        rows.append(
            {
                "name": r["Name"],
                **{
                    key: float(r[key])
                    for key in ("zFG", "zBG", "Re", "sigma", "e_sigma")
                },
                **{key: float(ring[key]) for key in ("bSIE", "bLTM")},
            }
        )
    archived = json.loads((HERE.parent / "receipts/stellar_lensing.json").read_text())
    if [r["name"] for r in rows] != [v[0] for v in archived["objects"]]:
        raise ValueError("The two approved catalogue selections differ")
    return rows


class Tracer:
    """Radial moments in units of a projected effective radius."""

    def __init__(self, size=901, order=96):
        self.radius = r = np.geomspace(1e-7, 1e4, size)
        self.logr = np.log(r)
        x, w = roots_legendre(order)
        self.angle = (x + 1) * np.pi / 4
        self.weights = w * np.pi / 4
        self.cos = np.cos(self.angle)
        self.sin = np.sin(self.angle)
        self.b = gammaincinv(8, 0.5)
        s = r[:, None] / self.cos
        minus_dI = self.b * s ** (-0.75) * np.exp(-self.b * s**0.25) / 4
        self.nu = (minus_dI / self.cos) @ self.weights / np.pi
        mass = cumulative_simpson(self.nu * r**3, x=self.logr, initial=0.0)
        self.fraction = np.minimum(1.0, mass / mass[-1])

    def aperture(self, re):
        r = self.radius[:, None] * self.sin * re
        response = np.zeros_like(r)
        sel = r < slacs.FIBRE_RADIUS + 12 * slacs.PSF_SIGMA
        response[sel] = chndtr(
            (slacs.FIBRE_RADIUS / slacs.PSF_SIGMA) ** 2,
            2,
            (r[sel] / slacs.PSF_SIGMA) ** 2,
        )
        radial = response @ (self.weights * self.sin * self.cos**2)
        tangential = response @ (self.weights * self.sin**3)
        return radial, tangential

    def moment(self, aperture, vc2):
        radial, tangential = aperture
        density, r, log = self.nu, self.radius, self.logr
        # The upper Jeans radial pressure has r² p_r = integral_r^infty r nu vc² dr.
        integrand = density * vc2
        p_max = (
            cumulative_simpson((r**2 * integrand)[::-1], x=-log[::-1], initial=0)[::-1]
            / r**2
        )
        light = simpson(r**3 * density * (radial + tangential), x=log)
        lower = simpson(r**3 * integrand * tangential / 2, x=log) / light
        upper = simpson(r**3 * p_max * radial, x=log) / light
        return np.array([lower, upper])

    def cylinder(self, re_kpc, ring_kpc, a0, edge_re):
        ratio = ring_kpc / re_kpc
        f2 = gammainc(8, self.b * ratio**0.25)
        rays = ratio / self.cos
        admitted = np.minimum(rays, edge_re)
        f3 = np.interp(np.log(admitted), self.logr, self.fraction)
        dark = (
            ring_kpc
            / G
            * np.sqrt(G * MASS * a0)
            * ((np.sqrt(f3) * admitted / rays) @ self.weights)
        )
        return MASS * f2, dark


def positive_mass(y, a, b):
    """Positive m with a*m + b*sqrt(m) = y; rationalized root."""
    return (2 * y / (np.sqrt(b * b + 4 * a * y) + b)) ** 2


def calc(
    rows,
    size=901,
    order=96,
    *,
    a0_endpoint=None,
    edge_re=np.inf,
    re_scale=1.0,
    ring_scale=1.0,
):
    tracer = Tracer(size, order)
    source = np.loadtxt(
        HERE.parent / "data/slacs_weak_lenses.dat", dtype=str, skiprows=2
    )
    population = {r[0]: r[5:].astype(float) for r in source}
    source = json.loads(
        (HERE.parent / "receipts/sparc_deep_regime_diagnostic.json").read_text()
    )
    a0_si = source["deep_radial_acceleration"]["a0_fixed_exponent_m_s2"]
    if a0_endpoint is not None:
        a0_si = source["deep_radial_acceleration"][
            "a0_galaxy_cluster_bootstrap_95pct_m_s2"
        ][a0_endpoint]
    a0 = a0_si * KPC_KM / 1000
    fields = []
    components = []
    for row in rows:
        lens = slacs.comoving(row["zFG"])
        efficiency = 1 - lens / slacs.comoving(row["zBG"])
        geometry = 2 * np.pi * slacs.ARCSEC_PER_RADIAN * efficiency
        angular_re = re_scale * row["Re"]
        aperture = tracer.aperture(angular_re)
        k = tracer.moment(aperture, 1.0)
        shear = ring_scale * row["bSIE"] * slacs.C_KM_S**2 / geometry
        sigma, err = row["sigma"], row["e_sigma"]
        fields.append(
            dict(
                name=row["name"],
                k=k.tolist(),
                total_weight=(k * shear / sigma**2).tolist(),
                speed_and_angle_box=(
                    k
                    * shear
                    * np.array(
                        [0.98 / (sigma + 2 * err) ** 2, 1.02 / (sigma - 2 * err) ** 2]
                    )
                ).tolist(),
            )
        )
        scale = (
            slacs.C_KM_S / 70 * 1000 * lens / (1 + row["zFG"]) / slacs.ARCSEC_PER_RADIAN
        )
        re = angular_re * scale
        vb = G * MASS / re * tracer.fraction / tracer.radius
        admitted = np.minimum(tracer.radius, edge_re)
        f3 = np.interp(np.log(admitted), tracer.logr, tracer.fraction)
        va = np.sqrt(G * MASS * a0 * f3) * admitted / tracer.radius
        kb = tracer.moment(aperture, vb)
        ka = tracer.moment(aperture, va)
        stellar_box = None
        if row["name"] in population:
            chab, dchab, salp, dsalp = population[row["name"]]
            stellar_box = (10 ** (np.array([chab - dchab, salp + dsalp]) - 11)).tolist()
        kin_box = [
            positive_mass((sigma - 2 * err) ** 2, kb[1], ka[1]),
            positive_mass((sigma + 2 * err) ** 2, kb[0], ka[0]),
        ]
        lo, hi = kin_box
        if stellar_box is not None:
            lo, hi = max(lo, stellar_box[0]), min(hi, stellar_box[1])
        rings = {}
        for label in ("bSIE", "bLTM"):
            theta = ring_scale * row[label]
            radius = theta * scale
            mb, ma = tracer.cylinder(re, radius, a0, edge_re)
            measured = (
                slacs.C_KM_S**2
                / (4 * G)
                * radius
                * theta
                / slacs.ARCSEC_PER_RADIAN
                / efficiency
            )

            def q(m):
                return (measured - mb * m) / (ma * np.sqrt(m))

            rings[label] = dict(
                radius_kpc=radius,
                ring_mass_msun=measured,
                baryon_cylinder_msun=mb,
                anomaly_cylinder_msun=ma,
                q_interval=None if lo > hi else [q(hi), q(lo)],
                pressure_column_over_c2_msun_pc2=[
                    (mb * m + ma * np.sqrt(m) - measured)
                    / (np.pi * (1000 * radius) ** 2)
                    for m in (lo, hi)
                ]
                if lo <= hi
                else None,
                speed_only_q=[q(kin_box[1]), q(kin_box[0])],
                zero_slip_stellar_mass_1e11=positive_mass(measured, mb, ma),
            )
        components.append(
            dict(
                name=row["name"],
                stellar_box_1e11=stellar_box,
                speed_box_1e11=kin_box,
                mass_intersection_1e11=[lo, hi],
                rings=rings,
                baryon_speed2=kb.tolist(),
                anomaly_speed2=ka.tolist(),
            )
        )

    def intersection(boxes):
        return [max(b[0] for b in boxes), min(b[1] for b in boxes)]

    totals = [f["total_weight"] for f in fields]
    summaries = {}
    for key in ("bSIE", "bLTM"):
        summaries[key] = {}
        for group, readout in (
            ("53_spectra", "speed_only_q"),
            ("18_spectra", "speed_only_q"),
            ("18_spectra_and_populations", "q_interval"),
        ):
            members = [
                r
                for r in components
                if group == "53_spectra" or r["stellar_box_1e11"] is not None
            ]
            boxes = [(r["name"], r["rings"][key][readout]) for r in members]
            valid = [(n, b) for n, b in boxes if b is not None]
            lower = max(valid, key=lambda v: v[1][0])
            upper = min(valid, key=lambda v: v[1][1])
            summaries[key][group] = dict(
                empty_mass_boxes=[n for n, b in boxes if b is None],
                zero_slip_outside=[n for n, b in valid if not b[0] <= 1 <= b[1]],
                maximum_lower=[lower[0], lower[1][0]],
                minimum_upper=[upper[0], upper[1][1]],
            )
    samples = [r for r in components if r["stellar_box_1e11"] is not None]
    return dict(
        a0_m_s2=a0_si,
        n_spectra=len(fields),
        n_matching_stellar_masses=len(samples),
        n_population_boxes_narrowed_by_spectra=sum(
            r["mass_intersection_1e11"] != r["stellar_box_1e11"] for r in samples
        ),
        total_flat_field=dict(
            common_nominal_interval=intersection(totals),
            includes_unit=sum(b[0] <= 1 <= b[1] for b in totals),
            two_speed_errors_and_two_percent_ring=intersection(
                [f["speed_and_angle_box"] for f in fields]
            ),
        ),
        extended_square_root=summaries,
        total_by_lens=fields,
        source_by_lens=components,
    )


def sensitivity(rows):
    scenarios = {
        "sparc_lower_scale": {"a0_endpoint": 0},
        "sparc_upper_scale": {"a0_endpoint": 1},
        "anomaly_ends_at_ten_Re": {"edge_re": 10},
        "effective_radii_minus_3.5_percent": {"re_scale": 0.965},
        "effective_radii_plus_3.5_percent": {"re_scale": 1.035},
        "ring_radii_minus_2_percent": {"ring_scale": 0.98},
        "ring_radii_plus_2_percent": {"ring_scale": 1.02},
    }
    summary = {}
    for name, options in scenarios.items():
        result = calc(rows, **options)
        groups = result["extended_square_root"]["bSIE"]
        summary[name] = {
            group: [entry["maximum_lower"][1], entry["minimum_upper"][1]]
            for group, entry in groups.items()
        }
        if any(entry["empty_mass_boxes"] for entry in groups.values()):
            raise ValueError("Separate empty mass boxes require reporting")
    return summary


def dump(obj):
    def short(value):
        if isinstance(value, (float, np.floating)):
            return float(f"{value:.9g}")
        if isinstance(value, dict):
            return {k: short(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [short(v) for v in value]
        return value

    return json.dumps(short(obj), indent=2, allow_nan=False) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = HERE / "sample.json"
    rows = extract(args.catalogue) if args.catalogue else json.loads(data.read_text())
    if args.catalogue and args.write:
        data.write_text(dump(rows))
    if rows != json.loads(data.read_text()):
        raise ValueError("archived photometry, spectra or imaging differ")
    result = calc(rows)
    result["one_at_a_time_controls"] = sensitivity(rows)
    results = dump(result)
    if args.write:
        (HERE / "results.json").write_text(results)
    elif json.loads(results) != json.loads((HERE / "results.json").read_text()):
        raise ValueError("archived virial comparison differs")
    result = json.loads(results)
    print(
        result["n_spectra"], result["total_flat_field"], result["extended_square_root"]
    )
