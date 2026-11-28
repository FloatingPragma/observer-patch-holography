"""A pressure-history bound from conserved energy and observed radial BAO."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from statistics import NormalDist

import numpy as np

HERE = Path(__file__).resolve().parent
PREFIX = "desi_gaussian_bao_ALL_GCcomb_"
C = 299792.458


def endpoint_test(s, ratios, sigma, cap):
    """Monotone endpoint inequalities in the published distance coordinates."""
    zscore = NormalDist().inv_cdf(1 - 0.01 / 4)
    # H increases as D_H/r_d decreases. Take opposite corners and the
    # longest ruler to minimize the inferred expansion difference.
    low_d = ratios[0] - zscore * sigma[0]
    high_d = ratios[1] + zscore * sigma[1]
    shape = 1 / high_d**2 - 1 / low_d**2
    assert shape > 0
    ruler = 165
    omega_b, omega_nu, omega_r = 0.02273, 0.0007, 0.00005
    hdiff = (C / (100 * ruler)) ** 2 * shape
    d2, d3, d4 = [s[1] ** k - s[0] ** k for k in (2, 3, 4)]
    total = sum(s**3)
    known = (omega_b + omega_nu) * d3 + omega_r * d4
    extra = hdiff - known
    return dict(
        distance_intervals_99=np.column_stack(
            (ratios - zscore * sigma, ratios + zscore * sigma)
        ).tolist(),
        largest_ruler_Mpc=ruler,
        minimum_delta_H_squared_over_100_squared=hdiff,
        maximum_known_difference=known,
        minimum_constant_cold_equivalent_per_baryon=extra / (omega_b * d3),
        minimum_absolute_endpoint_density_per_baryon=extra / (omega_b * total),
        minimum_open_curvature_omega_k=(extra - omega_b * cap * total) / d2,
        minimum_flat_ruler_Mpc=C
        / 100
        * np.sqrt(shape / (known + omega_b * cap * total)),
    )


def run():
    inv = HERE / "closed_inventory.json"
    source = json.loads(inv.read_bytes())
    mean, cov = [HERE / (PREFIX + s + ".txt") for s in ("mean", "cov")]
    data = np.loadtxt(mean, dtype=str, skiprows=1)
    ids = np.flatnonzero(data[:, 2] == "DH_over_rs")
    ids = ids[np.argsort(data[ids, 0].astype(float))][[0, -1]]
    redshift, dist = data[ids, :2].astype(float).T
    sigma = np.sqrt(np.diag(np.loadtxt(cov))[ids])
    s = 1 + redshift
    target = 12000 / 2237
    caps = {
        "four_owners": source["maximum_extra_per_ordinary"],
        "all_eligible": source["all_accepted_stock_ratio_upper"],
    }
    out = {}
    for name, cap in caps.items():
        # a^(-1)e changes by the integral of A d(1/a), where e=rho_x/rho_b
        # and A=(rho_x+3p_x)/rho_b. Use independent signs at both endpoints.
        active = cap * sum(s) / (s[1] - s[0])
        out[name] = dict(
            energy_cap=cap,
            maximum_mean_active_per_baryon=active,
            maximum_scale_factor_ratio_at_target=(target + cap) / (target - cap),
            endpoints=endpoint_test(s, dist, sigma, cap),
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "PAIR_PRESSURE.md", inv, mean, cov)
        },
        radial_redshifts=redshift.tolist(),
        observed_DH_over_rd=dist.tolist(),
        quoted_sigma=sigma.tolist(),
        scale_factor_ratio=s[1] / s[0],
        planck_base_cold_ratio=target,
        populations=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "pair_pressure.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from .fluctuation_cuts import rp

        rp.galaxy.check(json.loads(dst.read_bytes()), result)
    print(json.dumps(result["populations"], indent=2))
