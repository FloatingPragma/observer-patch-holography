"""Infer the ordinary energy-to-area response on a stationary survey and test its transfer."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from matched_einstein_limit.source_implication.full_algebra import enclosure
from matched_einstein_limit.source_implication.shared_radar import measure as radar

from . import returned_sources as scores
from . import six_station_read as six

HERE = six.HERE
iv = enclosure.iv


def run():
    iv.dps = 35
    source = HERE / "returned_head_lab.json"
    cones = HERE / "erasure_balls.json"
    certificate = Path(enclosure.__file__).with_name("certificate.json")
    observed = Path(radar.__file__).with_name("results.json")
    files = [source, cones, certificate, observed]
    raw = json.loads(source.read_bytes())

    # Parity is the energy outcome on the parked two-head instrument (16,17).
    def read(label):
        return scores.observed(
            raw["tapes"][label]["incoming"], raw["shots_per_setting"], False
        )

    X, e = read("reference")
    ordinary, f = read("ordinary")
    dO = ordinary - X
    odd = float(dO[1].trace().real)
    energy = iv.mpf([odd - e[1] - f[1], odd + e[1] + f[1]])
    entropy = iv.mpf(0)
    for b in (0, 1):
        eig, U = np.linalg.eigh(X[b])
        dx = 2 * e[b]
        K = (U * (-np.log(eig))) @ U.conj().T
        coeff = np.abs(np.einsum("ijkl,lk->ij", scores.heads.PRODUCTS, K)).sum() / 4
        ed = f[b] + e[b]
        norm = np.linalg.norm(dO[b], ord="nuc") + 4 * ed
        width = coeff * ed + norm * dx / (eig[0] - dx) + 1e-13
        value = float(np.trace(K @ dO[b]).real)
        entropy += iv.mpf([value - width, value + width])
    cert = json.loads(certificate.read_bytes())
    eta = iv.mpf(cert["gain_interval"])
    stationary = {}
    for name, bounds in (("ordinary_entropy", entropy), ("ordinary_energy", energy)):
        known = iv.mpf(cert["signed_transfers"][name].strip("[]").split(","))
        assert known.a <= bounds.b and bounds.a <= known.b

    # The interval certificate locates an exact nearby constrained
    # equilibrium. The count comparison is on its numerical representative.
    # Enclose both for an energy-to-area prediction that covers that difference.
    def hull(x, name):
        other = iv.mpf(cert["signed_transfers"][name].strip("[]").split(","))
        return iv.mpf([min(x.a, other.a), max(x.b, other.b)])

    beta = hull(entropy, "ordinary_entropy") / hull(energy, "ordinary_energy")
    prediction = -beta / eta
    records = json.loads(observed.read_bytes())
    for label, v in records["charges"].items():
        chosen = v["exposures"][str(2**50)]["matched"]["ordinary"]
        val = chosen["survey_area_slope"]
        err = chosen["simultaneous_sampling_halfwidth"][0]
        slope = iv.mpf([val - err, val + err]) / energy
        stationary[label] = six.bounds(slope)
        assert slope.b < 0 and slope.a <= prediction.b and prediction.a <= slope.b
    balls = json.loads(cones.read_bytes())["balls"]
    # At equal reference radius ell, V is fixed across source arms within
    # each comparison. A smooth, bounded stress requires the area slope
    # to scale as ell^4 as ell tends to zero.
    results = []
    for row in balls:
        area = iv.mpf(row["area_slopes"][0])
        assert not area.a <= 0 <= area.b
        rates = []
        for ell in (1 / 2, 1 / 4, 1 / 8, 1 / 16):
            small = area * (ell / 0.5) ** 2
            rates.append(
                {
                    "radius": ell,
                    "fixed_volume_area_slope": six.bounds(small),
                    "slope_per_radius_fourth": six.bounds(small / ell**4),
                }
            )
        results.append(
            {"triangle": row["triangle"], "port": row["center_port"], "scales": rates}
        )
    return {
        "sha256": {
            str(p.relative_to(HERE.parents[3])): sha256(p.read_bytes()).hexdigest()
            for p in files
        },
        "parked_energy_slope_from_counts": six.bounds(energy),
        "entropy_slope_from_counts": six.bounds(entropy),
        "entropy_per_parked_energy": six.bounds(beta),
        "ordinary_area_gain": six.bounds(eta),
        "fixed_volume_area_per_parked_energy": six.bounds(prediction),
        "independent_radar_area_per_parked_energy": stationary,
        "contracting_metric_balls": sum(r["area_slopes"][0][1] < 0 for r in balls),
        "expanding_metric_balls": sum(r["area_slopes"][0][0] > 0 for r in balls),
        "piecewise_euclidean_ball_refinement": results,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dst = HERE / "ordinary_stationarity.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dst.read_bytes()), result)
    print(
        {
            k: v
            for k, v in result.items()
            if k not in ("sha256", "piecewise_euclidean_ball_refinement")
        }
    )
