"""Energy and photon-number budgets of a specified pre-recombination work reservoir."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import numpy as np
from scipy import constants as si
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock import (
    trail,
)

HERE = Path(__file__).parent
# Standard post-BBN benchmarks: Fixsen (2009), DESI DR2 II eq. 14,
# Planck VI base parameter summary and N_eff extension, Froustey et al.
T0, DT = 2.72548, 0.00057
BC, DBC = 0.0224, 0.0001
BB, DBB = 0.02218, 0.00055
NSTAR = 3.044
NEFF, DNEFF = 2.99, 0.17


def photon_density(temperature):
    h100 = 100000 / (1e6 * si.parsec)
    rho100 = 3 * h100**2 / (8 * np.pi * si.G)
    return 4 * si.Stefan_Boltzmann * temperature**4 / (si.c**3 * rho100)


def history(q, baryon, redshifts, shares, photons):
    a = 1 / (1 + np.asarray(redshifts, dtype=float))
    f = np.array(shares, dtype=float)
    if min(f) < 0 or not np.isclose(f.sum(), 1):
        raise ValueError("Unnormalized stock history")
    lost = q * baryon * float(a @ f) / photons
    return dict(
        lost_over_final_comoving_photon_energy=lost,
        original_fraction_removed=lost / (1 + lost),
        baryon_to_photon_jump=(1 + lost) ** 0.75,
        neff_after=NSTAR * (1 + lost),
    )


def run():
    photons = photon_density(T0)
    # Parameter-sensitivity boxes, not a joint confidence limit for a new cosmology.
    photons_hi = photon_density(T0 + 2 * DT)
    bb_lo = BB - 2 * DBB
    bc_hi = BC + 2 * DBC
    max_by_eta = (bc_hi / bb_lo) ** (4 / 3) - 1
    max_by_neff = (NEFF + 2 * DNEFF) / NSTAR - 1
    out = []
    for b in json.loads((trail.HERE / "trail.json").read_text())["boards"]:
        q = float(Fraction(b["net_work"]) / Fraction(b["initial_energy"]))
        # For the baryon jump the smallest z comes at high final baryon density:
        # B / ((B/B0)^(4/3)-1) is decreasing in B for B>B0.
        bounds = {}
        for k, maximum, after in (
            ("baryon_ratio", max_by_eta, bc_hi),
            ("neutrino_ratio", max_by_neff, BC - 2 * DBC),
        ):
            upper = maximum * photons_hi / (q * after)
            bounds[k] = dict(
                allowed_energy_over_final=maximum,
                max_work_weighted_scale_factor=upper,
                minimum_single_burst_redshift=1 / upper - 1,
                max_fraction_since_z10000=min(1.0, 10001 * upper),
            )
        out.append(
            dict(
                heads=b["live_seams"],
                ratio=q,
                calibration_bounds=bounds,
                snapshots={
                    str(z): history(q, BC, [z], [1], photons)
                    for z in (1100, 1e4, 5e4, 1e5, 1e6)
                },
                mixed_epoch_control=history(q, BC, [2000, 1e5], [0.1, 0.9], photons),
            )
        )
    return dict(
        omega_gamma=photons,
        parameters=dict(
            temperature_K=T0,
            bbn_baryons=BB,
            cmb_baryons=BC,
            neff_standard=NSTAR,
            planck_bao_effective=NEFF,
        ),
        boards=out,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    r = run()
    p = HERE / "results.json"
    if args.write:
        p.write_text(json.dumps(r, indent=2) + "\n")
    else:
        from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture.measure import (
            check,
        )

        check(json.loads(p.read_text()), r)
