"""Photometric enclosure and speed requirements for the AGC 114905 transfer."""
import argparse
import hashlib
import json
import math

from scipy.integrate import quad
from scipy.optimize import brentq

import fit

# Mancera Piña et al., A&A 689, A344, sections 3 and 5, fiducial geometry.
KPC_PER_ARCSEC = 78_700 * math.pi / (180 * 3600)
GAS_FRACTION = 0.94
HALF_STELLAR_MASS_KPC = 2.80
G_KPC_KMS2_PER_MSUN = 1.3271244e20 / (fit.base.KPC_M * 1e6)


def gas_annulus(angle):
    # Common factors in 2 pi R Sigma(R) cancel in an enclosed fraction.
    return angle * math.exp(-angle / 0.676 + 559.231 * math.log1p(angle / 375.220))


def gas_fraction_inside(radius):
    total = quad(gas_annulus, 0, math.inf)[0]
    return quad(gas_annulus, 0, radius / KPC_PER_ARCSEC)[0] / total


def exterior_orbits(a0, spacing, training_galaxies):
    archive = fit.ROOT / 'outer_orbit_transfer'
    selection = json.loads((archive / 'results.json').read_text())
    inputs = (archive / 'inputs.json').read_bytes()
    if hashlib.sha256(inputs).hexdigest() != selection['input_sha256']:
        raise ValueError("Outer-orbit input differs from the enclosure study")
    results = {}
    for row in json.loads(inputs)['galaxies']:
        name = row['name']
        region = selection['enrollment'][name]
        if not (region['radial_pass'] and region['mass80_pass'] and region['deep']):
            continue
        if name in training_galaxies:
            raise ValueError("Outer test galaxy appears in the rotation training data")
        mass = 0.8*row['stellar_unit_msun'] + row['helium_hi_msun']
        radius = row['radius_kpc']
        gm = G_KPC_KMS2_PER_MSUN * mass
        excess = math.sqrt(gm * a0 * fit.base.KPC_M / 1e6)
        baryonic = gm / radius
        counts = (math.floor(radius/spacing), math.ceil(radius/spacing))
        results[name] = dict(mass90_pass=region['mass90_pass'], radius_kpc=radius,
            enclosure_lower_bound=region['enclosed_fraction_lower'],
            observed_kms=row['speed_kms'], observed_error_kms=row['speed_error_kms'],
            smooth_point_mass_kms=math.sqrt(baryonic+excess),
            phase_conditioned_speed_range_kms=[
                math.sqrt(baryonic+excess*n*spacing/radius) for n in counts])
    return results


def run():
    bank = json.loads((fit.HERE / 'results.json').read_text())
    spacing = bank['selected_spacing_kpc']
    a0 = next(s['a0_m_s2'] for s in bank['full_sample_profile'] if s['spacing_kpc'] == spacing)
    panels = json.loads((fit.ROOT / 'udg_transfer/rings.json').read_text())['panels']
    outer = {}
    for panel, inclination in [('figure8', 31.0), ('figure7', 35.86)]:
        ring = panels[panel]['rings'][-1]
        r, v, b = (ring[k] for k in ('radius_kpc', 'circular_kms', 'baryonic_kms'))
        smooth_excess = math.sqrt(b*b*r*fit.base.KPC_M*a0/1e6)
        least_excess = math.floor(r/spacing)*spacing/r*smooth_excess
        min_speed = math.sqrt(b*b + least_excess)
        needed_excess = v*v - b*b
        outer[panel] = dict(radius_kpc=r, baryonic_kms=b, circular_kms=v,
            min_phase_model_kms=min_speed, central_excess_fraction_of_min=needed_excess/least_excess,
            implied_to_fitted_smooth_a0=(needed_excess/smooth_excess)**2,
            projected_inclination_deg=math.degrees(math.asin(v*math.sin(math.radians(inclination))/min_speed)),
            spherical_equivalent_required_excess_msun=r*needed_excess/G_KPC_KMS2_PER_MSUN,
            spherical_equivalent_min_excess_msun=r*least_excess/G_KPC_KMS2_PER_MSUN)
    outer_radius = outer['figure8']['radius_kpc']
    enclosed = gas_fraction_inside(outer_radius)
    a0_kpc = a0 * fit.base.KPC_M / 1e6
    return dict(a0_m_s2=a0, spacing_kpc=spacing, outer=outer,
        exterior_orbits=exterior_orbits(a0, spacing, bank['held_log_density_by_galaxy']),
        outer_radius_stellar_half_mass_units=outer_radius/HALF_STELLAR_MASS_KPC,
        gas_mass_fraction_in_outer_cylinder=enclosed,
        baryonic_mass_fraction_upper_in_outer_cylinder=1-GAS_FRACTION*(1-enclosed),
        gas_90pct_radius_kpc=brentq(lambda radius: gas_fraction_inside(radius)-0.90, 1, 50),
        compact_asymptotic_anomalous_speed_kms=(G_KPC_KMS2_PER_MSUN*1.47e9*a0_kpc)**0.25)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    results = run()
    path = fit.HERE / 'disk_domain.json'
    if args.write:
        path.write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    else:
        from lensing_profile import check
        check(json.loads(path.read_text()), results)
    print(json.dumps(results, indent=2))
