"""Equal-proper-volume rim measurements at the common initial scale factor."""
import argparse
import json
from math import pi

import numpy as np

from . import measure
from ..quantum_scalar_radar import balls


def metric_at_preparation(energy):
    # The exact expanding dust solution is independent of the radar's ODE integration.
    kappa, vacuum = measure.radar.COUPLING, measure.radar.VACUUM
    omega = np.sqrt(3*vacuum)/2
    initial = np.sqrt(measure.V0)
    growth = np.sqrt(measure.V0+kappa*energy/vacuum)

    def source(times):
        z = omega*(np.asarray(times)-balls.APPOINTMENT)
        y = initial*np.cosh(z) + growth*np.sinh(z)
        if np.any(y <= 0):
            raise ValueError('Rim survey crossed the singularity')
        volume = y*y
        q = np.log(volume/measure.V0)/3
        expansion = 2*omega*(initial*np.sinh(z)+growth*np.cosh(z))/(3*y)
        momentum = -2*volume*expansion/kappa
        return np.stack([q, q, q, momentum, momentum, momentum])
    return source


def run(primary):
    arm = primary['arms']['interacting']
    kappa = arm['comparisons'][-1]['ordinary_coupling']
    high, low = [balls.survey_source(metric_at_preparation(e)) for e in arm['source_energies']]
    control = balls.survey_source(metric_at_preparation(primary['selective_gravity_control']['gravitating_energy']))
    rows = []
    for a, b, sham in zip(high, low, control, strict=True):
        radius = a['reference']
        assert radius == b['reference'] == sham['reference']
        factor = (8*pi*pi/15)*radius**4/measure.V0
        area = a['area']-b['area']
        population_charge = factor*arm['exact_overlap_energy']
        counted_charge = factor*arm['overlap_meter']
        rows.append(dict(radius=radius, area_difference=area,
                         measured_curvature=a['estimated_Guu']-b['estimated_Guu'],
                         population_weighted_charge=population_charge,
                         counted_weighted_charge=counted_charge,
                         forecast_area_from_counts=-kappa*counted_charge/(2*pi),
                         forecast_counting_se=kappa*factor*arm['overlap_counting_se']/(2*pi),
                         inferred_coupling=-2*pi*area/population_charge,
                         high=a, low=b, constant_charge_common_survey=sham))
    return dict(common_cell_volume=measure.V0, ordinary_radar_coupling=kappa, rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    results = run(json.loads((measure.HERE/'results.json').read_text()))
    target = measure.HERE/'rim.json'
    if parser.parse_args().write:
        target.write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    else:
        from cosmology.rar_deep_regime.lensing_profile import check
        check(json.loads(target.read_text()), results)
    for row in results['rows']:
        print(row['radius'], row['area_difference'], row['inferred_coupling'])
