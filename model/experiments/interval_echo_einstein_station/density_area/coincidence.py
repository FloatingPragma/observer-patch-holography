"""Calibrate the curvature of boundary log multiplicity with two disjoint tagged rays."""

import argparse
from copy import deepcopy
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from .measure import HERE, class_law, conditional_mean_energy, weighted_score
from .. import assert_same_replay


TAG_SHOTS, INFO_SHOTS = 2**38, 2**34


def pair_laws(k, m, p):
    single = p*(m-k)/m
    double = p*p*(m-k)*(m-k-1)/(m*(m-1))
    return np.column_stack([1-2*single+double, single-double, single-double, double])


def exposure(v, m, p, rng, meter):
    k, law, choices = class_law(v, m, v//2)
    flags = pair_laws(k, m, p)
    width = max(1,math.isqrt(m)//2)
    tests = []
    for weights in ([.5,-1.,.5], [1/8,-4/8,6/8,-4/8,1/8]):
        half = len(weights)//2
        indices = np.searchsorted(k, m//2 + width*np.arange(-half,half+1))
        transfer = np.zeros(len(k))
        transfer[indices] = weights
        step = min(law[indices])/4
        arms = []
        for sign in (-1,1):
            mixture = law + sign*step*transfer
            flag_prob = mixture @ flags
            hits = rng.multinomial(TAG_SHOTS,flag_prob)
            interior = rng.multinomial(INFO_SHOTS,mixture)
            bs, se = weighted_score(interior,-choices)
            arms.append(dict(flag_counts=hits.tolist(), interior_counts=interior.tolist(),
                             bulk_score=bs, bulk_score_se=se))
        change = (arms[1]['flag_counts'][3]-arms[0]['flag_counts'][3])/(2*step*TAG_SHOTS)
        factor = m*(m-1)/p**2
        se = math.sqrt(sum((prob:=a['flag_counts'][3]/TAG_SHOTS)*(1-prob)
                           for a in arms)/TAG_SHOTS)/(2*step)
        energy_slope = sum(Fraction(w)*conditional_mean_energy(v,m,v//2,int(index),meter)
                           for w,index in zip(weights,k[indices],strict=True))
        tests.append(dict(exact_mean_energy_slope=str(energy_slope), classes=k[indices].tolist(), weights=weights, displacement=float(step), arms=arms,
                          flag_slopes=(transfer@flags).tolist(),
                          exact_variance_slope=float(transfer@(k-m/2)**2),
                          measured_variance_slope=change*factor, variance_se=se*factor,
                          exact_bulk_slope=float(-transfer@choices),
                          bulk_slope_from_tag_variance=2*change*factor/m,
                          bulk_prediction_se=2*se*factor/m,
                          measured_bulk_slope=(arms[1]['bulk_score']-arms[0]['bulk_score'])/(2*step),
                          measured_bulk_se=math.hypot(*(a['bulk_score_se'] for a in arms))/(2*step)))
    return tests


def run():
    original = json.loads((HERE/'results.json').read_text())
    panels = []
    for pane in original['panels']:
        parent, depth = pane['parent_level'], pane['microscopic_subdivisions']
        for index, row in enumerate(pane['geometries']):
            source = next(s for s in row['sources'] if s['filling']=='half')
            tests = exposure(source['v'],source['m'],original['isolated_tag_probability'],
                             np.random.default_rng([761937,parent,depth,index]), source['energy_meter'])
            panels.append(dict(parent=parent, depth=depth, shape=row['name'], boundary_contacts=source['m'],
                               variance_transfer=tests[0], fixed_two_ray_law=tests[1]))
    return dict(seed=761937, tag_trials_per_arm=TAG_SHOTS, interior_trials_per_arm=INFO_SHOTS,
                panels=panels)


def verify_replay(actual, stored):
    """Check both finite differences at their own recorded scores, then compare every input."""
    def primitive_readings(record):
        record = deepcopy(record)
        for panel in record["panels"]:
            for arm in ("variance_transfer", "fixed_two_ray_law"):
                test = panel[arm]
                left, right = test["arms"]
                slope = (right["bulk_score"]-left["bulk_score"])/(2*test["displacement"])
                if not math.isfinite(slope) or slope != test.pop("measured_bulk_slope"):
                    raise ValueError("Reported bulk slope differs from its two sampled scores")
        return record

    assert_same_replay(primitive_readings(actual), primitive_readings(stored))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    encoded = json.dumps(run(), separators=(',', ':'), allow_nan=False)
    with args.output.open('x') as f:
        f.write(encoded+'\n')
