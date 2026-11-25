"""Exact entropy decomposition for repeated surveys, with a two-board native lottery check."""

from pathlib import Path
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import xlogy

from . import station


def head_probabilities(q, z):
    weights = np.full(4, 1/32)
    weights[2] += 7*z/8
    weights[3] += 7*(1-z)/8
    return np.r_[(1-q)*weights, q*weights]


def classical_entropy(probabilities):
    return float(-xlogy(probabilities, probabilities).sum())


def measure():
    folder = Path(__file__).parent
    single = json.loads((folder/'results.json').read_text())
    equilibrium = np.array(single['x_real']) + 1j*np.array(single['x_imag'])
    effects = station.radar_menus()[6]
    diagonal = effects.diagonal(axis1=-2, axis2=-1).real.reshape(12, 8)
    eta, target = single['ordinary_eta'], single['target_volume']

    def survey(probabilities):
        return station.geometry(diagonal @ probabilities)

    # The two controllers share one fair selector, and use independent native lotteries thereafter.
    q1, q2 = 1/64, 1/32
    z1, spread = .77, .1
    def joint_lottery(second_control):
        return sum(np.outer(head_probabilities(q1, z1+sign*spread),
                            head_probabilities(q2, second_control+sign*spread))/2 for sign in (-1,1))

    p = head_probabilities(q1, z1)
    a1, v1 = survey(p)
    z2 = brentq(lambda z: survey(head_probabilities(q2, z))[1] + v1 - 2*target, .1, .9, xtol=1e-14)
    q = head_probabilities(q2, z2)
    a2, v2 = survey(q)
    joint = joint_lottery(z2)
    equilibrium_area = station.geometry(station.pair(effects, equilibrium))[0]
    source_entropy = station.entropy(equilibrium)
    gap = 2*(source_entropy+eta*equilibrium_area) - classical_entropy(joint) - eta*(a1+a2)
    total_correlation = classical_entropy(p)+classical_entropy(q)-classical_entropy(joint)
    separate_gaps = sum(source_entropy-classical_entropy(w)+eta*(equilibrium_area-a)
                        for w, a in ((p, a1), (q, a2)))
    log_equilibrium = station.log_blocks(equilibrium).diagonal(axis1=-2,axis2=-1).real.ravel()
    relative_entropy = -classical_entropy(joint) - (p+q) @ log_equilibrium
    assert abs(v1+v2-2*target) < 2e-14 and abs(v1-target) > 5e-4
    assert max(np.max(abs(joint.sum(axis=1)-p)), np.max(abs(joint.sum(axis=0)-q))) < 1e-14
    assert total_correlation > 0 and abs(gap-total_correlation-separate_gaps) < 1e-12
    assert gap > .28*relative_entropy
    pooled_z = brentq(lambda z: survey((p+head_probabilities(q2, z))/2)[1]-target, .1, .9, xtol=1e-14)
    q_pooled = head_probabilities(q2, pooled_z)
    mean = (p+q_pooled)/2
    pooled_area, pooled_volume = survey(mean)
    pooled_joint = joint_lottery(pooled_z)
    pooled_entropy = classical_entropy(pooled_joint)
    pooled_gap = 2*(source_entropy+eta*(equilibrium_area-pooled_area)) - pooled_entropy
    pooled_relative = -pooled_entropy - 2*mean @ log_equilibrium
    pooled_correlation = 2*classical_entropy(mean)-pooled_entropy
    predicted_gap = pooled_correlation + 2*(source_entropy-classical_entropy(mean)
                                              +eta*(equilibrium_area-pooled_area))
    separate_offset = v1+survey(q_pooled)[1]-2*target
    assert abs(pooled_volume-target) < 2e-14 and abs(separate_offset) > 1e-7
    assert abs(pooled_gap-predicted_gap) < 1e-12 and pooled_gap > .28*pooled_relative > 0
    return dict(pooled_volume_check=dict(second_mean_control=pooled_z, volume_offset=pooled_volume-target,
                                         separate_total_volume_offset=separate_offset,
                                         objective_gap_nats=pooled_gap, lower_bound_nats=.28*pooled_relative),
                odd_fractions=[q1, q2], mean_controls=[z1, z2], shared_control_step=spread,
                individual_volume_offsets=[v1-target,v2-target], total_volume_offset=v1+v2-2*target,
                total_correlation_nats=total_correlation, objective_gap_nats=gap,
                relative_entropy_to_equilibrium_product_nats=float(relative_entropy),
                lower_bound_nats=float(.28*relative_entropy))


if __name__ == '__main__':
    print(json.dumps(measure(), indent=2))
