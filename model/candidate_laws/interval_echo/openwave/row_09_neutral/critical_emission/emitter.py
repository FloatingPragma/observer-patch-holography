"""Isolated material emitter, joint local depletion and distal absorption."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import betaincinv

from ...row_20_electric import experiment as prior
from ...row_20_electric.schedule import CYCLES, DELTA
from .apparatus import Guide, aperture, end_words, leakage
from .measure import draw

SHOTS = 2**26


def reflected(c, v):
    # Conjugation by the linear Majorana with the real unit coefficient v.
    d = np.eye(len(c)) - 2*np.outer(v, v)
    return d @ c @ d.T


def launch(cold, kind):
    if kind == 'sham':
        return cold.c.copy()
    if kind == 'parent':
        return reflected(cold.c, np.eye(2*cold.n)[0])
    if kind == 'packet':
        # Z on bath dimer 2 includes the first dimer's X string in global JW coordinates.
        # That first dimer is in the definite X=+ state of the decoupled medium.
        return reflected(cold.c, np.eye(2*cold.n)[2])
    raise ValueError(kind)


def selected_parent(c, sign):
    norm = 1 + sign*c[0, 1]
    if norm < 1e-12:
        return 0., None
    ci, cj = c[:, 0], c[:, 1]
    result = c + sign*(np.outer(cj, ci) - np.outer(ci, cj))/norm
    result[:2, :] = 0
    result[:, :2] = 0
    result[0, 1], result[1, 0] = sign, -sign
    np.testing.assert_allclose(result @ result, -np.eye(len(c)), atol=3e-10, rtol=0)
    return norm/2, result


def pure_overlap(ca, cb):
    product = (np.eye(len(ca)) - ca @ cb)/2
    _, logdet = np.linalg.slogdet(product)
    return float(np.exp(logdet/2))


def absorber(g, target):
    angle = target*g.period
    ann = g.receiver*aperture(angle-g.phi, 2048)
    create = g.receiver.conj()*aperture(angle+g.phi, 2048)
    vector = (g.modes.conj() @ ann + g.modes @ create)/2
    weight = float(np.vdot(ann, ann).real)
    if weight < 1e-28:
        return vector, 0., g.c
    v = (g.modes @ (ann.conj()/np.sqrt(weight))).real
    np.testing.assert_allclose(np.vdot(v, v), 1, atol=2e-11)
    return vector, weight, reflected(g.c, v)


def three_exits(c, operator, even):
    vector, weight, post = operator
    all_pass = float(np.vdot(vector, (np.eye(len(c))-1j*c) @ vector).real)
    # The product of end parity and the terminal local phase changes only an overall sign on each fixed sector.
    vac = 0. if even else weight*pure_overlap(c, post)
    p = np.array([vac, all_pass-vac, 1-all_pass])
    if p.min() < -2e-11:
        raise ValueError('assembled absorption and vacuum probabilities disagree')
    return np.maximum(0, p)


def acquisition():
    rng = np.random.default_rng(90128)
    clock = prior.phase_reader(DELTA, rng, CYCLES)
    panels = []
    for n in (23, 45):
        cold = Guide(n, 'native', 0)
        target = round(256*np.pi/((2*n+1)*(1+2/65)))/256
        for factor in (0., .25, .5, 1.):
            g = Guide(n, 'native', factor)
            operator = absorber(g, target)
            hist = []
            for kind in ('parent', 'sham', 'packet'):
                start = launch(cold, kind)
                local = []
                for tick in (0, 4, 8, 16, 32, 64, 4*n, 9*n):
                    r = np.linalg.matrix_power(g.r, tick)
                    c = r @ start @ r.T
                    local.append(dict(time=tick, joint_end_counts=draw(end_words(c), SHOTS, rng)))
                r = np.linalg.matrix_power(g.r, 4*n)
                c = r @ start @ r.T
                joint, conditional = [], []
                for s in (1, -1):
                    p, collapsed = selected_parent(c, s)
                    exits = three_exits(collapsed, operator, kind == 'sham') if p > 0 else np.array([0., 0., 1.])
                    joint.extend(p*exits)
                    conditional.append(exits.tolist())
                direct = three_exits(c, operator, kind == 'sham')
                counts = draw(joint, SHOTS, rng)
                k = counts[0]
                lower = float(betaincinv(k, SHOTS-k+1, .01/24)) if k else 0.
                hist.append(dict(launch=kind, local_records=local, joint_counts=counts,
                                 joint_lower_99=lower, uninspected_counts=draw(direct, SHOTS, rng),
                                 joint_probabilities=joint, conditional_exit_probabilities=conditional))
            gap = clock['energy']/g.period
            bound = leakage((.1*(clock['energy']-clock['simultaneous_error'])/g.period-target)*g.period, 2048)
            panels.append(dict(dimers=n, source_link_fraction=factor, target=target, parent_isolated_gap=gap,
                               lowest_coupled_energy=float(g.phi[0]/g.period),
                               bright_only_tag_upper=bound, sources=hist))
    return dict(seed=90128, shots=SHOTS, clock=clock, joint_order='parent X=+, then X=-; in each: temporal and vacuum, temporal only, temporal refused',
                source_launches=['isolated X=- emitter','same isolated X=+ emitter','isolated X=+ emitter and bath-site Z packet'],
                panels=panels)


if __name__ == '__main__':
    args = argparse.ArgumentParser()
    args.add_argument('--output', type=Path, required=True)
    dest = args.parse_args().output
    dest.write_text(json.dumps(acquisition(), indent=1, allow_nan=False)+'\n')
