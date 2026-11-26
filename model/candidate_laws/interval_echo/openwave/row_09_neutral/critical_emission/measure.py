"""Frozen local source, end reads and energy-tagged remote detector cohorts."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import betaincinv

from model.executable_model.carrier import federation
from ...row_20_electric import experiment as prior
from ...row_20_electric.schedule import CYCLES, DELTA
from .apparatus import Guide

LOCAL_SHOTS, CLOCK_SHOTS = 2**24, 2**26


def draw(weights, n, rng):
    p = np.asarray(weights)
    if p.min() < -2e-12 or abs(sum(p) - 1) > 3e-12:
        raise ValueError(f'invalid complete detector distribution {p}')
    p = np.maximum(p, 0)
    return rng.multinomial(n, p/p.sum()).tolist()


def acquisition():
    board = federation(2)
    _, path, _ = prior.guide(board)
    j, h = prior.couplings(8)
    rng = np.random.default_rng(90127)
    cal = prior.phase_reader(2*h*DELTA, rng, CYCLES)
    source_tag = {str(q): prior.sign_read(q) for q in (0, 1, 2)}
    panels = []
    for n in (23, 45):
        ports = path[1:1+2*n]
        if len(ports) != 2*n or any(board.degree[a] != board.degree[b] for a, b in zip(ports[::2], ports[1::2])):
            raise ValueError('missing equal-reference physical dimers')
        for mode in ('native', 'half', 'cut'):
            g = Guide(n, mode)
            target = round(256*np.pi/((2*n+1)*(1+2/65))) / 256
            threshold = .1*(cal['energy']-cal['simultaneous_error'])/g.period
            local = []
            for t in (0, 1, n//2, 4*n, 7*n, 8*n, 9*n, 10*n, 12*n):
                rec = [draw(g.end_read(t, sham), LOCAL_SHOTS, rng) for sham in (False, True)]
                local.append(dict(periods=t, counts=rec,
                                  parent_survival=float(abs(abs(g.source)**2 @ np.exp(-1j*t*g.phi))**2)))
            frequency = []
            for start, probe in [(0, i/256) for i in range(32)] + [(8*n, target)]:
                counts, ideal = [], []
                for sham in (False, True):
                    weights, vac = g.absorb(probe, start, sham)
                    counts.append(draw(weights, CLOCK_SHOTS, rng))
                    ideal.append(vac)
                frequency.append(dict(start_period=start, probe_energy=probe, counts=counts, exact_vacuum=ideal))
            bank = next(r for r in frequency if r['probe_energy'] == target and r['start_period'] == 0)
            k = bank['counts'][0][0]
            # All six nominated settings plus the six delay holds, simultaneous lower bound.
            ci = float(betaincinv(k, CLOCK_SHOTS-k+1, .01/12)) if k else 0.
            densities = (abs(g.modes[::2, 0])**2 + abs(g.modes[1::2, 0])**2)/2
            panels.append(dict(dimers=n, law=mode, ports=list(map(int, ports)), period=g.period,
                 normal_energies=(g.phi/g.period).tolist(), parent_spectral_weights=(abs(g.source)**2).tolist(),
                 lowest_mode_dimer_weights=densities.tolist(), independent_nominal_target=target,
                 bright_cutoff=threshold, entire_bright_only_tag_upper=g.bright_bound(target, threshold),
                 nominated_probability_lower_99=ci, local_end_banks=local, absorption_banks=frequency))
    return dict(seed=90127, local_shots=LOCAL_SHOTS, absorption_shots=CLOCK_SHOTS,
                single_dimer_clock=cal, clock_shots_per_quadrature=prior.SHOTS, native_constants=[j, h],
                calibrated_single_port_words=source_tag, launch_order=['source','sham'],
                joint_end_order=['parent ground, far ground','parent ground, far excited',
                                 'parent excited, far ground','parent excited, far excited'],
                filter_order=['aperture and vacuum tag','aperture accepted and vacuum tag refused','aperture refused'],
                panels=panels)


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--output', type=Path, required=True)
    args = cli.parse_args()
    args.output.write_text(json.dumps(acquisition(), indent=1, allow_nan=False) + '\n')
