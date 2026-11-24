"""Held small-radius scatter and distance controls for the finite-cut hypothesis."""

import argparse
import json

import numpy as np
from scipy.special import logsumexp, ndtri

import fit


def paired_score(scores):
    rng = np.random.default_rng(112419)
    return dict(mean_nats_per_galaxy=float(np.mean(scores)),
                paired_bootstrap_95=np.quantile(rng.choice(scores,(5000,len(scores))).mean(axis=1),
                                               [.025,.975]).tolist())


def theta(record):
    return np.array([np.log(record['a0_m_s2']/fit.UNIT), record['ring_fractional_scatter']**2,
                     record['galaxy_fractional_scatter']**2])


def at_distances(group, parameters, spacing, relative_error, nodes):
    # Equal-probability quadrature for a common positive distance in the catalogue frame.
    z = ndtri((np.arange(nodes)+.5)/nodes)
    width = np.sqrt(np.log1p(relative_error**2))
    distances = np.exp(width*z-width**2/2)
    galaxies = {}
    for i,d in enumerate(distances):
        datum = group.copy()
        datum[:,:2] *= d
        galaxies[i] = datum
    scores = fit.GalaxyLikelihood(galaxies, spacing).log_scores(parameters)
    return float(logsumexp(scores)-np.log(nodes))


def run():
    bank = json.loads((fit.HERE/'results.json').read_text())
    sample = fit.sample()
    names = list(sample)
    geometry = {g['name']:g['e_dist_mpc']/g['dist_mpc'] for g in fit.base.read_table1()}
    engine = fit.GalaxyLikelihood(sample, 0, radial_scatter=True)
    smooth = np.array([bank['held_log_density_by_galaxy'][k][0] for k in names])
    grid = np.array([bank['held_log_density_by_galaxy'][k][1] for k in names])
    other = np.empty(len(names))
    folds = np.array([fit.held.fold_for(k) for k in names])
    integrated = {n:np.zeros(len(names)) for n in (33,65)}
    rows = []
    for fold in bank['folds']:
        hold = folds == fold['fold']
        x = engine.fit(~hold)
        other[hold] = engine.log_scores(x)[hold]
        rows.append(dict(fold=fold['fold'], scatter_length_kpc=float(np.sqrt(x[3])),
                         **fit.parameters(x)))
        for i in np.flatnonzero(hold):
            name = names[i]
            for nodes, deltas in integrated.items():
                results = [at_distances(sample[name], theta(fold[k]), s, geometry[name], nodes)
                           for k,s in [('selected',fold['spacing_kpc']),('smooth',0.)]]
                deltas[i] = results[0]-results[1]
    return dict(folds=rows,smooth_radial_scatter_minus_smooth=paired_score(other-smooth),
                cuts_minus_smooth_radial_scatter=paired_score(grid-other),
                distance_averaged_cuts_minus_smooth={str(n):paired_score(s) for n,s in integrated.items()},
                largest_individual_33_to_65_node_change_nats=float(max(abs(integrated[33]-integrated[65]))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    output = run()
    dest = fit.HERE/'controls.json'
    if args.write:
        dest.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    else:
        from lensing_profile import check
        check(json.loads(dest.read_text()),output)
    print(output['cuts_minus_smooth_radial_scatter'],output['distance_averaged_cuts_minus_smooth'])
