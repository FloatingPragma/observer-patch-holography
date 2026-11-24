"""Exact phase-marginalized galaxy likelihood for RER's regularly spaced cuts."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
import held_galaxies as held
import rar_deep_regime_diagnostic as base
import sparc_full_rar_calibration as catalogue
from lensing_profile import FOUR_G

SPACINGS = (0., .1, .25, .5, 1., 2., 4., 8., 12., 16.)
UNIT = 1e-10


def sample():
    catalogue.verify_source_snapshot()
    meta = {g["name"]: g for g in base.read_table1()}
    galaxies = {}
    for ring in base.read_table2():
        name = ring["name"]
        galaxy = meta[name]
        radius, obs, err = (ring[k] for k in ("rad_kpc", "vobs_kms", "e_vobs_kms"))
        if galaxy["quality"] > 2 or galaxy["inclination_deg"] < 30 or radius <= 0 or err <= 0:
            continue
        parts = np.array([ring[k] for k in ("vgas_kms", "vdisk_kms", "vbul_kms")])
        vb2 = (parts * abs(parts)) @ [1., base.UPSILON_DISK, base.UPSILON_BULGE]
        gbar = vb2 * 1e6 / (base.KPC_M * radius)
        if 0 < gbar < .1 * base.DEEP_CUT_REFERENCE_A0_M_S2:
            galaxies.setdefault(name, []).append((radius, vb2, obs, err))
    return {name: np.array(galaxies[name]) for name in sorted(galaxies)}


def phase_law(radius, spacing):
    if spacing == 0:
        return np.ones(1), np.ones((1, len(radius)))
    x = radius / spacing
    edges = np.unique(np.r_[0., 1 - (x % 1), 1.])
    midpoints = (edges[1:] + edges[:-1]) / 2
    multiplier = np.floor(x[None, :] + midpoints[:, None]) / x
    return np.diff(edges), multiplier


class GalaxyLikelihood:
    def __init__(self, galaxies, spacing, radial_scatter=False):
        self.names = list(galaxies)
        self.spacing = spacing
        self.radial_scatter = radial_scatter
        probabilities, sizes, groups, radial = [], [], [], []
        for name, data in galaxies.items():
            r, vb2, v, err = data.T
            weights, multiplier = phase_law(r, spacing)
            smooth = np.sqrt(vb2 * r * base.KPC_M * UNIT / 1e6)
            probabilities.extend(weights)
            groups.append(len(weights))
            sizes.extend([len(r)] * len(weights))
            radial.append(np.column_stack([np.tile(vb2, len(weights)),
                                           (smooth * multiplier).ravel(),
                                           np.tile(v, len(weights)), np.tile(err, len(weights)),
                                           np.tile(r, len(weights))]))
        self.vb2, self.amplitude, self.obs, self.err, self.radius = np.concatenate(radial).T
        self.log_weights = np.log(probabilities)
        self.phase_starts = np.r_[0, np.cumsum(sizes)[:-1]]
        self.galaxy_starts = np.r_[0, np.cumsum(groups)[:-1]]
        self.groups = groups
        self.sizes = np.array(sizes)

    def log_scores(self, x):
        log_a, within_var, between_var = x[:3]
        if self.radial_scatter:
            within_var = within_var + x[3] / self.radius**2
        mean = np.sqrt(self.vb2 + self.amplitude * np.exp(log_a/2))
        noise = self.err**2 + within_var * mean**2
        residual = self.obs - mean
        def add(values):
            return np.add.reduceat(values, self.phase_starts)
        # Matrix determinant lemma and a rank-one inverse on each common phase.
        denom = 1 + between_var * add(mean**2 / noise)
        quad = add(residual**2 / noise) - between_var * add(mean*residual/noise)**2 / denom
        terms = self.log_weights - .5*(quad + add(np.log(noise)) + np.log(denom)
                                       + self.sizes*np.log(2*np.pi))
        peak = np.maximum.reduceat(terms, self.galaxy_starts)
        return peak + np.log(np.add.reduceat(np.exp(terms - np.repeat(peak, self.groups)),
                                             self.galaxy_starts))

    def fit(self, train):
        def objective(x):
            return -self.log_scores(x)[train].sum()
        bounds = [(np.log(.01), np.log(10.)), (0., 1.), (0., 1.)]
        starts = ([-.2, .004, .03], [0., .016, .008])
        if self.radial_scatter:
            bounds.append((0., 16.))
            starts = [s + [.01] for s in starts]
        runs = [minimize(objective, x, method="L-BFGS-B", bounds=bounds,
                         options={"ftol": 1e-12, "maxiter": 180}) for x in starts]
        viable = [r for r in runs if r.success and np.isfinite(r.fun)]
        if not viable:
            raise RuntimeError([r.message for r in runs])
        return min(viable, key=lambda r: r.fun).x


def parameters(x):
    return dict(a0_m_s2=float(UNIT*np.exp(x[0])), ring_fractional_scatter=float(np.sqrt(x[1])),
                galaxy_fractional_scatter=float(np.sqrt(x[2])))


def external(x, spacing):
    a0 = UNIT*np.exp(x[0])
    rings = json.loads((ROOT/'udg_transfer/rings.json').read_text())['panels']['figure8']['rings']
    data = np.array([(r['radius_kpc'], r['baryonic_kms']**2, r['circular_kms'],
                      r['plotted_error_kms']) for r in rings])
    fit = GalaxyLikelihood({'AGC114905':data}, spacing)
    r, vb2, obs, err = data.T
    weight, ratio = phase_law(r, spacing)
    smooth = np.sqrt(vb2 * r * base.KPC_M * a0 / 1e6)
    means = np.sqrt(vb2 + smooth * ratio)
    deviation = np.sqrt(err**2 + (x[1]+x[2])*means**2)
    udg = dict(full_curve_log_density=float(fit.log_scores(x)[0]), outer_radius_kpc=r[-1],
               outer_speed_kms=obs[-1], outer_error_kms=err[-1],
               smooth_outer_mean_kms=float(np.sqrt(vb2[-1]+smooth[-1])),
               min_phase_outer_mean_kms=float(means[:,-1].min()),
               outer_lower_tail_probability=float(weight @ ndtr((obs[-1]-means[:,-1])/deviation[:,-1])),
               necessary_spacing_kpc_at_outer_plus_two_errors=float(max(
                   0., r[-1]*(1-((obs[-1]+2*err[-1])**2-vb2[-1])/smooth[-1]))))
    profiles = json.loads((ROOT/'receipts/lensing_profile.json').read_text())['profiles']
    shear = {}
    for name in ('KiDS-bright', 'GAMA'):
        profile = profiles[name]
        mask = np.array(profile['lower_bin_edge_m_s2']) >= 1e-13
        g = np.array(profile['g_bar_m_s2'])[mask]
        predicted = (4*g/np.pi + np.sqrt(g*a0))/FOUR_G
        residual = np.array(profile['tangential_esd'])[mask] - predicted
        cov = np.array(profile['covariance'])[np.ix_(mask, mask)]
        shear[name] = dict(bins=len(g), chi_square=float(residual @ np.linalg.solve(cov,residual)),
                           mean_esd=predicted.tolist())
    return dict(AGC114905=udg, shear=shear)


def run():
    galaxies = sample()
    engines = [GalaxyLikelihood(galaxies, s) for s in SPACINGS]
    n = len(galaxies)
    folds = np.array([held.fold_for(name) for name in galaxies])
    out = np.zeros((n,2))
    results, profile = [], []
    for fold in [0, 1, 2, 3, 4, None]:
        train = np.ones(n, bool) if fold is None else folds != fold
        fits = [eng.fit(train) for eng in engines]
        scores = np.array([eng.log_scores(x) for eng,x in zip(engines,fits,strict=True)])
        train_score = scores[:,train].sum(axis=1)
        best = int(np.argmax(train_score))
        if fold is None:
            profile = [dict(spacing_kpc=s, log_density=float(sc), **parameters(x))
                       for s, sc, x in zip(SPACINGS,train_score,fits,strict=True)]
            transfer = {title:external(fits[i], SPACINGS[i]) for title,i in [('smooth',0),('selected',best)]}
            selected = best
        else:
            out[~train] = scores[[0,best]][:,~train].T
            results.append(dict(fold=fold, training=int(train.sum()), held=int((~train).sum()),
                                spacing_kpc=SPACINGS[best], selected=parameters(fits[best]),
                                smooth=parameters(fits[0]), held_log_density_delta=float(
                                    (scores[best,~train]-scores[0,~train]).sum())))
    diff = out[:,1]-out[:,0]
    rng = np.random.default_rng(112419)
    ci = np.quantile(rng.choice(diff,(5000,n)).mean(axis=1),[.025,.975])
    return dict(galaxies=n,rings=sum(map(len,galaxies.values())),
                selected_spacing_kpc=SPACINGS[selected],full_sample_profile=profile, folds=results,
                held_log_density_nats_per_galaxy_delta=float(diff.mean()),
                paired_fixed_prediction_interval=ci.tolist(),
                held_log_density_by_galaxy={k:v.tolist() for k,v in zip(galaxies,out,strict=True)},
                external=transfer)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2,allow_nan=False)+'\n')
    else:
        from lensing_profile import check
        check(json.loads(path.read_text()), result)
    print(result['galaxies'], result['rings'], result['selected_spacing_kpc'],
          result['held_log_density_nats_per_galaxy_delta'])
