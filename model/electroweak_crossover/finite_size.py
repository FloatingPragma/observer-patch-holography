"""Matched physical-size cutoff fits and finite-volume response from retained lattice sets."""
from math import log, sqrt
import numpy as np


def fit_cutoff(points, powers=(1,)):
    """Fixed physical size, weighted fit in explicit powers of 1/beta_G; no infinite-volume limit."""
    if len(points) < len(powers) + 1 or len({round(p['size'], 8) for p in points}) != 1:
        raise ValueError('a cutoff fit needs enough sets at one matched physical size')
    beta = np.array([p['beta_G'] for p in points])
    if len(set(beta)) != len(beta) or np.any(beta <= 0):
        raise ValueError('distinct positive lattice couplings required')
    errors = np.array([p['err'] for p in points])
    if not np.isfinite(errors).all() or np.any(errors <= 0):
        raise ValueError('finite positive uncertainties required')
    values = np.array([p['value'] for p in points])
    design = np.column_stack([np.ones(len(beta))] + [beta**(-power) for power in powers])
    a = design/errors[:, None]
    if np.linalg.matrix_rank(a) != a.shape[1]:
        raise ValueError('cutoff design is rank deficient')
    cov = np.linalg.inv(a.T @ a)
    coefs = np.linalg.lstsq(a, values/errors, rcond=None)[0]
    residual = values-design@coefs
    return dict(value=float(coefs[0]), stat=float(sqrt(cov[0, 0])), coefficients=coefs.tolist(),
                covariance=cov.tolist(), chi2=float(sum((residual/errors)**2)),
                dof=len(points)-len(coefs), powers=list(powers), beta_G=beta.tolist(), size=points[0]['size'])


def effective_exponent(left, right):
    """Log susceptibility ratio over log spatial-side ratio, for separately sampled volumes."""
    n1, c1, e1 = left; n2, c2, e2 = right
    if not (n2 > n1 > 0 and c1 > 0 and c2 > 0 and e1 >= 0 and e2 >= 0):
        raise ValueError('positive susceptibilities and ordered volume sizes required')
    denom = log(n2/n1)
    return dict(L=[n1, n2], value=log(c2/c1)/denom,
                stat=sqrt((e1/c1)**2+(e2/c2)**2)/denom,
                assumption='independently sampled volume ensembles; delta-method statistical error')


def analyse(production):
    sets = {key:s for key,s in production['chain_point'].items() if s['yscan']['x'] == 0.296}
    result = dict(scope='Finite design sensitivity; no bound on missing higher cutoff powers or infinite volume.')
    result['matched_size_cutoff'] = {}
    for observable in ('y_max', 'T_max'):
        points=[]
        for key,s in sets.items():
            size=4*s['L']/s['bG']
            if size != 8: continue
            q=s['yscan'] if observable=='y_max' else s['tscan']['declared']
            points.append(dict(key=key,size=size,beta_G=s['bG'],value=q[observable],err=q[observable+'_err']))
        points.sort(key=lambda p:p['beta_G'])
        variants={'linear':fit_cutoff(points), 'quadratic_only':fit_cutoff(points,(2,)),
                  'linear_and_quadratic':fit_cutoff(points,(1,2))}
        for point in points:
            variants['omit_'+point['key']]=fit_cutoff([p for p in points if p is not point])
        primary=variants['linear']['value']
        result['matched_size_cutoff'][observable]=dict(points=points,variants=variants,
            sensitivity_envelope=max(abs(v['value']-primary) for v in variants.values()),
            zero_dof_variants=[name for name,v in variants.items() if v['dof']==0])
    groups={}
    for key,s in sets.items():
        if s['bG']==8:
            q=s['yscan'];groups[s['L']]=(s['L'],q['chi3_max'],q['chi3_max_err'])
    ordered=sorted(groups.values())
    result['crossover_volume_exponents']=[effective_exponent(a,b) for a,b in zip(ordered[:-1],ordered[1:])]
    controls=production['validation']['first_order']['sets']
    ordered=sorted((int(L),q['chi3_max'],q['chi3_max_err']) for L,q in controls.items())
    result['first_order_volume_exponents']=[effective_exponent(a,b) for a,b in zip(ordered[:-1],ordered[1:])]
    largest=sorted(sets.values(), key=lambda s:s['L'])
    fixed=[s for s in largest if s['bG']==8]
    result['finite_volume_peak_shifts']=[dict(L=[a['L'],b['L']],
        delta_y=b['yscan']['y_max']-a['yscan']['y_max'],
        stat_y=float(np.hypot(a['yscan']['y_max_err'],b['yscan']['y_max_err'])),
        delta_T_GeV=b['tscan']['declared']['T_max']-a['tscan']['declared']['T_max']) for a,b in zip(fixed[:-1],fixed[1:])]
    return result
