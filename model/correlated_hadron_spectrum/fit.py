"""Full-time covariance fits; empirical mass values never choose a fit rule."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from model.quenched_hadron_binding import analysis as old

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TABLE = ROOT / 'model/quenched_hadron_binding/recorded_tables.json'
OLD_RESULT = ROOT / 'model/quenched_hadron_binding/results.json'
RULE = dict(
    lattice='B', covariance='full time covariance from paired stream means; reestimate after each outer deletion',
    rank='full numerical rank; eigenvalue exceeds 1e-12 times largest correlation eigenvalue',
    max_correlation_condition=1e8, minimum_stream_surplus=2,
    singular_policy='reject complete fit window; no pseudoinverse, shrinkage, diagonal fallback or replica deletion',
    windows='both original primary and alternate windows, both SP and SS sinks, all four kappas',
    mass_interval=[0.02, 3.0], grid_points=150,
    minimum_selection='all grid-bracketed local minima and both endpoints; smallest profile objective',
    chain='original pion positivity and <=25percent relative-error admission; all replica labels remain paired',
    model='single cosh for pion/rho, single exponential for nucleon',
    chiral_regression='retain original diagonal mass-error weighted chiral regressions as a separate approximation',
    scale='original paired Wilson-potential fits; no new potential-time correlation model',
    target_mass_used_to_choose_rules=False)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stream_means(mean, jackknife):
    """Invert delete-one-stream means, retaining all stream labels."""
    mean = np.asarray(mean, float)
    jackknife = np.asarray(jackknife, float)
    require(jackknife.ndim == 2 and jackknife.shape[1:] == mean.shape,
            'mean and stream-deletion shape differ')
    require(np.isfinite(mean).all() and np.isfinite(jackknife).all(), 'nonfinite retained data')
    n = len(jackknife)
    scale = max(float(np.max(abs(mean))), np.finfo(float).tiny)
    require(np.max(abs(jackknife.mean(0)-mean)) <= 1e-11*scale,
            'retained means do not obey the delete-one identity')
    # The retained decimal central mean and mean of jackknife rows are
    # rounded independently. A common shift cancels from every covariance.
    samples = n*jackknife.mean(0)-(n-1)*jackknife
    require(np.allclose((samples.sum(0)-samples)/(n-1), jackknife,
                        rtol=1e-11, atol=scale*1e-13), 'stream reconstruction failed')
    return samples


def mean_covariance(samples):
    x = np.asarray(samples, float)
    require(x.ndim == 2 and len(x) >= x.shape[1]+RULE['minimum_stream_surplus'],
            'too few independent streams for the fixed covariance rule')
    d = x-x.mean(0)
    return d.T@d/(len(x)*(len(x)-1))


def covariance_metric(covariance):
    c = np.asarray(covariance, float)
    require(c.ndim == 2 and c.shape[0] == c.shape[1] and np.isfinite(c).all(),
            'finite square covariance required')
    require(np.all(np.diag(c) > 0), 'zero marginal variance is rejected')
    sd = np.sqrt(np.diag(c))
    corr = c/np.outer(sd, sd)
    require(np.max(abs(corr-corr.T)) < 1e-12, 'asymmetric covariance')
    eig, frame = np.linalg.eigh((corr+corr.T)/2)
    require(eig[0] > 1e-12*eig[-1], 'singular covariance rejected')
    condition = float(eig[-1]/eig[0])
    require(condition <= RULE['max_correlation_condition'], 'ill-conditioned covariance rejected')
    inverse = (frame/eig)@frame.T
    return sd, inverse, dict(rank=len(eig), correlation_eigenvalues=eig.tolist(),
                             correlation_condition=condition)


def fit_mass(y, covariance, times, temporal_extent, kind):
    """Profile out amplitude, fit mass with the frozen bounded search."""
    y = np.asarray(y, float)
    times = np.asarray(times, float)
    require(kind in {'cosh', 'exp'} and y.shape == times.shape, 'unsupported fit coordinates')
    sd, precision, diagnostic = covariance_metric(covariance)
    z = y/sd
    t0 = float(times.min())

    def profile(mass, with_amplitude=False):
        shape = np.exp(-mass*(times-t0))
        if kind == 'cosh':
            shape += np.exp(-mass*(temporal_extent-times-t0))
        f = shape/sd
        amplitude = float(f@precision@z/(f@precision@f))
        residual = z-amplitude*f
        chi2 = float(residual@precision@residual)
        return (chi2, float(amplitude*np.exp(mass*t0))) if with_amplitude else chi2

    grid = np.linspace(*RULE['mass_interval'], RULE['grid_points'])
    values = np.array([profile(m) for m in grid])
    candidates = [(float(grid[0]), float(values[0])), (float(grid[-1]), float(values[-1]))]
    brackets = []
    for k in range(1, len(grid)-1):
        if values[k] <= values[k-1] and values[k] <= values[k+1]:
            answer = minimize_scalar(profile, bounds=(grid[k-1], grid[k+1]),
                                     method='bounded', options={'xatol': 1e-13, 'maxiter': 200})
            require(answer.success, 'fixed mass minimizer failed')
            candidates.append((float(answer.x), float(answer.fun)))
            brackets.append([float(grid[k-1]), float(grid[k+1])])
    mass, _ = min(candidates, key=lambda item: (item[1], item[0]))
    chi2, amplitude = profile(mass, True)
    require(np.isfinite(chi2) and np.isfinite(amplitude), 'nonfinite fitted result')
    return dict(mass=mass, amplitude=amplitude, chi2=chi2, dof=len(times)-2,
                boundary_minimum=mass in RULE['mass_interval'], positive_amplitude=amplitude > 0,
                numerical_search='fixed grid brackets, not a global optimization proof',
                local_brackets=brackets, covariance=diagnostic)


def fit_window(mean, jackknife, times, temporal_extent, kind):
    samples = stream_means(mean, jackknife)
    t = np.asarray(times, int)
    full = samples[:, t]
    cov = mean_covariance(full)
    central = fit_mass(np.asarray(mean)[t], cov, t, temporal_extent, kind)
    replicas = []
    for omitted in range(len(samples)):
        reduced = np.delete(full, omitted, axis=0)
        replica = fit_mass(np.asarray(jackknife)[omitted, t], mean_covariance(reduced),
                           t, temporal_extent, kind)
        replicas.append(replica)
    masses = np.array([row['mass'] for row in replicas])
    return dict(central=central, replica_masses=masses.tolist(),
                retained_mean_minus_jackknife_mean=(np.asarray(mean)-np.asarray(jackknife).mean(0)).tolist(),
                replica_amplitudes=[row['amplitude'] for row in replicas],
                replica_chi2=[row['chi2'] for row in replicas],
                replica_covariance_conditions=[row['covariance']['correlation_condition'] for row in replicas],
                boundary_replica_indices=[i for i,row in enumerate(replicas) if row['boundary_minimum']],
                nonpositive_amplitude_replica_indices=[i for i,row in enumerate(replicas) if not row['positive_amplitude']],
                statistical_error=float(old.jk_err(masses)), central_time_covariance=cov.tolist(),
                central_covariance_equals_original_jackknife=bool(np.allclose(cov, old.jk_cov(np.asarray(jackknife)[:,t]), rtol=1e-9, atol=1e-25)),
                replica_count=len(replicas), time_indices=t.tolist())


def gev_row(chain, ratio, rho_anchor, string_anchor):
    rho = float(chain['mrho_phys'])
    nucleon = float(chain['mN_phys'])
    sigma = float(chain['sqrt_sigma'])
    require(rho > 0 and nucleon > 0 and sigma > 0, 'positive fitted scale chain required')
    return dict(pion_rho_anchor=ratio*rho_anchor, rho_rho_anchor=rho_anchor,
                nucleon_rho_anchor=rho_anchor*nucleon/rho,
                pion_string_anchor=ratio*string_anchor*rho/sigma,
                rho_string_anchor=string_anchor*rho/sigma,
                nucleon_string_anchor=string_anchor*nucleon/sigma)


def fit_all(dump):
    fits, rejected = {}, []
    for ik, kappa in enumerate(dump['kappas']):
        for channel, kind in old.KINDS.items():
            for sink in ('SP','SS'):
                data = dump['correlators'][f'{channel}_{sink}_{ik}']
                for tag, window in [('primary',dump['windows'][channel]), ('alt',dump['windows']['alt'][channel])]:
                    key = f'{channel}/{sink}/{tag}/{ik}'
                    try:
                        fits[key] = fit_window(data['mean'], data['jk'],
                            np.arange(window[0],window[1]+1),dump['lattice'][3],kind)
                    except ValueError as error:
                        rejected.append(dict(key=key,kappa=kappa,reason=str(error)))
    return fits, rejected


def mass_chain(dump, fits, sigma, sigma_j, r0, r0_j, sink, window, selected):
    getter = lambda ch,i: fits[f'{ch}/{sink}/{window}/{i}']
    central = {ch:np.array([getter(ch,i)['central']['mass'] for i in selected]) for ch in old.KINDS}
    replicas = {ch:np.array([getter(ch,i)['replica_masses'] for i in selected]) for ch in old.KINDS}
    errors = {ch:np.array([getter(ch,i)['statistical_error'] for i in selected]) for ch in old.KINDS}
    epi2 = np.array([old.jk_err(np.array(getter('pion',i)['replica_masses'])**2) for i in selected])
    inv_kappa = 1/np.array(dump['kappas'])[selected]
    def chain(m,sig,rzero):
        return old.chiral_chain(m['pion'],m['rho'],m['nucleon'],epi2,errors['rho'],errors['nucleon'],
                                inv_kappa,sig,rzero,dump['ratio_phys_input'])
    cen = chain(central,sigma,r0)
    js = [chain({ch:v[:,j] for ch,v in replicas.items()},sigma_j[j],r0_j[j])
          for j in range(dump['streams'])]
    require(all(np.isfinite(v) for row in [cen,*js] for v in row.values()), 'nonfinite chiral chain')
    return cen,js


def source_bindings():
    paths = [TABLE,OLD_RESULT,ROOT/'model/quenched_hadron_binding/analysis.py',
             ROOT/'model/quenched_hadron_binding/PROOF.md',ROOT/'model/pion_mass_readout/results.json']
    paths += [p for p in HERE.iterdir() if p.suffix in {'.py','.md'}]
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def certificate():
    pins = source_bindings()
    tables = json.loads(TABLE.read_text())
    dump = tables['lattices'][RULE['lattice']]
    archived = json.loads(OLD_RESULT.read_text())
    rho_anchor = archived['measured']['m_rho_MeV']/1000
    string_anchor = archived['convention']['sqrt_sigma_MeV']/1000
    fits,rejected = fit_all(dump)
    require(not rejected, 'at least one fixed window rejected; inspect fit_all instead of silently changing rules')
    potential,(sigma,sigma_j,r0,r0_j) = old.potential(dump)
    selected=[]
    for i in range(len(dump['kappas'])):
        f=fits[f'pion/SP/primary/{i}'];w=dump['windows']['pion']
        mean=dump['correlators'][f'pion_SP_{i}']['mean']
        if min(mean[w[0]:w[1]+1])>0 and f['statistical_error']/f['central']['mass']<=.25:selected.append(i)
    require(len(selected)>=3,'insufficient predeclared pion-admitted points')
    choices=[('primary','SP','primary',selected),('alt_window','SP','alt',selected),
             ('smeared_sink','SS','primary',selected)]
    if len(selected)>=4:choices.append(('lightest_three','SP','primary',selected[1:]))
    variants={};raw={};gev={}
    for name,sink,window,sel in choices:
        cen,replicas=mass_chain(dump,fits,sigma,sigma_j,r0,r0_j,sink,window,sel)
        raw[name]=(cen,replicas)
        variants[name]={k:[float(v),float(old.jk_err([row[k] for row in replicas]))]for k,v in cen.items()}
        g=gev_row(cen,dump['ratio_phys_input'],rho_anchor,string_anchor)
        gj=[gev_row(row,dump['ratio_phys_input'],rho_anchor,string_anchor)for row in replicas]
        names=list(g);matrix=np.array([[row[k]for k in names]for row in gj])
        gev[name]=dict(names=names,central=g,replicas=matrix.tolist(),
                       statistical={k:float(old.jk_err(matrix[:,j]))for j,k in enumerate(names)},
                       covariance=old.jk_cov(matrix).tolist())
    sensitivity=old.chiral_systematics(variants,potential['fits'])
    diag=old.analyse(dump)
    primary=gev['primary'];cross=[]
    for name,(cen,_)in raw.items():
        for pname,p in potential['fits'].items():
            changed=dict(cen,sqrt_sigma=float(np.sqrt(p['sigma'])),r0=float(p['r0']))
            cross.append(dict(hadron_choice=name,potential_choice=pname,
                              values=gev_row(changed,dump['ratio_phys_input'],rho_anchor,string_anchor)))
    shifts={k:max(abs(row['values'][k]-v)for row in cross)for k,v in primary['central'].items()}
    baseline=archived['production']['B']['chiral']['primary']
    require(all(np.allclose(diag['chiral']['primary'][k],v,rtol=1e-9,atol=1e-12)for k,v in baseline.items()),
            'independent diagonal control does not replay the archived fit')
    diagonal_gev=gev_row({k:v[0]for k,v in baseline.items()},dump['ratio_phys_input'],rho_anchor,string_anchor)
    mass_comparison=[]
    for i,kappa in enumerate(dump['kappas']):
        channels={}
        for channel in old.KINDS:
            correlated=fits[f'{channel}/SP/primary/{i}']
            control=diag['hadrons'][i]['mass_'+channel]
            channels[channel]=dict(correlated_mass=correlated['central']['mass'],
                correlated_statistical=correlated['statistical_error'],
                diagonal_mass=control[0],diagonal_statistical=control[1],
                central_shift=correlated['central']['mass']-control[0])
        mass_comparison.append(dict(kappa=kappa,channels=channels))
    checks=dict(all_fixed_windows_admitted=len(fits)==48,all_replicas_retained=all(f['replica_count']==32 for f in fits.values()),
        covariance_reconstructed=all(f['central_covariance_equals_original_jackknife']for f in fits.values()),
        paired_scale_chain=True,no_covariance_fallback=True,no_target_selected_rule=True,
        diagonal_control_replayed=True,all_original_four_kappa_rows_retained=True)
    require(all(checks.values()),'correlated fit checks failed')
    result=dict(schema='oph.correlated_hadron_spectrum.v1',source_bindings=pins,rule=RULE,
        all_checks_pass=True,checks=checks,source_lattice={k:dump[k]for k in ['label','lattice','beta','streams','configurations']},
        fits=fits,rejected_windows=rejected,kappas_used=[dump['kappas'][i]for i in selected],
        inherited_unresolved_nucleon_kappas=archived['production']['B']['nucleon_not_resolved_at'],
        chiral_variants=variants,chiral_crossed_fit_sensitivity=sensitivity,
        paired_primary_joint=old.joint_summary(*raw['primary']),gev=gev,
        gev_crossed_fit_choices=cross,gev_fit_sensitivity=shifts,
        empirical_calibrations=dict(light_quark_ratio=dump['ratio_phys_input'],rho_anchor_GeV=rho_anchor,
                                    alternative_string_anchor_GeV=string_anchor,used_only_after_mass_fit=True),
        diagonal_control=dict(chiral_primary=diag['chiral']['primary'],gev=diagonal_gev,
                              all_four_kappa_mass_comparisons=mass_comparison,
                              correlated_minus_diagonal_GeV={k:primary['central'][k]-v for k,v in diagonal_gev.items()}),
        covariance_sampling='finite32-stream plug-in covariance, reestimated inside each outer deletion; no chi-square coverage claim',
        scope=dict(raw_gauge_configurations_available=False,independent_particle_prediction=False,
                   physical_source_clock_attachment=False,sea_quark_volume_continuum_errors_bounded=False,
                   pion_rho_anchored_mass_is_calibration=True))
    require(source_bindings()==pins,'source changed during calculation')
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result=certificate();data=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:(HERE/'results.json').write_text(data)
    if args.check:require((HERE/'results.json').read_text()==data,'correlated hadron receipt differs')
    print(json.dumps(dict(all_checks_pass=result['all_checks_pass'],fits=len(result['fits']),
                         GeV=result['gev']['primary']['central'])))


if __name__=='__main__':main()
