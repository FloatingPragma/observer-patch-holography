"""Fit chain from jackknife samples to the potential, the hadron masses and the ratios.

Input: one entry per lattice of recorded_tables.json with the delete-one-stream jackknife samples of the smeared
Wilson loops W(r, t) and of the folded pion, rho and nucleon correlators. Every number of the result is recomputed
on every jackknife sample. Fit windows and forms are those fixed before the runs (PROOF.md, section 5).
"""
import numpy as np

KINDS = dict(pion='cosh', rho='cosh', nucleon='exp')


def jk_err(samples):
    samples = np.asarray(samples)
    n = samples.shape[0]
    return np.sqrt((n - 1) / n * ((samples - samples.mean(0)) ** 2).sum(0))


def jk_cov(samples):
    """Covariance of jointly deleted streams; the first axis identifies the same omitted stream."""
    a = np.asarray(samples, dtype=float)
    if a.ndim != 2 or len(a) < 2 or not np.isfinite(a).all():
        raise ValueError('joint jackknife needs at least two finite vector replicas')
    centered = a - a.mean(axis=0)
    return (len(a) - 1) / len(a) * centered.T @ centered


def joint_summary(central, replicas):
    """Joint scale/mass covariance and contrasts without independent-error assumptions."""
    names = ['mrho_chiral', 'mN_chiral', 'sqrt_sigma', 'r0',
             'mrho_over_sqrtsigma_chiral', 'mN_over_sqrtsigma_chiral',
             'mrho_r0_chiral', 'mN_r0_chiral', 'mN_over_mrho_chiral']
    samples = np.array([[row[key] for key in names] for row in replicas])
    covariance = jk_cov(samples)
    sd = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(sd, sd)
    contrasts = {}
    for label, rho, nucleon in [('lattice', 'mrho_chiral', 'mN_chiral'),
                               ('over_sqrt_sigma', 'mrho_over_sqrtsigma_chiral', 'mN_over_sqrtsigma_chiral'),
                               ('times_r0', 'mrho_r0_chiral', 'mN_r0_chiral')]:
        i, j = names.index(rho), names.index(nucleon)
        variance = covariance[i, i] + covariance[j, j] - 2 * covariance[i, j]
        contrasts[label] = dict(value=float(central[nucleon] - central[rho]),
                               stat=float(jk_err(samples[:, j] - samples[:, i])),
                               variance_from_covariance=float(variance),
                               independent_stat_control=float(np.hypot(sd[i], sd[j])))
    return dict(names=names, replicas=samples.tolist(), covariance=covariance.tolist(),
                correlation=correlation.tolist(), contrasts=contrasts,
                scale_cancelled_ratio_max_error=float(np.max(abs(
                    samples[:, 5] / samples[:, 4] - samples[:, 8]))),
                sampling_unit='delete-one-stream replicas, paired in their retained order')


def _golden(f, lo, hi, tol=1e-13):
    """Minimum of a unimodal function on [lo, hi] by golden-section search."""
    g = 0.5 * (np.sqrt(5.0) - 1.0)
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return 0.5 * (a + b)


def fit_mass(y, s, t, big_t, kind, steps=60):
    """Uncorrelated least squares of A [exp(-m t) + exp(-m (T - t))] (cosh) or A exp(-m t) (exp).

    The amplitude is eliminated. A grid brackets the minimum of chi2(m), and bisection on the derivative of chi2
    locates it, so the result does not depend on comparisons of nearly equal chi2 values. Returns m and chi2(m).
    """
    w = 1.0 / s ** 2
    back = big_t - t

    def shape(m):
        f = np.exp(-m * t)
        df = -t * f
        if kind == 'cosh':
            b = np.exp(-m * back)
            f, df = f + b, df - back * b
        return f, df

    def chi2(m):
        f, _ = shape(m)
        a = (w * y * f).sum() / (w * f * f).sum()
        return float((w * (y - a * f) ** 2).sum())

    def slope(m):
        f, df = shape(m)
        a = (w * y * f).sum() / (w * f * f).sum()
        return float(-a * (w * (y - a * f) * df).sum())
    grid = np.linspace(0.02, 3.0, 150)
    f = np.exp(-grid[:, None] * t[None, :])
    if kind == 'cosh':
        f = f + np.exp(-grid[:, None] * back[None, :])
    a = (w * y * f).sum(1) / (w * f * f).sum(1)
    k = int(np.argmin((w * (y - a[:, None] * f) ** 2).sum(1)))
    lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
    if slope(lo) < 0.0 < slope(hi):
        for _ in range(steps):
            mid = 0.5 * (lo + hi)
            if slope(mid) < 0.0:
                lo = mid
            else:
                hi = mid
        m = 0.5 * (lo + hi)
    else:
        m = _golden(chi2, lo, hi)
    return float(m), chi2(m)


def effective_mass(c, big_t, kind):
    out = np.full(big_t // 2, np.nan)
    for t in range(big_t // 2):
        a, b = c[t], c[t + 1]
        if not (a > 0 and b > 0):
            continue
        if kind == 'exp':
            out[t] = np.log(a / b)
            continue
        f = lambda m: np.cosh(m * (t - big_t / 2)) / np.cosh(m * (t + 1 - big_t / 2)) - a / b
        lo, hi = 1e-6, 8.0
        if f(lo) * f(hi) > 0:
            continue
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(lo) * f(mid) <= 0:
                hi = mid
            else:
                lo = mid
        out[t] = 0.5 * (lo + hi)
    return out


def linfit(x, y, s):
    w = 1.0 / s ** 2
    s0, s1, s2 = w.sum(), (w * x).sum(), (w * x * x).sum()
    t0, t1 = (w * y).sum(), (w * x * y).sum()
    det = s0 * s2 - s1 * s1
    c0 = (s2 * t0 - s1 * t1) / det
    c1 = (s0 * t1 - s1 * t0) / det
    return c0, c1, float((w * (y - c0 - c1 * x) ** 2).sum())


def potential_fit(r, v, s, fix_e=None):
    """V = V0 + sigma r - e / r by weighted linear least squares. Returns V0, sigma, e, chi2."""
    w = 1.0 / s
    if fix_e is None:
        a = np.stack([np.ones_like(r), r, -1.0 / r], 1) * w[:, None]
        p, *_ = np.linalg.lstsq(a, v * w, rcond=None)
        return p[0], p[1], p[2], float((((v - p[0] - p[1] * r + p[2] / r) * w) ** 2).sum())
    a = np.stack([np.ones_like(r), r], 1) * w[:, None]
    p, *_ = np.linalg.lstsq(a, (v + fix_e / r) * w, rcond=None)
    return p[0], p[1], fix_e, float((((v - p[0] - p[1] * r + fix_e / r) * w) ** 2).sum())


def potential(dump):
    r = np.array(dump['r'])
    wm = np.array(dump['W_mean'])
    wj = np.array(dump['W_jk'])
    out = dict(r=r.tolist(), separations=dump['separations'])
    est = {}
    for t in (1, 2, 3, 4):
        with np.errstate(invalid='ignore', divide='ignore'):
            vm = np.log(wm[:, t - 1] / wm[:, t])
            vj = np.log(wj[:, :, t - 1] / wj[:, :, t])
        est[t] = (vm, vj, jk_err(vj))
        out[f'V_t{t}'] = [None if not np.isfinite(x) else float(x) for x in vm]
        out[f'V_t{t}_err'] = [None if not np.isfinite(x) else float(x) for x in est[t][2]]
    fits = {}
    keep = {}

    def run(tag, t, sel, fix_e=None):
        vm, vj, ve = est[t]
        p = potential_fit(r[sel], vm[sel], ve[sel], fix_e)
        pj = np.array([potential_fit(r[sel], vj[i][sel], ve[sel], fix_e)[:3] for i in range(vj.shape[0])])
        pe = jk_err(pj)
        with np.errstate(invalid='ignore'):
            r0 = np.sqrt((1.65 - p[2]) / p[1])
            r0j = np.sqrt((1.65 - pj[:, 2]) / pj[:, 1])
        fits[tag] = dict(estimator_t=t, r_used=r[sel].tolist(), V0=float(p[0]), V0_err=float(pe[0]), sigma=float(p[1]),
                         sigma_err=float(pe[1]), e=float(p[2]), e_err=float(pe[2]), chi2=p[3],
                         dof=int(sel.sum() - (3 if fix_e is None else 2)), r0=float(r0), r0_err=float(jk_err(r0j)))
        keep[tag] = (p, pj, r0, r0j)
    ok2 = np.isfinite(est[2][0]) & np.isfinite(est[2][2])
    run('F1_primary', 2, ok2 & (r >= 1.41))
    run('F2_all_r', 2, ok2 & (r >= 0.99))
    run('F3_e_fixed_pi12_r_ge_2', 2, ok2 & (r >= 1.99), fix_e=np.pi / 12)
    with np.errstate(invalid='ignore'):
        ok3 = np.isfinite(est[3][0]) & np.isfinite(est[3][2]) & (r >= 1.41) & (est[3][2] < 0.10 * np.abs(est[3][0]))
    if ok3.sum() >= 4:
        run('F4_t3', 3, ok3)
    sig = [f['sigma'] for f in fits.values()]
    r0s = [f['r0'] for f in fits.values()]
    out.update(fits=fits, sigma=fits['F1_primary']['sigma'], sigma_stat=fits['F1_primary']['sigma_err'],
               sigma_syst=float(0.5 * (max(sig) - min(sig))), r0=fits['F1_primary']['r0'],
               r0_stat=fits['F1_primary']['r0_err'], r0_syst=float(0.5 * (np.nanmax(r0s) - np.nanmin(r0s))))
    out['sigma_over_err'] = float(out['sigma'] / np.hypot(out['sigma_stat'], out['sigma_syst']))
    vm, vj, ve = est[2]
    axis = [i for i, x in enumerate(r) if abs(x - round(x)) < 1e-9 and np.isfinite(vm[i])]
    out['on_axis_increments'] = [dict(r_from=float(r[a]), r_to=float(r[b]),
        value=float(vm[b] - vm[a]), stat=float(jk_err(vj[:, b] - vj[:, a])),
        independent_stat_control=float(np.hypot(ve[a], ve[b]))) for a, b in zip(axis[:-1], axis[1:])]

    p, pj, r0, r0j = keep['F1_primary']
    return out, (p[1], pj[:, 1], r0, r0j)


def chiral_chain(mpi, mrho, mn, epi2, erho, en, inv_kappa, sigma, r0, ratio_phys):
    c0, c1, chi_pi = linfit(inv_kappa, mpi ** 2, epi2)
    inv_kc = -c0 / c1
    x = mpi ** 2
    a0, a1, chi_r = linfit(x, mrho, erho)
    n0, n1, chi_n = linfit(x, mn, en)
    c = ratio_phys
    qa, qb, qc = c * c * a1 * a1, 2 * c * c * a0 * a1 - 1.0, c * c * a0 * a0
    xph = (-qb - np.sqrt(qb * qb - 4 * qa * qc)) / (2 * qa) if qa != 0 else -qc / qb
    mrho_ph, mn_ph = a0 + a1 * xph, n0 + n1 * xph
    ss = np.sqrt(sigma)
    return dict(kappa_c=1.0 / inv_kc, slope_mpi2=c1, chi2_pi=chi_pi, mrho_chiral=a0, mrho_slope=a1, chi2_rho=chi_r,
                mN_chiral=n0, mN_slope=n1, chi2_N=chi_n, mpi2_phys=xph, kappa_phys=1.0 / (inv_kc + xph / c1),
                mrho_phys=mrho_ph, mN_phys=mn_ph, mN_over_mrho_chiral=n0 / a0, mN_over_mrho_phys=mn_ph / mrho_ph,
                mrho_over_sqrtsigma_chiral=a0 / ss, mN_over_sqrtsigma_chiral=n0 / ss,
                mrho_over_sqrtsigma_phys=mrho_ph / ss, mN_over_sqrtsigma_phys=mn_ph / ss,
                mrho_r0_chiral=a0 * r0, mN_r0_chiral=n0 * r0, mrho_r0_phys=mrho_ph * r0, mN_r0_phys=mn_ph * r0,
                sqrt_sigma=ss, r0=r0)


def chiral_systematics(variants, potential_fits):
    """Maximum central-value displacement over the crossed hadron and potential fit choices.

    Every hadron variant carries its primary potential scale. Rescaling its normalized masses to each
    potential fit keeps the two fit choices independent. These finite envelopes are sensitivity estimates,
    not confidence intervals or bounds on quenching, finite-volume or continuum effects.
    """
    primary = variants['primary']
    out = {k: 0.0 for k in primary}
    for variant in variants.values():
        for fit in potential_fits.values():
            sigma, r0 = fit['sigma'], fit['r0']
            if not (np.isfinite(sigma) and sigma > 0 and np.isfinite(r0) and r0 > 0):
                raise ValueError('potential fit must have a finite positive string tension and Sommer scale')
            for key, (value, _) in variant.items():
                if '_over_sqrtsigma_' in key:
                    value *= variant['sqrt_sigma'][0] / np.sqrt(sigma)
                elif '_r0_' in key:
                    value *= r0 / variant['r0'][0]
                elif key == 'sqrt_sigma':
                    value = np.sqrt(sigma)
                elif key == 'r0':
                    value = r0
                out[key] = max(out[key], float(abs(value - primary[key][0])))
    return out


def analyse(dump):
    """Full chain for one lattice. Returns a dict with the keys potential, hadrons, chiral."""
    big_t = dump['lattice'][3]
    kappas = dump['kappas']
    win = dump['windows']
    out = dict(label=dump['label'], lattice=dump['lattice'], beta=dump['beta'], configurations=dump['configurations'],
               streams=dump['streams'])
    out['potential'], (sigma_m, sigma_j, r0_m, r0_j) = potential(dump)
    masses = {}
    rows = []
    for ik, kap in enumerate(kappas):
        row = dict(kappa=kap)
        for ch in KINDS:
            for ty in ('SP', 'SS'):
                d = dump['correlators'][f'{ch}_{ty}_{ik}']
                cm, cj = np.array(d['mean']), np.array(d['jk'])
                ce = jk_err(cj)
                if ty == 'SP':
                    row[f'meff_{ch}'] = [None if not np.isfinite(v) else float(v) for v in effective_mass(cm, big_t, KINDS[ch])]
                    row[f'correlator_positive_{ch}'] = bool((cm[win[ch][0]:win[ch][1] + 1] > 0).all())
                for tag, w in (('primary', win[ch]), ('alt', win['alt'][ch])):
                    t = np.arange(w[0], w[1] + 1)
                    m0, c2 = fit_mass(cm[t], ce[t], t, big_t, KINDS[ch])
                    mj = np.array([fit_mass(c[t], ce[t], t, big_t, KINDS[ch])[0] for c in cj])
                    masses[(ch, ty, tag, ik)] = (m0, mj)
                    row[f'mass_{ch}_{ty}_{tag}'] = [m0, float(jk_err(mj)), c2]
        for ch in KINDS:
            m0, mj = masses[(ch, 'SP', 'primary', ik)]
            syst = max(abs(m0 - masses[(ch, 'SP', 'alt', ik)][0]), abs(m0 - masses[(ch, 'SS', 'primary', ik)][0]))
            row[f'mass_{ch}'] = [m0, float(jk_err(mj)), float(syst)]
        for name, a, b in (('mN_over_mrho', 'nucleon', 'rho'), ('mpi_over_mrho', 'pion', 'rho')):
            a0, aj = masses[(a, 'SP', 'primary', ik)]
            b0, bj = masses[(b, 'SP', 'primary', ik)]
            alt = masses[(a, 'SP', 'alt', ik)][0] / masses[(b, 'SP', 'alt', ik)][0]
            ss = masses[(a, 'SS', 'primary', ik)][0] / masses[(b, 'SS', 'primary', ik)][0]
            row[name] = [a0 / b0, float(jk_err(aj / bj)), float(max(abs(a0 / b0 - alt), abs(a0 / b0 - ss)))]
        rows.append(row)
    out['hadrons'] = rows
    inv_k = np.array([1.0 / k for k in kappas])
    use, dropped = [], []
    for i in range(len(kappas)):
        m0, mj = masses[('pion', 'SP', 'primary', i)]
        rel = float(jk_err(mj) / m0) if m0 > 0 else np.inf
        if rows[i]['correlator_positive_pion'] and rel <= 0.25:
            use.append(i)
        else:
            dropped.append(dict(kappa=kappas[i], pion_relative_error=rel, correlator_positive=rows[i]['correlator_positive_pion']))

    primary_joint = {}

    def variant(ty, tag, sel):
        m = {ch: np.array([masses[(ch, ty, tag, i)][0] for i in sel]) for ch in KINDS}
        mj = {ch: np.array([masses[(ch, ty, tag, i)][1] for i in sel]) for ch in KINDS}
        e = {ch: np.array([jk_err(masses[(ch, ty, tag, i)][1]) for i in sel]) for ch in KINDS}
        e2 = np.array([jk_err(masses[('pion', ty, tag, i)][1] ** 2) for i in sel])
        x = inv_k[sel]
        cen = chiral_chain(m['pion'], m['rho'], m['nucleon'], e2, e['rho'], e['nucleon'], x, sigma_m, r0_m,
                           dump['ratio_phys_input'])
        js = [chiral_chain(mj['pion'][:, i], mj['rho'][:, i], mj['nucleon'][:, i], e2, e['rho'], e['nucleon'], x,
                           sigma_j[i], r0_j[i], dump['ratio_phys_input']) for i in range(mj['pion'].shape[1])]
        if ty == 'SP' and tag == 'primary' and sel == use:
            primary_joint.update(joint_summary(cen, js))
        return {k: [float(cen[k]), float(jk_err(np.array([j[k] for j in js])))] for k in cen}
    ch = dict(kappas_used=[kappas[i] for i in use], kappas_dropped=dropped, primary=variant('SP', 'primary', use),
              alt_window=variant('SP', 'alt', use), smeared_sink=variant('SS', 'primary', use))
    names = ['alt_window', 'smeared_sink']
    if len(use) >= 4:
        ch['lightest_three'] = variant('SP', 'primary', use[1:])
        names.append('lightest_three')
    variants = {name: ch[name] for name in ['primary', *names]}
    ch['systematic_hadron_only'] = {
        k: float(max(abs(ch['primary'][k][0] - ch[v][k][0]) for v in names)) for k in ch['primary']}
    ch['systematic_potential_only'] = chiral_systematics({'primary': ch['primary']}, out['potential']['fits'])
    ch['systematic'] = chiral_systematics(variants, out['potential']['fits'])
    ch['joint'] = primary_joint
    out['chiral'] = ch
    return out
