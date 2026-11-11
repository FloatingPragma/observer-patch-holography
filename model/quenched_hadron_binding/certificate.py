"""Receipt of the package: the fast gates of the code and the analysis of the recorded production tables."""
import hashlib
import json
from pathlib import Path

import numpy as np

from . import analysis as an
from . import gates as gt
from .literature import CONVENTION, LITERATURE, MEASURED

TABLES = Path(__file__).with_name('recorded_tables.json')
HEADLINE = 'B'


def recorded():
    return json.loads(TABLES.read_text())


def resolved(mass):
    """A mass [value, statistical error, systematic error] is resolved when it exceeds five combined errors."""
    return bool(mass[0] > 5 * np.hypot(mass[1], mass[2]))


def production(tables):
    """Analysis of every lattice of the recorded tables."""
    out = {}
    for label, dump in sorted(tables['lattices'].items()):
        res = an.analyse(dump)
        res['summary'] = dump['summary']
        res['windows'] = dump['windows']
        used = res['chiral']['kappas_used']
        res['nucleon_not_resolved_at'] = [row['kappa'] for row in res['hadrons']
                                          if row['kappa'] in used and not resolved(row['mass_nucleon'])]
        out[label] = res
    return out


def production_checks(prod):
    checks = {}
    for label, res in prod.items():
        p = res['potential']
        r = np.array(p['r'])
        v = np.array([np.nan if x is None else x for x in p['V_t2']])
        e = np.array([np.nan if x is None else x for x in p['V_t2_err']])
        on_axis = [k for k in range(len(r)) if abs(r[k] - round(r[k])) < 1e-9 and np.isfinite(v[k])]
        rise = (len(on_axis) >= 2 and len(p['on_axis_increments']) == len(on_axis) - 1
                and all(row['value'] > 3 * row['stat'] for row in p['on_axis_increments']))
        checks[f'5_potential.{label}.rises_at_every_on_axis_step'] = bool(rise)
        joint = res['chiral']['joint']
        cov = np.array(joint['covariance'])
        checks[f'8_covariance.{label}.paired_matrix_positive_semidefinite'] = bool(
            np.allclose(cov, cov.T, atol=1e-14) and np.linalg.eigvalsh(cov).min() > -1e-12)
        checks[f'8_covariance.{label}.shared_scale_cancels_in_mass_ratio'] = bool(
            joint['scale_cancelled_ratio_max_error'] < 1e-12)
        checks[f'8_covariance.{label}.contrasts_match_joint_variance'] = all(
            abs(c['stat'] ** 2 - c['variance_from_covariance']) < 1e-12 for c in joint['contrasts'].values())

        checks[f'5_potential.{label}.string_tension_five_sigma'] = bool(p['sigma_over_err'] > 5)
        checks[f'5_potential.{label}.every_fit_has_positive_string_tension'] = bool(
            all(f['sigma'] > 5 * f['sigma_err'] for f in p['fits'].values()))
        pl = res['summary']['plaquette']
        checks[f'5_gauge.{label}.plaquette_matches_literature'] = bool(
            abs(pl['value'] - LITERATURE['plaquette_beta_5.7']['value']) < 5e-4 + 3 * pl['err'])
        ac = res['summary']['autocorrelation']
        checks[f'5_gauge.{label}.configurations_decorrelated'] = bool(ac['plaquette_tau_int_sweeps'] < 10
                                                                     and abs(ac['lag1_plaquette'][0]) < 3 * ac['lag1_plaquette'][1])
        checks[f'6_solver.{label}.true_residual_below_1e-9'] = bool(
            all(s['true_residual_max'] < 1e-9 for s in res['summary']['solver']))
        checks[f'6_solver.{label}.no_configuration_excluded'] = bool(res['summary']['excluded'] == [])
        used = res['chiral']['kappas_used']
        for row in res['hadrons']:
            if row['kappa'] not in used:
                continue
            k = f"{row['kappa']:.4f}"
            mpi, mrho, mn = row['mass_pion'], row['mass_rho'], row['mass_nucleon']
            checks[f'7_hadrons.{label}.kappa_{k}.pion_and_rho_resolved'] = resolved(mpi) and resolved(mrho)
            checks[f'7_hadrons.{label}.kappa_{k}.pion_below_rho'] = bool(mpi[0] < mrho[0])
            if resolved(mn):
                checks[f'7_hadrons.{label}.kappa_{k}.nucleon_resolved_and_above_rho'] = bool(mrho[0] < mn[0])
        checks[f'7_hadrons.{label}.nucleon_not_resolved_at_most_at_the_lightest_quark_mass'] = bool(
            set(res['nucleon_not_resolved_at']) <= {max(used)})
    if HEADLINE in prod:
        res = prod[HEADLINE]
        c, sy = res['chiral']['primary'], res['chiral']['systematic']
        checks[f'8_chiral.{HEADLINE}.four_hopping_parameters'] = bool(len(res['chiral']['kappas_used']) == 4)
        checks[f'8_chiral.{HEADLINE}.rho_and_nucleon_massive_in_chiral_limit'] = bool(
            c['mrho_chiral'][0] > 5 * np.hypot(c['mrho_chiral'][1], sy['mrho_chiral'])
            and c['mN_chiral'][0] > 5 * np.hypot(c['mN_chiral'][1], sy['mN_chiral']))
        checks[f'8_chiral.{HEADLINE}.nucleon_heavier_than_rho'] = bool(
            c['mN_over_mrho_chiral'][0] - 1 > 3 * np.hypot(c['mN_over_mrho_chiral'][1], sy['mN_over_mrho_chiral']))
        checks[f'8_chiral.{HEADLINE}.kappa_c_above_every_hopping_parameter'] = bool(
            c['kappa_c'][0] - 5 * np.hypot(c['kappa_c'][1], sy['kappa_c']) > max(res['chiral']['kappas_used']))
    if 'A' in prod and 'B' in prod:
        a, b = prod['A']['potential'], prod['B']['potential']
        d = abs(a['sigma'] - b['sigma'])
        s = np.hypot(np.hypot(a['sigma_stat'], a['sigma_syst']), np.hypot(b['sigma_stat'], b['sigma_syst']))
        checks['9_volume.string_tension_agrees_between_lattices'] = bool(d < 3 * s)
    return checks


def weighted_slope(r, f, e, r_min=2.0):
    """Slope of the weighted straight line through the points with r >= r_min."""
    r, f, e = (np.asarray(x, float) for x in (r, f, e))
    sel = r >= r_min
    w = 1.0 / e[sel] ** 2
    x, y = r[sel], f[sel]
    s0, s1, s2 = w.sum(), (w * x).sum(), (w * x * x).sum()
    t0, t1 = (w * y).sum(), (w * x * y).sum()
    return float((s0 * t1 - s1 * t0) / (s0 * s2 - s1 * s1))


def melting(tables):
    """Slope of F/T = -ln C(r) below and above the transition, refitted from the recorded points."""
    out = {}
    for tag, q in tables['melting'].items():
        out[tag] = dict(beta=q['beta'], lattice=q['lattice'], streams=q['streams'], measurements=q['measurements'],
                        abs_L=[q['abs_L'], q['abs_L_err']], r_max_resolved=q['r_max_resolved'], points=len(q['r']),
                        slope_refit=weighted_slope(q['r'], q['F_over_T'], q['F_over_T_err']),
                        slope=q['slope_F_over_T_r_ge_2'], on_axis=q['on_axis'])
    return out


def melting_checks(melt):
    lo, hi = melt['below'], melt['above']
    checks = {}
    for tag, q in melt.items():
        checks[f'10_melting.{tag}.recorded_slope_equals_refit'] = bool(abs(q['slope_refit'] - q['slope'][0]) < 1e-9)
    f, e = np.array(lo['on_axis']['F_over_T']), np.array(lo['on_axis']['err'])
    steps = [f[k + 1] - f[k] > 3 * np.hypot(e[k], e[k + 1]) for k in range(5)]
    checks['10_melting.below.free_energy_rises_at_every_on_axis_step_to_r_6'] = bool(all(steps))
    checks['10_melting.below.slope_twenty_sigma'] = bool(lo['slope'][0] > 20 * lo['slope'][1])
    f = np.array(hi['on_axis']['F_over_T'])
    checks['10_melting.above.free_energy_flat_beyond_r_4'] = bool(np.ptp(f[3:]) < 0.01)
    checks['10_melting.slope_above_below_one_per_cent_of_slope_below'] = bool(hi['slope'][0] < 0.01 * lo['slope'][0])
    checks['10_melting.loop_small_below_large_above'] = bool(hi['abs_L'][0] > 10 * lo['abs_L'][0])
    return checks


def headline(prod, melt):
    """The numbers of the README table, from the largest lattice."""
    res = prod[HEADLINE]
    c, sy, p = res['chiral']['primary'], res['chiral']['systematic'], res['potential']

    def tot(k):
        return dict(value=c[k][0], stat=c[k][1], syst=sy[k], err=float(np.hypot(c[k][1], sy[k])))
    return dict(lattice=res['lattice'], beta=res['beta'], configurations=res['configurations'], streams=res['streams'],
                a2_sigma=dict(value=p['sigma'], stat=p['sigma_stat'], syst=p['sigma_syst'],
                              err=float(np.hypot(p['sigma_stat'], p['sigma_syst']))),
                r0_over_a=dict(value=p['r0'], stat=p['r0_stat'], syst=p['r0_syst'], err=float(np.hypot(p['r0_stat'], p['r0_syst']))),
                kappa_c=tot('kappa_c'), a_m_rho_chiral=tot('mrho_chiral'), a_m_N_chiral=tot('mN_chiral'),
                mN_over_mrho_chiral=tot('mN_over_mrho_chiral'), mN_over_mrho_phys=tot('mN_over_mrho_phys'),
                mrho_over_sqrtsigma_chiral=tot('mrho_over_sqrtsigma_chiral'), mN_over_sqrtsigma_chiral=tot('mN_over_sqrtsigma_chiral'),
                mrho_over_sqrtsigma_phys=tot('mrho_over_sqrtsigma_phys'), mN_over_sqrtsigma_phys=tot('mN_over_sqrtsigma_phys'),
                mrho_r0_chiral=tot('mrho_r0_chiral'), mN_r0_chiral=tot('mN_r0_chiral'),
                kappas_used=res['chiral']['kappas_used'], nucleon_not_resolved_at=res['nucleon_not_resolved_at'],
                melting_slope_below=melt['below']['slope'], melting_slope_above=melt['above']['slope'],
                mN_over_mrho_measured_not_blind=MEASURED['mN_over_mrho'])


def certificate():
    tables = recorded()
    checks, numbers = gt.gates()
    prod = production(tables)
    melt = melting(tables)
    checks.update(production_checks(prod))
    checks.update(melting_checks(melt))
    gates = tables['production_gates']
    checks['11_recorded_gates.every_gate_of_the_compiled_kernels_passed'] = bool(all(g['ok'] for g in gates.values()))
    return dict(schema='oph.quenched_hadron_binding.v1',
                scope=('Quenched SU(3) lattice gauge theory with Wilson quarks at one coupling: static potential, pion, '
                       'rho and nucleon masses at four quark masses, chiral limit, and the Polyakov-loop correlator '
                       'below and above the deconfinement transition. Quark loops are dropped. Masses are in lattice '
                       'units and in units of the string tension of the same lattice. No mass in GeV is derived.'),
                all_checks_pass=bool(all(checks.values())), checks=checks, gate_numbers=numbers, production=prod,
                melting=melt, production_gates=gates, headline=headline(prod, melt), measured=MEASURED,
                convention=CONVENTION, literature=LITERATURE,
                input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                              [TABLES, *sorted(TABLES.parent.glob('*.py'))]})


def close(a, b, rtol=1e-6, atol=1e-9, path=''):
    """Recursive comparison of two receipts: exact for booleans, strings and integers, tolerant for floats."""
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return [f'{path}: keys differ']
        return [x for k in a for x in close(a[k], b[k], rtol, atol, f'{path}/{k}')]
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f'{path}: lengths differ']
        return [x for i, (u, v) in enumerate(zip(a, b)) for x in close(u, v, rtol, atol, f'{path}[{i}]')]
    if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, str) or a is None or b is None:
        return [] if a == b else [f'{path}: {a!r} != {b!r}']
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return [] if abs(a - b) <= atol + rtol * max(abs(a), abs(b)) else [f'{path}: {a!r} != {b!r}']
    return [] if a == b else [f'{path}: {a!r} != {b!r}']
