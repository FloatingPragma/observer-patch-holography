"""All checks of the package: exact centre statements, kernel gates on small lattices, the recorded lattice tables
and the scale arithmetic."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import centre, scale, su3

HERE = Path(__file__).parent
BETA_C_LITERATURE = 5.6925          # Boyd et al. 1996, N_t = 4, infinite volume
PLAQUETTE_LITERATURE = {'5.7': 0.5492, '6.0': 0.5937}   # Lepage and Mackenzie 1993


def stable(obj, digits=11):
    """Round every float to a fixed number of significant digits, so that the receipt is the same on every platform."""
    if isinstance(obj, float):
        return float(f'{obj:.{digits}g}')
    if isinstance(obj, dict):
        return {k: stable(v, digits) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [stable(v, digits) for v in obj]
    return obj


def recorded():
    return json.loads((HERE / 'lattice_tables.json').read_text())


def strong_coupling():
    rows = {}
    for beta in (0.5, 1.0, 2.0, 3.0):
        u = su3.one_plaquette_average(beta)
        rows[str(beta)] = dict(one_plaquette=round(u, 12), series=round(beta / 18 + beta ** 2 / 216 - 5 * beta ** 4 / 93312, 12),
                               with_cube_term=round(u + 4 * u ** 5, 12))
    moments = {str(k): round(su3.haar_moment_re_trace(k), 12) + 0.0 for k in range(1, 6)}
    return dict(haar_moments_re_trace=moments, one_plaquette=rows)


def kernel_gates(seed=11):
    """Monte Carlo gates of the numpy kernel on 4^4. Returns measured numbers; they depend on the random stream."""
    out = {}
    lat = su3.Lattice((4, 4, 4, 4), 5.7, seed=seed, start='hot')
    lat.compound_sweeps(3, 2)
    out['unitarity_defect'] = max(lat.unitarity_defect())
    p = lat.plaquette()
    lat.sweep(False)
    out['overrelaxation_plaquette_change'] = abs(lat.plaquette() - p)
    lat = su3.Lattice((4, 4, 4, 4), 1.0, seed=seed + 1, start='hot')
    lat.compound_sweeps(10, 1)
    ps = []
    for _ in range(60):
        lat.compound_sweeps(1, 1)
        ps.append(lat.plaquette())
    u = su3.one_plaquette_average(1.0)
    out['strong_coupling_beta_1'] = dict(measured=float(np.mean(ps)), err=float(np.std(ps, ddof=1) / math.sqrt(len(ps))),
                                         expected=u + 4 * u ** 5)
    return out


def register_growth(reg):
    """Weighted fit of ln C_max against ln L over the sizes whose runs pass between the branches."""
    out = {}
    for name in ('Z3', 'Z2'):
        sizes = [s for s in reg[name]['sizes'] if not s.get('metastable')]
        lx = np.log([s['L'] for s in sizes])
        ly = np.log([s['C_max'] for s in sizes])
        w = np.array([(s['C_max'] / s['C_max_err']) ** 2 for s in sizes])
        A = np.vstack([lx, np.ones_like(lx)]).T
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        sol = cov @ (A.T @ (w * ly))
        dof = len(sizes) - 2
        chi2 = float((w * (A @ sol - ly) ** 2).sum())
        out[name] = dict(L=[s['L'] for s in sizes], exponent=float(sol[0]),
                         err=float(math.sqrt(cov[0, 0]) * math.sqrt(max(chi2 / dof, 1.0))), chi2=chi2, dof=dof)
    return out


def table_checks(rec):
    """Internal consistency of the recorded lattice tables and their comparison with the imported values."""
    c = {}
    lat = rec['su3']
    g = rec['gates']
    c['gates_all_pass'] = all(v['ok'] for v in g.values() if isinstance(v, dict))
    for b, ref in PLAQUETTE_LITERATURE.items():
        row = rec['plaquette_symmetric'][b]
        c[f'plaquette_beta_{b}_within_tolerance'] = abs(row['value'] - ref) < 5e-4 + 3 * row['err']
    sizes = lat['sizes']
    ns = [s['N_s'] for s in sizes]
    c['three_volumes'] = ns == [8, 12, 16]
    c['chi_max_grows_with_volume'] = all(sizes[i + 1]['chi_max'] > sizes[i]['chi_max'] + 3 * math.hypot(
        sizes[i + 1]['chi_max_err'], sizes[i]['chi_max_err']) for i in range(len(sizes) - 1))
    # refit of the volume exponent and of the infinite-volume coupling from the recorded per-size values
    lx = np.log([s['N_s'] for s in sizes])
    ly = np.log([s['chi_max'] for s in sizes])
    w = np.array([(s['chi_max'] / s['chi_max_err']) ** 2 for s in sizes])
    A = np.vstack([lx, np.ones_like(lx)]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    sol = cov @ (A.T @ (w * ly))
    c['volume_exponent_refit'] = abs(sol[0] - lat['chi_max_exponent']) < 1e-6
    x = np.array([(lat['N_t'] / s['N_s']) ** 3 for s in sizes])
    bc = np.array([s['beta_c'] for s in sizes])
    be = np.array([max(s['beta_c_err'], 1e-5) for s in sizes])
    A = np.vstack([np.ones_like(x), -x]).T
    cov = np.linalg.inv(A.T @ (A / be[:, None] ** 2))
    sol = cov @ (A.T @ (bc / be ** 2))
    c['beta_c_extrapolation_refit'] = abs(sol[0] - lat['beta_c']) < 1e-6
    c['beta_c_within_one_per_mille_of_published_value'] = abs(lat['beta_c'] - BETA_C_LITERATURE) < 1e-3 * BETA_C_LITERATURE
    c['two_peaks_at_largest_volume'] = bool(sizes[-1]['two_peaks'])
    ok = []
    for s in sizes:
        hi = [p['absL'] for p in s['points'] if p['beta'] >= 5.75]
        lo = [p['absL'] for p in s['points'] if p['beta'] <= 5.60]
        if hi and lo:
            ok.append(min(hi) > 3 * max(lo))
    c['loop_small_below_large_above'] = bool(ok) and all(ok)
    sc = rec['scatter_sectors']
    c['three_sectors_above'] = sorted(sc['sectors_seen']) == [0, 1, 2]
    reg = rec['register']
    growth = register_growth(reg)
    z3, z2 = growth['Z3'], growth['Z2']
    c['register_Z3_growth_above_Z2_all_sizes'] = z3['exponent'] > z2['exponent'] + 3 * math.hypot(z3['err'], z2['err'])
    dips = {name: [h['dip_over_second_peak'] for h in reg[name]['hist_dip']] for name in ('Z3', 'Z2')}
    c['register_Z3_two_peaks_deepen_with_size'] = all(b < a for a, b in zip(dips['Z3'], dips['Z3'][1:])) and dips['Z3'][-1] < 0.25
    c['register_Z2_single_peak_at_every_size'] = min(dips['Z2']) > 0.9
    return {k: bool(v) for k, v in c.items()}


def certificate(monte_carlo=True):
    rec = recorded()
    checks = {}
    # 1 colour factor of the five-trace plaquette action
    checks['1_action.colour_factor_is_wilson_with_three_fifths'] = centre.colour_factor_identity() < 1e-13
    checks['1_action.weak_factor_is_two_fifths'] = abs(
        centre.five_trace_cost(V2=np.diag([np.exp(0.7j), np.exp(-0.7j)])) - 0.4 * (1 - math.cos(0.7))) < 1e-13
    thin = centre.thin_flux_cost()
    checks['1_action.thin_centre_flux_costs'] = [thin[str(f)] for f in range(6)] == ['0', '17/10', '9/10', '4/5', '9/10', '17/10']
    # 2 centre classes
    ft = centre.flux_table()
    checks['2_centre.classes_2_4_act_on_colour_only'] = all(ft[f]['acts_on_colour'] and not ft[f]['acts_on_weak'] for f in ('2', '4'))
    checks['2_centre.class_3_acts_on_weak_only'] = ft['3']['acts_on_weak'] and not ft['3']['acts_on_colour']
    checks['2_centre.colour_flux_rotates_quark_classes'] = ft['2']['classes_rotated'] == ft['4']['classes_rotated'] == [1, 2, 4, 5]
    checks['2_centre.loop_phase_is_triality'] = all(ft['2']['loop_phase_turns'][str(k)] == str((k % 3) * centre.Fraction(1, 3) % 1) for k in range(6))
    # 3 centre transformation on the lattice
    lc = centre.lattice_checks()
    for k, v in lc.items():
        checks['3_lattice.' + k] = v < 1e-12
    vc = centre.vortex_checks()
    checks['3_lattice.region_flux_equals_surface'] = vc['flux_plaquettes'] == vc['region_surface_links'] and vc['all_plaquettes_central']
    checks['3_lattice.region_flux_on_time_like_plaquettes_of_the_slice'] = vc['spatial_flux_plaquettes'] == 0 and vc['flux_only_on_the_slice']
    checks['3_lattice.flux_surface_is_closed'] = vc['cube_closure_violation'] == 0
    checks['3_lattice.loops_rotate_inside_only'] = vc['loops_inside'][0] < 1e-12 and vc['loops_outside'][0] < 1e-12
    # 4 strong coupling
    sc = strong_coupling()
    checks['4_strong_coupling.haar_moments'] = [sc['haar_moments_re_trace'][str(k)] for k in range(1, 6)] == [0.0, 0.5, 0.25, 0.75, 0.9375]
    checks['4_strong_coupling.series_matches_exact_at_small_beta'] = abs(sc['one_plaquette']['0.5']['one_plaquette'] - sc['one_plaquette']['0.5']['series']) < 1e-6
    # 5 kernel gates
    replay = {}
    if monte_carlo:
        kg = kernel_gates()
        replay['kernel_gates'] = kg
        checks['5_kernel.unitarity'] = kg['unitarity_defect'] < 1e-12
        checks['5_kernel.overrelaxation_keeps_action'] = kg['overrelaxation_plaquette_change'] < 1e-12
        s = kg['strong_coupling_beta_1']
        checks['5_kernel.strong_coupling_beta_1'] = abs(s['measured'] - s['expected']) < 4 * s['err'] + 3e-4
    # 6 recorded tables
    for k, v in table_checks(rec).items():
        checks['6_tables.' + k] = v
    # 7 scale
    st = scale.table()
    checks['7_scale.g_star'] = scale.g_star(18, 50) == scale.G_ABOVE and scale.g_star(5, 14) == scale.G_BELOW
    checks['7_scale.ratio'] = abs(st['primary']['ratio'] - 156.5 / 338.0) < 1e-15
    checks['7_scale.temperature'] = abs(st['primary']['T_MeV'] - 156.5 * 340.0 / 338.0) < 1e-9
    checks['7_scale.time_coefficient'] = abs(scale.cosmic_time_seconds(1.0, 1.0) - 2.4203) < 1e-3
    checks['7_scale.ticks'] = abs(st['ticks_g_above'] * st['tick_seconds'] - st['t_seconds_g_above']) < 1e-18
    checks = {k: bool(v) for k, v in checks.items()}
    files = ['su3.py', 'su3_numba.py', 'centre.py', 'scale.py', 'certificate.py', 'lattice_tables.json']
    return dict(schema='colour_deconfinement.v1',
                scope='Pure SU(3) lattice gauge theory at finite temperature as the colour factor of the link system. '
                      'Exact centre statements, recorded Monte Carlo tables, and the placement of the transition on the '
                      'coupling chain. Quarks are absent from the simulation. The scale uses imported lattice QCD results.',
                checks=checks, all_checks_pass=all(checks.values()),
                centre=dict(flux_table=ft, thin_flux_cost=thin,
                            vortex=dict(flux_plaquettes=vc['flux_plaquettes'], region_surface_links=vc['region_surface_links'],
                                        spatial_flux_plaquettes=vc['spatial_flux_plaquettes'], cube_closure_violation=vc['cube_closure_violation'],
                                        loops_inside_turned=vc['loops_inside'][0] < 1e-12, loops_outside_unchanged=vc['loops_outside'][0] < 1e-12)),
                strong_coupling=stable(sc), scale=stable(st), lattice=rec['su3'], gates=rec['gates'],
                plaquette_symmetric=rec['plaquette_symmetric'], scatter_sectors=rec['scatter_sectors'],
                string_tension=rec.get('string_tension'), coupling_scan=rec.get('coupling_scan'), N_t_6=rec.get('N_t_6'),
                register=rec['register'], register_growth=stable(register_growth(rec['register']), 6),
                statements=rec['statements'], references=rec['references'],
                input_sha256={f: hashlib.sha256((HERE / f).read_bytes()).hexdigest() for f in files},
                replay=replay)


def result():
    return json.loads((HERE / 'results.json').read_text())
