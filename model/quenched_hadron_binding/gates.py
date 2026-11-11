"""Fast gates of the lattice code (plain numpy, small lattices, fixed seeds)."""
import numpy as np

from . import lattice as lt


def one_plaquette_average(beta, n=240):
    """Exact <Re tr U / 3> with weight exp((beta / 3) Re tr U) by Weyl integration over the maximal torus."""
    t = 2 * np.pi * (np.arange(n) + 0.5) / n
    a, b = np.meshgrid(t, t, indexing='ij')
    c = -a - b
    e = [np.exp(1j * a), np.exp(1j * b), np.exp(1j * c)]
    haar = np.abs(e[0] - e[1]) ** 2 * np.abs(e[0] - e[2]) ** 2 * np.abs(e[1] - e[2]) ** 2
    retr = np.cos(a) + np.cos(b) + np.cos(c)
    w = haar * np.exp(beta / 3.0 * retr)
    return float((w * retr / 3.0).sum() / w.sum())


def gates():
    """Returns (checks, numbers): named booleans and the measured values behind them."""
    checks, num = {}, {}
    rng = np.random.default_rng(5701)
    g, g5 = lt.gammas()
    c = lt.charge_conjugation()

    # 1. gamma algebra
    cl = max(np.abs(g[m] @ g[n] + g[n] @ g[m] - 2 * (m == n) * np.eye(4)).max() for m in range(4) for n in range(4))
    he = max(np.abs(x - x.conj().T).max() for x in g + [g5])
    a5 = max(np.abs(g5 @ x + x @ g5).max() for x in g)
    cc = max(np.abs(c @ x @ np.linalg.inv(c) + x.T).max() for x in g)
    checks['1_gamma.clifford_algebra'] = bool(cl < 1e-14 and he < 1e-14)
    checks['1_gamma.gamma5_anticommutes'] = bool(a5 < 1e-14 and np.abs(g5 @ g5 - np.eye(4)).max() < 1e-14)
    checks['1_gamma.charge_conjugation'] = bool(cc < 1e-14)

    # 2. gauge sector: strong coupling, invariance, overrelaxation
    beta = 1.0
    u = lt.hot_links((4, 4, 4, 4), rng)
    vals = []
    for n in range(140):
        lt.heatbath_sweep(u, beta, rng, overrelax=0)
        if n >= 40:
            vals.append(lt.plaquette(u))
    vals = np.array(vals)
    blocks = vals.reshape(10, -1).mean(1)
    mean, err = float(vals.mean()), float(blocks.std(ddof=1) / np.sqrt(10))
    u1 = one_plaquette_average(beta)
    ref = u1 + 4 * u1 ** 5
    num['2_gauge.strong_coupling'] = dict(beta=beta, lattice=[4, 4, 4, 4], plaquette=mean, err=err, exact_one_plaquette=u1,
                                          with_cube_term=ref)
    checks['2_gauge.strong_coupling_plaquette'] = bool(abs(mean - ref) < 5e-4 + 4 * err)
    dims = (4, 4, 4, 8)
    u = lt.hot_links(dims, rng)
    for _ in range(12):
        lt.heatbath_sweep(u, 5.7, rng, overrelax=1)
    p0 = lt.plaquette(u)
    uo = u.copy()
    lt.heatbath_sweep(uo, 5.7, rng)                       # one more heat-bath sweep changes the action
    num['2_gauge.rough_lattice'] = dict(lattice=list(dims), beta=5.7, sweeps=12, plaquette=float(p0))
    unit = max(np.abs(u @ lt.dag(u) - np.eye(3)).max(), np.abs(np.linalg.det(u) - 1).max())
    checks['2_gauge.links_stay_in_su3'] = bool(unit < 1e-12)
    gt = lt.random_su3(dims, rng)
    ug = lt.gauge_transform(u, gt)
    v = lt.ape_smear(u, 0.5, 6)
    vg = lt.ape_smear(ug, 0.5, 6)
    w = [lt.wilson_loop(u, v, r, t) for r in (1, 2) for t in (1, 2)]
    wg = [lt.wilson_loop(ug, vg, r, t) for r in (1, 2) for t in (1, 2)]
    checks['2_gauge.plaquette_gauge_invariant'] = bool(abs(lt.plaquette(ug) - p0) < 1e-12)
    checks['2_gauge.smeared_loops_gauge_invariant'] = bool(np.abs(np.array(w) - np.array(wg)).max() < 1e-12)
    checks['2_gauge.ape_links_in_su3_and_covariant'] = bool(
        np.abs(v @ lt.dag(v) - np.eye(3)).max() < 1e-12 and np.abs(vg - lt.gauge_transform(v, gt)).max() < 1e-11)
    # the 1 x 1 loop with thin links is the time-like plaquette
    pt = 0.0
    for i in range(3):
        p = u[i] @ np.roll(u[3], -1, axis=i) @ lt.dag(np.roll(u[i], -1, axis=3)) @ lt.dag(u[3])
        pt += np.trace(p, axis1=-2, axis2=-1).real.mean() / 9.0
    checks['2_gauge.loop_1x1_is_plaquette'] = bool(abs(lt.wilson_loop(u, u, 1, 1) - pt) < 1e-13)
    num['2_gauge.loops'] = dict(W=[float(x) for x in w])

    # 3. Wilson-Dirac operator
    uf = lt.fermion_links(u)
    kap = 0.125
    phi = rng.normal(size=dims + (4, 3)) + 1j * rng.normal(size=dims + (4, 3))
    chi = rng.normal(size=dims + (4, 3)) + 1j * rng.normal(size=dims + (4, 3))
    lhs, rhs = np.vdot(phi, lt.apply_mdag(uf, kap, chi)), np.vdot(lt.apply_m(uf, kap, phi), chi)
    num['3_dirac.gamma5_hermiticity'] = float(abs(lhs - rhs) / abs(rhs))
    checks['3_dirac.gamma5_hermiticity'] = bool(abs(lhs - rhs) / abs(rhs) < 1e-12)
    s, res = lt.propagator(uf, kap)
    sg, resg = lt.propagator(lt.fermion_links(ug), kap)
    g0 = gt[0, 0, 0, 0]
    expect = np.einsum('...ac,...scrd,bd->...sarb', gt, s, np.conj(g0))
    num['3_dirac.true_residual'] = float(max(res, resg))
    checks['3_dirac.solver_true_residual'] = bool(max(res, resg) < 1e-9)
    dcov = float(np.abs(sg - expect).max() / np.abs(s).max())
    num['3_dirac.gauge_covariance'] = dcov
    checks['3_dirac.propagator_gauge_covariant'] = bool(dcov < 1e-7)
    # free field against the momentum-space form
    fdims = (4, 4, 4, 16)
    cold = lt.cold_links(fdims)
    sf, resf = lt.propagator(lt.fermion_links(cold), 0.120, tol=1e-12)
    sx = lt.free_propagator(fdims, 0.120)
    ref = np.zeros_like(sf)
    for a in range(3):
        ref[..., :, a, :, a] = sx
    dfree = float(np.abs(sf - ref).max() / np.abs(ref).max())
    pion, rho = lt.meson_correlators(sf)
    pion_ref = 3.0 * (np.abs(sx) ** 2).sum(axis=(0, 1, 2, 4, 5))
    rho_ref = 0.0
    for i in range(3):
        a = g5 @ g[i]
        rho_ref = rho_ref + np.einsum('ij,xyztjk,kl,xyztil->t', a, sx, a.conj().T, np.conj(sx)).real
    num['3_dirac.free_field'] = dict(lattice=list(fdims), kappa=0.120, propagator=dfree,
                                     pion=float(np.abs(pion / pion_ref - 1).max()), rho=float(np.abs(rho / rho_ref - 1).max()))
    checks['3_dirac.free_propagator'] = bool(dfree < 1e-8)
    checks['3_dirac.free_pion_correlator'] = bool(np.abs(pion / pion_ref - 1).max() < 1e-7)
    checks['3_dirac.free_rho_correlator'] = bool(np.abs(rho / rho_ref - 1).max() < 1e-7)

    # 4. hadron correlators
    pi0, rho0 = lt.meson_correlators(s)
    pig, rhog = lt.meson_correlators(sg)
    cp, cm = lt.nucleon_correlators(s)
    cpg, cmg = lt.nucleon_correlators(sg)
    dinv = float(max(np.abs(pi0 / pig - 1).max(), np.abs(rho0 / rhog - 1).max(), np.abs(cp - cpg).max() / np.abs(cp).max(),
                     np.abs(cm - cmg).max() / np.abs(cm).max()))
    num['4_hadrons.gauge_invariance_point'] = dinv
    checks['4_hadrons.correlators_gauge_invariant'] = bool(dinv < 1e-7)
    rs = rng.normal(size=(1, 1, 1, 2, 4, 3, 4, 3)) + 1j * rng.normal(size=(1, 1, 1, 2, 4, 3, 4, 3))
    ap, am = lt.nucleon_correlators(rs)
    dw = 0.0
    for t in range(2):
        ref_site = lt.nucleon_site_reference(rs[0, 0, 0, t])
        dw = max(dw, abs(ap[t] - (ref_site[0, 0] + ref_site[1, 1])) / abs(ap[t]),
                 abs(am[t] - (ref_site[2, 2] + ref_site[3, 3])) / abs(am[t]))
    num['4_hadrons.nucleon_wick_sum'] = float(dw)
    checks['4_hadrons.nucleon_contraction_equals_wick_sum'] = bool(dw < 1e-12)
    checks['4_hadrons.meson_correlators_positive'] = bool((pi0 > 0).all() and (rho0 > 0).all())
    ss, _ = lt.propagator(uf, kap, v=v, smear=(2.0, 8))
    ssg, _ = lt.propagator(lt.fermion_links(ug), kap, v=vg, smear=(2.0, 8))
    a1, a2 = lt.nucleon_correlators(ss), lt.nucleon_correlators(ssg)
    m1, m2 = lt.meson_correlators(ss), lt.meson_correlators(ssg)
    dsm = float(max(np.abs(m1[0] / m2[0] - 1).max(), np.abs(m1[1] / m2[1] - 1).max(),
                    np.abs(a1[0] - a2[0]).max() / np.abs(a1[0]).max()))
    num['4_hadrons.gauge_invariance_smeared'] = dsm
    checks['4_hadrons.smeared_correlators_gauge_invariant'] = bool(dsm < 1e-7)
    # free nucleon: three quarks at rest, mass 3 ln(1 + m_q); analytic propagator on a long lattice
    ndims = (6, 6, 6, 32)
    kn = 0.110
    sxn = lt.free_propagator(ndims, kn)
    sn = np.zeros(ndims + (4, 3, 4, 3), complex)
    for a in range(3):
        sn[..., :, a, :, a] = sxn
    cpn, cmn = lt.nucleon_correlators(sn)
    big_t = ndims[3]
    meff = np.log(np.abs(cpn[1:big_t // 2].real / cpn[2:big_t // 2 + 1].real))
    mq = np.log(1.0 + 1.0 / (2 * kn) - 4.0)
    num['4_hadrons.free_nucleon'] = dict(lattice=list(ndims), kappa=kn, three_quark_mass=float(3 * mq),
                                         effective_mass_t10=float(meff[9]), effective_mass_t12=float(meff[11]))
    checks['4_hadrons.free_nucleon_is_three_quarks'] = bool(abs(meff[9] / (3 * mq) - 1) < 0.03 and abs(meff[11] / (3 * mq) - 1) < 0.03)
    refl = float(np.abs(cpn[1:big_t].real + cmn[big_t - 1:0:-1].real).max() / abs(cpn[1].real))
    checks['4_hadrons.nucleon_time_reflection'] = bool(refl < 1e-10 and cpn[3].real > 0)
    pin, _ = lt.meson_correlators(sn)
    mpi = float(np.log(pin[10] / pin[11]))
    num['4_hadrons.free_pion'] = dict(effective_mass_t10=mpi, two_quark_mass=float(2 * mq))
    checks['4_hadrons.free_pion_is_two_quarks'] = bool(abs(mpi / (2 * mq) - 1) < 0.02)
    return checks, num
