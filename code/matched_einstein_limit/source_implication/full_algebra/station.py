"""Full parity-block entropy maximum on the native directed-range instrument."""

from itertools import permutations
from pathlib import Path
import json

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import brentq, linprog, root

from .. import full_station, intrinsic_area
from ..geometry_gain import distance_check, simplex_area_volume
from ..isometric_collar import measure as iso
from model.experiments.interval_echo_einstein_station.collar_pairing.phase_geometry import measure as phase

ORDER = tuple(zip(*permutations(range(4), 2)))


def radar_menus():
    """Both parity sectors include all final misses at the cap and a fair choice of charge."""
    board = iso.pairing.formation.federation(1)
    route = iso.pairing.formation.region(board)
    records = []
    for odd in (False, True):
        by_charge = []
        for charge in (-1, 1):
            h = iso.ranging.generator(board, route, charge, 8, odd)
            kernels = iso.ranging.kernels(h, route, .5, 12)
            by_charge.append(iso.effects(kernels, .5))
        records.append(np.mean(by_charge, axis=0))
    return {cap: np.stack([panel[i] for panel in records], axis=1) + (cap/2)*np.eye(4)
            for i, cap in enumerate((2, 4, 6))}


def directions():
    paulis = iso.basis(False)
    d = np.zeros((31, 2, 4, 4), complex)
    d[:15, 0], d[15:30, 1] = paulis, paulis
    d[30] = np.stack((np.eye(4), -np.eye(4))) / np.sqrt(8)
    return d


def square(r):
    m = np.zeros((4, 4), dtype=r.dtype)
    m[ORDER] = r
    return m


def geometry(r):
    return np.array(simplex_area_volume(square(r)))


def range_gradient(r):
    return np.array([geometry(r + 1e-16j*e).imag / 1e-16 for e in np.eye(12)]).T


def pair(a, b):
    return np.einsum('...bij,bji->...', a, b).real


def log_blocks(x):
    val, u = np.linalg.eigh(x)
    if val.min() <= 0:
        raise ValueError("interior full-algebra calculation requires both blocks faithful")
    return (u * np.log(val)[:, None, :]) @ u.conj().transpose(0, 2, 1)


def entropy(x):
    return -float(pair(log_blocks(x), x))


def gibbs(k):
    val, u = np.linalg.eigh(k)
    w = np.exp(val-val.max())
    w /= w.sum()
    return (u*w[:, None, :]) @ u.conj().transpose(0, 2, 1)


class Survey:
    def __init__(self, effects):
        self.effects = effects
        self.base = effects.trace(axis1=2,axis2=3).sum(axis=1).real/8
        self.centered = effects - self.base[:,None,None,None]*np.eye(4)
        self.directions = directions()
        self.range_map = np.einsum('pbjk,dbkj->pd', self.centered, self.directions).real

    def state(self, ranges, eta, price):
        area, volume = range_gradient(ranges)
        return gibbs(np.einsum('p,pbij->bij', eta*area-price*volume, self.centered))

    def solve(self, eta, target, steps=5):
        z = np.r_[self.base, 0.]
        tracial_volume = geometry(self.base)[1]
        for fraction in np.linspace(0, 1, steps+1)[1:]:
            wanted = tracial_volume + fraction*(target-tracial_volume)
            def equations(args):
                ranges, price = args[:-1], args[-1]
                x = self.state(ranges, fraction*eta, price)
                return np.r_[pair(self.effects, x)-ranges, geometry(ranges)[1]-wanted]
            fit = root(equations, z, tol=2e-11)
            if np.linalg.norm(equations(fit.x)) > 1e-10:
                raise ArithmeticError("range and volume equations failed their residual check")
            z = fit.x
        return self.state(z[:-1], eta, z[-1]), float(z[-1])

    def data(self, x, eta):
        r = pair(self.effects, x)
        area, vol = range_gradient(r) @ self.range_map
        ds = -pair(self.directions, log_blocks(x))
        price = float((ds+eta*area) @ vol/(vol @ vol))
        residual = ds + eta*area - price*vol
        return r, area, vol, ds, price, residual

    def hessian(self, x, eta, step):
        r, _, dv, _, price, _ = self.data(x, eta)
        lam, u = np.linalg.eigh(x)
        h = np.zeros((31,31))
        for b in (0,1):
            w = lam[b]
            diff = w[:,None] - w
            logdiff = np.log(w[:,None]) - np.log(w)
            frechet = np.divide(logdiff, diff, out=np.broadcast_to(1/w, diff.shape).copy(), where=diff != 0)
            local = u[b].conj().T @ self.directions[:,b] @ u[b]
            h -= np.einsum('mij,nij,ij->mn', local.conj(), local, frechet).real
        def geometric(ranges):
            a, v = range_gradient(ranges)
            return eta*a - price*v
        second = np.array([(geometric(r+step*e)-geometric(r-step*e))/(2*step) for e in np.eye(12)])
        h += self.range_map.T @ ((second + second.T)/2) @ self.range_map
        tangent = null_space(dv[None])
        restricted = np.linalg.eigvalsh(tangent.T @ h @ tangent)
        return np.linalg.eigvalsh(h), restricted


def preparation(x):
    """Atom = local Pauli product kept even, or the same atom dephased into two equal blocks."""
    projectors = [(np.eye(2)+sign*a)/2 for a in phase.AXES[1:] for sign in (-1,1)]
    states = np.array([np.kron(a,b) for a in projectors for b in projectors])
    obs = np.concatenate([np.eye(4)[None]/2, iso.basis(False)])
    table = np.einsum('kij,pji->kp', obs, states).real
    def solve_weights(block):
        moments = np.einsum('kij,ji->k', obs, block).real
        result = linprog(np.r_[np.zeros(36),-1], A_eq=np.c_[table, np.zeros(16)], b_eq=moments,
                         A_ub=np.c_[-np.eye(36), np.ones(36)], b_ub=np.zeros(36),
                         bounds=(None,None), method='highs')
        if not result.success:
            raise ArithmeticError(result.message)
        weights = result.x[:36]
        restored = np.einsum('i,ijk->jk', weights, states)
        return weights, restored
    kept, even = solve_weights(x[0]-x[1])
    dephased, mixed = solve_weights(2*x[1])
    reconstructed = np.stack((even + mixed/2, mixed/2))
    return dict(min_kept=float(kept.min()), min_dephased=float(dephased.min()),
                max_error=float(np.max(abs(x-reconstructed))),
                weight_sum=float((kept+dephased).sum()), kept=kept.tolist(), dephased=dephased.tolist())


def held(survey, x, eta, formed):
    """A previously published range-invisible even transfer and a separate regular volume-matched chord."""
    ops = iso.basis(False)
    rows = iso.components(survey.effects[:,0],ops)
    hidden = null_space(rows)
    direction = hidden @ (hidden.T @ iso.components(-iso.pairing.log_density(formed),ops))
    even = np.einsum('k,kij->ij', direction/np.linalg.norm(direction),ops)
    blind = np.stack([even,np.zeros_like(even)])
    read = survey.data(x, eta)
    _, da, dv, ds, _, _ = read
    shape = da - dv*(dv@da)/(dv@dv)
    shape /= np.linalg.norm(shape)
    physical = np.einsum('k,kbij->bij',shape,survey.directions)
    adjust = np.einsum('k,kbij->bij',dv/(dv@dv),survey.directions)
    target = geometry(pair(survey.effects,x))[1]
    tests = []
    for eps in (1e-3, 5e-4):
        records = []
        for sign in (-1,1):
            y = x+sign*eps*physical
            shift = brentq(lambda h: geometry(pair(survey.effects,y+h*adjust))[1]-target,
                           -.001,.001,xtol=1e-15)
            y = y+shift*adjust
            area, vol = geometry(pair(survey.effects,y))
            records.append((area,vol,entropy(y)))
        low, high = np.array(records)
        A,V,S = (high-low)/(2*eps)
        tests.append(dict(step=eps, entropy=S, area=A, residual=S+eta*A,
                          max_volume_error=max(abs(v[1]-target) for v in records)))
    return dict(hidden_range_error=float(np.max(abs(pair(survey.effects,blind)))),
                hidden_entropy=float(-pair(blind,log_blocks(x))),
                shape_entropy=float(ds@shape),shape_area=float(da@shape),finite=tests)


def run():
    menus = radar_menus()
    survey = Survey(menus[6])
    product, difference = intrinsic_area.head_source()
    prior = full_station.Station(intrinsic_area.survey_corners(intrinsic_area.radar_banks(),1,6),
                                 product,difference)
    eta = prior.read(np.array([.25,.5,.125+1/256]))['eta']
    target = prior.geometry(np.array([.25,.5,.125]))[0][1]
    x, price = survey.solve(eta, target)
    results = {}
    for cap in (4,6):
        instrument = Survey(menus[cap])
        r, _, _, _, own_price, residual = instrument.data(x,eta)
        geom = distance_check(square(r))
        results[cap] = dict(area_volume=geometry(r).tolist(),price=own_price,
                            max_residual=float(np.linalg.norm(residual)),
                            cycles=geom['cycle_ratios'],necessary_log_error_floor=geom['necessary_log_error_floor'],
                            directional_log_residual=geom['log_residual_max'])
    h, constrained = survey.hessian(x,eta,1e-4)
    hc = survey.hessian(x,eta,1e-5)[1]
    p = preparation(x)
    return dict(ordinary_eta=float(eta),target_volume=float(target),price=price,
                odd=float(np.trace(x[1]).real), spectrum=np.linalg.eigvalsh(x).tolist(),
                x_real=x.real.tolist(),x_imag=x.imag.tolist(), reports=results, preparation=p,
                full_hessian_extremes=h[[0,-1]].tolist(),constrained_eigenvalues=constrained.tolist(),
                hessian_check=float(max(abs(hc-constrained))),held=held(survey,x,eta,product+difference))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    report = json.dumps(run(),indent=2,allow_nan=False) + '\n'
    if parser.parse_args().write:
        Path(__file__).with_name("results.json").write_text(report)
    else:
        print(report,end='')
