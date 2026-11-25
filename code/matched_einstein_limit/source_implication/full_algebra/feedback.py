"""Select the six-unit ensemble by feedback from tracial records and ordinary calibration.

For B(X)=Gibbs(gradient T(X)), put U=B(X), W=B(Y). Symmetrized relative
entropy is Tr[(U-W)(gradient T(X)-gradient T(Y))] >= ||U-W||_1^2 by
Pinsker with natural logs. The certified bilinear Hessian gives an upper
bound .72 ||U-W||_1 ||X-Y||_1 on the same expression. Banach iteration
therefore selects the unique fixed point at any certified price, from any
state. It maximizes S+T globally by the relative-entropy gap in PROOF.md.
Adding the optimality inequalities for two prices makes V(X_price)
nonincreasing. Continuity follows from uniform contraction. If two prices
attain the target they share the same maximizer, by strict optimality.
Subtracting their Gibbs equations then forces dV=0 there. The certificate
below checks dV != 0 on the entire enclosed target-state ball.

The twelve effect tables are calibrated on known probes. Affine range
expectations can then be evaluated on trial matrices without deploying
them. The final positive lottery loads the selected ensemble. Direct mean
readings instead require separate copies and all timeouts counted. For m
fixed-price updates let w_i=(upper_i-lower_i)/2 from the certified spectral
box, and clip each empirical range to its own box. With independent reads
in [0,3], N_i >= 9 log(24m/alpha)/(2 delta^2 w_i^2) per mean gives all
errors <= delta*w_i with probability >= 1-alpha by Hoeffding and a union
bound. On the segment of clipped ranges the certified bound sums
|T_ij| w_i w_j < .72. Subtract the spectral midpoint from each E_i. Then
the gradient error modulo scalars is <= .72 delta in operator norm, and
symmetrized Pinsker gives the same trace-norm Gibbs-update error. Including
preparation error, a per-update trace-norm error e gives a final error
<= .72^m times the initial error plus e*(1-.72^m)/(1-.72). This estimate
assumes the calibrated effects, gain and price. The exact sign checks below
use the independently enclosed calibration. The construction supplies an
ensemble controller. Its update differs from the free propagator.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import certificate as cert
from . import global_survey, station

PRICE = (50.98, 51.20)
Q = .72


def calibration():
    """Compare the specified ordinary interventions before selecting any source."""
    product, difference = station.intrinsic_area.head_source()
    corners = station.intrinsic_area.survey_corners(station.intrinsic_area.radar_banks(), 1, 6)
    ordinary = station.full_station.Station(corners, product, difference)
    ref = np.array([.25, .5, .125])
    target = ordinary.geometry(ref)[0][1]
    eta = ordinary.read(ref + [0, 0, 1/256])["eta"]
    return float(eta), float(target)


def fixed_price(survey, eta, price, atom_margins=None):
    if not PRICE[0] <= price <= PRICE[1]:
        raise ValueError("price outside the certified contraction domain")
    x = np.stack((np.eye(4), np.eye(4)))/8
    for n in range(1, 160):
        y = survey.state(station.pair(survey.effects, x), eta, price)
        change = sum(abs(np.linalg.eigvalsh(y-x)).flat)
        x = y
        if atom_margins is not None:
            prep = station.preparation(x)
            atom_margins.append(min(prep["min_kept"], prep["min_dephased"]))
        if change < 3e-13:
            return x, n, float(change)
    raise ArithmeticError("tracial Gibbs iteration exceeded the prescribed stopping residual")


def select(survey, eta, target):
    def trial(price):
        x, n, residual = fixed_price(survey, eta, price)
        difference = station.geometry(station.pair(survey.effects, x))[1] - target
        return x, n, residual, float(difference)

    lo, hi = PRICE
    endpoints = [trial(lo)[-1], trial(hi)[-1]]
    if not endpoints[0] > 0 > endpoints[1]:
        raise ValueError("complete range readings do not bracket the calibration")
    outer = 0
    while hi-lo > 2e-10:
        middle = (lo+hi)/2
        if trial(middle)[-1] > 0:
            lo = middle
        else:
            hi = middle
        outer += 1
    center = (lo+hi)/2
    x, n, residual, difference = trial(center)
    return x, center, dict(price=center, outer_steps=outer, final_tracial_iterations=n,
                           last_update_trace_norm=residual, volume_difference=difference,
                           original_endpoint_volume_differences=endpoints)


class ExactSurvey:
    """Check the fixed-point residual with rational effect centers and interval arithmetic."""

    def __init__(self):
        cert.iv.dps = 28
        nums, unit, error = cert.exact_effects()
        self.error = sp.Rational(error)
        self.effects = cert.rational_effects(nums, unit)
        self.eta, self.target = cert.calibration(self.effects, self.error)
        self.ops = [cert.prior.sym_matrix(h) for h in station.iso.basis(False)]
        self.sqrt8 = cert.iv.sqrt(8)
        self.width = np.array([cert.ri(sp.Rational(hi-lo, 20000))
                               for lo, hi in global_survey.BOUNDS], dtype=object)
        jac = []
        for pair in self.effects:
            local = [cert.ri_band(cert.inner(h,p), 2*self.error) for h in pair for p in self.ops]
            parity = cert.ri_band((pair[0]-pair[1]).trace(), 8*self.error)/self.sqrt8
            jac.append(local+[parity])
        self.jac = np.array(jac, dtype=object)
        assert self.eta.a > cert.ri("23.1042361") and self.eta.b < cert.ri("23.1042363")

    def check(self, state, prices):
        # Evaluate the log-gradient score on exactly trace-normalized Hermitian rational blocks.
        x = [cert.prior.sym_matrix(b) for b in state]
        x = [(b+b.H)/2 for b in x]
        tr = sum(b.trace() for b in x)
        x = [b/tr for b in x]
        ent = []
        for b in x:
            vals, err = cert.prior.entropy_coordinates(b, 0, sp.Rational(".020"), self.ops+[sp.eye(4)])
            ent.append([cert.ri_band(v, (2 if i < 15 else 4)*err) for i, v in enumerate(vals)])
        entropy = ent[0][:15]+ent[1][:15]+[(ent[0][-1]-ent[1][-1])/self.sqrt8]
        r = [cert.ri_band(sum(cert.inner(h,b) for h,b in zip(p,x)), self.error) for p in self.effects]
        a, v = cert.enc.geometry_jet(r)
        price = cert.iv.mpf(list(prices))
        residual = np.array(entropy, dtype=object) + (a.g@self.jac)*self.eta - (v.g@self.jac)*price
        score = cert.iv.sqrt(sum(max(abs(f.a),abs(f.b))**2 for f in residual))
        # Pinsker between B(X) and X bounds their trace distance by the norm of
        # the traceless score. The Hilbert-Schmidt score bounds its operator norm.
        radius = score/(1-cert.ri(".72"))
        ranges = [u+cert.iv.mpf([-1,1])*w*radius for u,w in zip(r,self.width)]
        _, at_fixed_point = cert.enc.geometry_jet(ranges)
        dv = at_fixed_point.g@self.jac
        vnorm = cert.iv.sqrt(sum(max(abs(f.a),abs(f.b))**2 for f in dv))
        change = cert.iv.mpf([-1,1])*radius*vnorm
        return v.v+change-self.target, radius.b, dv[-1], x


def run():
    eta, v0 = calibration()
    survey = station.Survey(station.radar_menus()[6])
    x, price, replay = select(survey, eta, v0)
    atom_margins = []
    fixed_price(survey, eta, price, atom_margins)
    replay["smallest_nominal_atom_weight_among_final_price_updates"] = min(atom_margins)
    exact = ExactSurvey()
    # Round endpoints out around the controller's independently located price.
    # The sign proof is independent of the sign decisions in floating-point bisection.
    lower = (int(price*1e5)-1)/1e5
    upper = (int(price*1e5)+2)/1e5
    signs = []
    for p in (lower, upper):
        y = fixed_price(survey, eta, p)[0]
        signature = exact.check(y, (str(p), str(p)))[0]
        signs.append(signature)
    assert signs[0].a > 0 > signs[1].b
    _, radius, volume_direction, rational_x = exact.check(x, (str(lower), str(upper)))
    assert volume_direction.a > .1
    weights = station.preparation(x)
    # Correct both target-block mixtures by the product Pauli frame, sigma_min=1.
    # Each target changes by at most twice the full trace-norm bound.
    rays = [(sp.eye(2)+s*p)/2 for p in (sp.Matrix([[0,1],[1,0]]),
            sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1)) for s in (-1,1)]
    atoms = [sp.kronecker_product(a,b) for a in rays for b in rays]
    margins = {}
    for k, goal in (("kept",rational_x[0]-rational_x[1]), ("dephased",2*rational_x[1])):
        ws = [sp.Rational(w) for w in weights[k]]
        err = goal-sum((w*b for w,b in zip(ws,atoms)),sp.zeros(4))
        margin = cert.ri(min(ws))-cert.iv.sqrt(cert.ri(cert.inner(err,err)))-2*radius
        assert margin.a > .003
        margins[k] = float(np.nextafter(float(margin.a), -np.inf))
    _, _, _, _, _, residual = survey.data(x, eta)
    held, held_price = survey.solve(eta, v0)
    agreement = dict(trace_norm=float(abs(np.linalg.eigvalsh(x-held)).sum()),
                     absolute_price_difference=float(abs(price-held_price)))
    return dict(calibration=dict(ordinary_gain=eta, comparison_volume=v0),
                feedback=replay, held_continuation_agreement=agreement, native_odd_weight=float(x[1].trace().real),
                physical_price_bracket=[lower, upper],
                true_state_error_center="Hermitian trace-normalized rational feedback iterate",
                certified_endpoint_volume_differences=[str(v) for v in signs],
                true_state_trace_distance_upper=float(np.nextafter(float(radius), np.inf)),
                volume_parity_gradient_lower=float(np.nextafter(float(volume_direction.a), -np.inf)),
                nominal_full_gradient_residual=float(np.linalg.norm(residual)),
                native_atoms=dict(count=72, nominal_weight_sum=weights["weight_sum"],
                    nominal_max_reconstruction_error=weights["max_error"],
                    corrected_minimum_weights_lower=margins,
                    kept=weights["kept"], dephased=weights["dephased"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = json.dumps(run(), indent=2, allow_nan=False)+"\n"
    if parser.parse_args().write:
        Path(__file__).with_name("feedback.json").write_text(result)
    else:
        print(result, end="")
