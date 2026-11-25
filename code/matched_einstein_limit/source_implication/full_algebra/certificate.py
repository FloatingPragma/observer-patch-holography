"""Certify an interior exact fixed-volume maximum near the published numerical state."""

from pathlib import Path
from itertools import combinations
import json

import numpy as np
import sympy as sp

from . import enclosure as enc
from . import global_survey
from .enclosure import iv
from .native_precision import exact_effects, prior
from .station import ORDER, iso


def ri(x):
    x = sp.Rational(sp.expand(x) if isinstance(x,sp.Basic) else x)
    return iv.mpf(int(x.p))/int(x.q)


def ri_band(center, error):
    return ri(center) + iv.mpf([-1,1])*ri(error)


def inner(a, b):
    return sp.expand(sp.re(sp.trace(a*b)))


def rational_effects(nums, unit):
    e = []
    for pair in nums:
        blocks = []
        for re, im in pair:
            h = (sp.Matrix(re)+sp.I*sp.Matrix(im))/unit
            blocks.append(3*sp.eye(4)+(h+h.H)/2)
        e.append(blocks)
    return e


def calibration(e, eps):
    from .native_precision import formed_state

    rho = formed_state()
    left = sp.Matrix(2,2,lambda i,j: sum(rho[2*i+b,2*j+b] for b in range(2)))
    right = sp.Matrix(2,2,lambda i,j: sum(rho[2*b+i,2*b+j] for b in range(2)))
    fresh = sp.kronecker_product(left,right)
    delta = rho-fresh
    zero = sp.zeros(4)
    dirs = ((-fresh,fresh),(delta,zero),(-delta,delta))
    u,s = sp.Rational(1,4),sp.Rational(1,2)
    rho_err = sp.Rational(1,10**18)
    def norm_trace_bound(block):
        return 2*iv.sqrt(ri(inner(block,block)))
    def blocks(v):
        return [(1-u)*fresh+(s-v)*delta, u*fresh+v*delta]
    def ranges(pair):
        # Each full range effect lies in [0, 3I], and the state trace-norm error <4*rho_err.
        return [ri_band(sum(inner(h,b) for h,b in zip(p,pair)),eps+12*rho_err) for p in e]
    baseline = blocks(u*s)
    _, vol = enc.geometry_jet(ranges(baseline))
    ref = blocks(u*s + sp.Rational(1,256))
    a,v = enc.geometry_jet(ranges(ref))
    jac = np.array([[ri(sum(inner(h,b) for h,b in zip(p,d)))+iv.mpf([-1,1])*(
        ri(eps)*sum(norm_trace_bound(b) for b in d)+ri(18*rho_err)) for d in dirs] for p in e],dtype=object)
    da,dv = a.g@jac, v.g@jac
    adjust = -dv[0]/dv[2]
    entropy = np.array([iv.mpf(0)]*3,dtype=object)
    for i in range(2):
        exact_dirs = [row[i] for row in dirs]
        vals,err = prior.entropy_coordinates(ref[i],4*rho_err,sp.Rational('.002'),exact_dirs)
        for k,b in enumerate(exact_dirs):
            # The approximate direction differs by at most 3*rho_err in each block.
            # The true block log has operator norm below 7 at the checked minimum >.002.
            entropy[k] += ri(vals[k])+iv.mpf([-1,1])*(ri(err)*norm_trace_bound(b)+ri(21*rho_err))
    gain = -(entropy[0]+adjust*entropy[2])/(da[0]+adjust*da[2])
    return gain, vol.v


def signed_transfers(x, e, eps, widths, ranges, radius):
    """Check an entropy-preserving, marginal-preserving volume adjuster at every enclosed state."""
    pauli = [sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1)]
    generators = [sp.kronecker_product(sp.eye(2),p) for p in pauli]
    bloch = [ri(sum(inner(b,p) for b in x))+iv.mpf([-1,1])*iv.sqrt(8)*radius for p in generators]
    length = iv.sqrt(sum(b**2 for b in bloch))
    assert length.a > 0
    tangent = []
    for h,w in zip(e,widths):
        coeff = []
        for g in generators:
            dirs = [(-sp.I/2*(g*b-b*g)).applyfunc(sp.expand) for b in x]
            coeff.append(ri(sum(inner(effect,d) for effect,d in zip(h,dirs))) +
                         iv.mpf([-1,1])*(w*radius+ri(eps)))
        tangent.append(sum(b*c for b,c in zip(bloch,coeff))/length)
    _,v = enc.geometry_jet(ranges)
    volume = sum(g*t for g,t in zip(v.g,tangent))
    assert volume.b < ri('-2e-5')
    total = x[0]+x[1]
    right = sp.Matrix(2,2,lambda i,j: sum(total[2*k+i,2*k+j] for k in range(2)))
    left = [sp.Matrix(2,2,lambda i,j: sum(block[2*i+k,2*j+k] for k in range(2))) for block in x]
    transfers = [[total/2 - b for b in x], [sp.kronecker_product(a,right)-b for a,b in zip(left,x)]]
    slopes = []
    for transfer in transfers:
        central = iv.mpf(0)
        for block,t in zip(x,transfer):
            vals,err = prior.entropy_coordinates(block,0,sp.Rational('.02'),[t])
            central += ri_band(vals[0],2*err)
        # T has trace norm <=2. Product-of-marginals minus identity changes by
        # <=3||delta state||_1, and dephasing minus identity by <=2||delta state||_1.
        uncertainty = (2/(ri('.02')-radius) + 3*iv.sqrt(8)*4)*radius
        slopes.append(central + iv.mpf([-1,1])*uncertainty)
    assert slopes[0].a > ri('.409') and slopes[1].a > ri('.231')
    energy = ri(sum(x[0][j,j]-x[1][j,j] for j in range(4))/2)+iv.mpf([-1,1])*iv.sqrt(8)*radius/2
    return dict(right_adjuster_volume=str(volume), ordinary_energy=str(energy),
                ordinary_entropy=str(slopes[0]), pairing_entropy=str(slopes[1]))


def matched_counterpreparation(effects, eps, target, eta):
    """A faithful native product lottery at the same volume with nonstationary entropy plus area."""
    diagonal = np.array([[[ri_band(h[j,j], eps) for j in range(4)] for h in pair]
                         for pair in effects], dtype=object)
    q = ri(sp.Rational(1, 64))
    fraction = ri(sp.Rational(7, 8))
    direction = np.array([ri(0), ri(0), fraction, -fraction], dtype=object)
    at_q = diagonal[:,0]*(1-q) + diagonal[:,1]*q

    def read(z):
        weights = np.array([ri(sp.Rational(1,32))]*4, dtype=object)
        weights[2] += fraction*z
        weights[3] += fraction*(1-z)
        area, volume = enc.geometry_jet(at_q @ weights)
        return weights, area, volume

    lower, upper = ri('.79363406'), ri('.79363408')
    assert read(lower)[2].v.a > target.b and read(upper)[2].v.b < target.a
    weights, area, volume = read(iv.mpf([lower.a, upper.b]))
    z_ranges = at_q @ direction
    q_ranges = (diagonal[:,1]-diagonal[:,0]) @ weights
    vz = volume.g @ z_ranges
    assert vz.b < ri('-.031')
    dz = -(volume.g @ q_ranges) / vz
    entropy_q = iv.ln((1-q)/q) + fraction*iv.ln(weights[3]/weights[2])*dz
    area_q = area.g @ (q_ranges + z_ranges*dz)
    slope = entropy_q + eta*area_q
    assert slope.a > ri('7.47') and slope.b < ri('7.49'), str(slope)
    return dict(odd_weight="1/64", preparation_interval=[".79363406", ".79363408"],
                volume_slope=str(vz), fixed_volume_objective_slope=str(slope))


def verify(native=None):
    iv.dps = 28
    if native is None:
        native = exact_effects()
    nums, unit, epsilon = native
    eps = sp.Rational(epsilon)
    e = rational_effects(nums,unit)
    record = json.loads(Path(__file__).with_name("results.json").read_text())
    x = [prior.sym_matrix(h) for h in np.array(record["x_real"])+1j*np.array(record["x_imag"])]
    x = [(h+h.H)/2 for h in x]
    x = [h / sum(m.trace() for m in x) for h in x]
    radius = ri('0.0001')
    # The semidefinite inequalities below use exact Hermitian leading principal minors.
    for h in x:
        for p in (h - sp.Rational('.02')*sp.eye(4), sp.Rational('.30')*sp.eye(4)-h):
            assert all(sp.re(p[:k,:k].det())>0 for k in range(1,5))
    vectors = [prior.sym_matrix(h) for h in iso.basis(False)]
    ent = []
    for h in x:
        coeff, le = prior.entropy_coordinates(h,0,sp.Rational('.02'),vectors+[sp.eye(4)])
        ent.append([ri_band(coeff[i], (2 if i<15 else 4)*le) for i in range(16)])
    sqrt8 = iv.sqrt(8)
    entropy_score = ent[0][:15]+ent[1][:15]+[(ent[0][-1]-ent[1][-1])/sqrt8]
    jac, r, widths = [], [], []
    for p in e:
        trace = (p[0].trace()+p[1].trace())/8
        width_sq = sum(inner(v-trace*sp.eye(4),v-trace*sp.eye(4)) for v in p)
        widths.append(iv.sqrt(ri(width_sq))+sqrt8*ri(eps))
        r.append(ri_band(sum(inner(h,v) for h,v in zip(x,p)), eps))
        jac.append([ri_band(inner(v,b), 2*eps) for v in p for b in vectors]
                   + [ri_band((p[0]-p[1]).trace(), 8*eps)/sqrt8])
    jac = np.array(jac,dtype=object)
    lam = ri('51.0886428')
    eta = iv.mpf(['23.1042361','23.1042363'])
    target = iv.mpf(['2.763471835','2.763471839'])
    physical_eta, physical_volume = calibration(e,eps)
    assert eta.a < physical_eta.a and physical_eta.b < eta.b, str(physical_eta)
    assert target.a < physical_volume.a and physical_volume.b < target.b, str(physical_volume)
    a,v = enc.geometry_jet(r)
    mismatch = abs(v.v-target)
    assert mismatch.b < ri('1e-8')
    grad = entropy_score + ((a.g*eta-v.g*lam)@jac)
    slope = iv.sqrt(sum(max(abs(g.a),abs(g.b))**2 for g in grad))
    assert slope.b < ri('1e-6')
    region = [a+iv.mpf([-1,1])*w*radius for a,w in zip(r,widths)]
    ar,vr = enc.geometry_jet(region)
    # All exact range directional derivatives obey |dr_i| <= widths[i]*||delta||_HS.
    oriented = (vr.g@jac[:,-1])/1
    assert oriented.a > ri('.10')
    broader = iv.mpf(['50.98','51.20'])
    hg = ar.h*eta-vr.h*broader
    curvature = sum(abs(v)*w for v,w in zip(hg.flat,np.outer(widths,widths).flat))
    assert curvature.b < ri('1.25')
    # Thus for center price, the Hessian lies between -52 and -2 on this ball.
    # 1/(.02-r)+1.25 < 52, 1/(.30+r)-1.25 > 2.
    assert 1/(ri('.02')-radius)+curvature.b < 52
    assert 1/(ri('.30')+radius)-curvature.b > 2
    near = ri('2e-8')
    assert oriented.a*near > mismatch.b
    assert slope.b*radius-radius**2 < -slope.b*near-26*near**2
    state_radius = ri('2e-7')
    assert 2*state_radius > slope.b and slope.b*state_radius-state_radius**2 < -slope.b*near-26*near**2
    true_ranges = [a+iv.mpf([-1,1])*w*state_radius for a,w in zip(r,widths)]
    log = {edge:iv.ln(q) for edge,q in zip(zip(*ORDER),true_ranges)}
    def skew(i,j):
        return log[i,j]-log[j,i]
    floor = max(abs(skew(i,j)+skew(j,k)+skew(k,i)).a for i,j,k in combinations(range(4),3))/6
    assert floor > ri('.00254')
    # At a constrained stationary point, the true multiplier stays in broader.
    drift = (slope.b+52*radius)/oriented.a
    assert lam.a-drift.b > broader.a and lam.b+drift.b < broader.b
    signed = signed_transfers(x,e,eps,widths,true_ranges,state_radius)
    # A 36-atom Pauli product frame has minimum singular value 1 for the HS norm.
    # Its pseudoinverse can correct each of the two source mixtures in this ball.
    rays = [prior.sym_matrix((np.eye(2)+s*p)/2) for p in
            (np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])) for s in (-1,1)]
    atoms = [sp.kronecker_product(a,b) for a in rays for b in rays]
    prep = record['preparation']
    for name, goal, perturb in [('kept',x[0]-x[1],2),('dephased',2*x[1],2)]:
        ws = [sp.Rational(float(w)) for w in prep[name]]
        err = goal-sum((w*p for w,p in zip(ws,atoms)),sp.zeros(4))
        assert ri(min(ws)) > iv.sqrt(ri(inner(err,err)))+perturb*radius
    return dict(neighborhood_HS_radius=1e-4, state_HS_error_upper=2e-7,
                necessary_log_error_floor_lower=float(floor.a), gradient_norm_upper=float(slope.b),
                volume_offset_upper=float(mismatch.b), minimum_volume_direction=float(oriented.a),
                geometric_curvature_bound=float(curvature.b), multiplier_drift_upper=float(drift.b),
                range_effect_error_upper=float(eps), gain_interval=[23.1042361,23.1042363],
                volume_interval=[2.763471835,2.763471839], signed_transfers=signed,
                global_survey=global_survey.verify(e, eps),
                matched_counterpreparation=matched_counterpreparation(e, eps, physical_volume, eta),
                calibrated_gain=str(physical_eta),calibrated_volume=str(physical_volume))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    body = json.dumps(verify(),indent=2,allow_nan=False)+"\n"
    if args.write:
        Path(__file__).with_name("certificate.json").write_text(body)
    else:
        print(body,end="")
