"""An entropy-increasing native tangent invisible to all twelve mean radar ranges."""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from .. import certificate as c
from .. import reuse_equilibrium, station
from . import measure


def frobenius(m):
    return c.iv.sqrt(c.ri(sum(x*x for x in m)))


def dyadic(m):
    """Use the numerical linear algebra as a proposal only."""
    return sp.Matrix(m.shape[0],m.shape[1],lambda i,j:
                     sp.Rational(round(float(m[i,j])*2**38), 2**38))


def stacked_error():
    # Existing per-channel bound is 79 (2*d+d^2) plus the final conversion,
    # where d bounds the operator error of *each* terminal Kraus branch.
    conversion = sp.Rational(1,10**13)
    budget = sp.Rational(measure.channel_error(3840,43))
    per_branch = (budget-conversion)/(2*79)
    # The direct sum of all terminal branches has error <= sqrt(79)*per_branch.
    difference = 9*per_branch
    return 2*difference + difference**2 + conversion


def tomography(r, rerror):
    # These six parity-blind mean products include the witness's YY control.
    # A nondemolition parked parity read selects the left transfer time. In both
    # branches the excess-to-rail gate is the same -ZZ SWAP, including its phase.
    bell = np.array([0,1,1,0])
    twice_projector = 2*np.eye(4,dtype=int)-np.outer(bell,bell)
    swap = np.eye(4,dtype=int)[[0,2,1,3]]
    zz = np.diag([1,-1,-1,1])
    for sign, rate in ((1,8),(-1,9)):
        blocks = station.iso.ranging.parity_blocks(1,8,sign)
        h = blocks[1:,1:].transpose(0,2,1,3).reshape(4,4)
        np.testing.assert_array_equal(2*h,rate*twice_projector)
    np.testing.assert_array_equal(np.eye(4)-twice_projector,-zz@swap)
    station.phase.rail_analyzers(1)
    labels = ("XX","XY","XZ","YX","YY","YZ","ZX","ZY","ZZ")
    chosen = ("XZ","YX","YY","YZ","ZY","ZZ")
    extra = sp.zeros(6,18)
    for i,label in enumerate(chosen):
        j = labels.index(label)
        extra[i,j] = extra[i,j+9] = 2
    augmented = r.col_join(extra)
    inverse = dyadic(np.linalg.inv(np.array(augmented,dtype=float)))
    inorm = frobenius(inverse)
    defect = frobenius(sp.eye(18)-augmented*inverse)+rerror*inorm
    assert defect.b < c.ri(".000003")
    lower = (1-defect)/inorm
    assert lower.a > c.ri("9.7e-9")
    return dict(parity_blind_products=list(chosen),minimal_additional_scalar_means=6,
                left_excess_transfer_dwell_even_odd=["pi/8","pi/9"],right_excess_transfer_dwell="pi/8",
                augmented_inverse_defect_upper=float(np.nextafter(float(defect.b),np.inf)),
                singular_value_lower=float(np.nextafter(float(lower.a),-np.inf)),
                HS_reconstruction_error_per_l2_readout_error_upper=float(
                    np.nextafter(float((1/lower).b),np.inf)))


def verify(bank):
    c.iv.dps=30
    error = stacked_error()
    upper = float(np.nextafter(float(error),np.inf))
    ops = [c.prior.sym_matrix(h) for h in station.iso.basis(False)]
    centers = [reuse_equilibrium.sector_center(bank[:,b].reshape(24,16,16),upper,ops) for b in (0,1)]
    radius = max(r for _,r in centers)
    # Orthonormal, weight-preserving, independent two-coin correlation directions in each block.
    corr = [c.prior.sym_matrix(h) for h in station.iso.basis(True)]
    assert sp.Matrix([[c.inner(a,b) for a in corr] for b in corr])==sp.eye(9)
    for m in corr:
        assert all(sum(m[2*i+b,2*j+b] for b in (0,1))==0 for i in (0,1) for j in (0,1))
        assert all(sum(m[2*b+i,2*b+j] for b in (0,1))==0 for i in (0,1) for j in (0,1))
    nums, unit, native_error = c.exact_effects()
    effects = c.rational_effects(nums,unit)
    r = sp.Matrix([[c.inner(b,m) for b in pair for m in corr] for pair in effects])
    # Each exact range coefficient has error <= ||Pauli/2||_1 * native_error.
    rerror = c.ri(30*sp.Rational(native_error))   # 2 sqrt(12*18) < 30
    nominal, log_errors = [], []
    for state, _ in centers:
        vals, bounded = c.prior.entropy_coordinates(state,0,sp.Rational(".18"),corr)
        nominal.extend(vals)
        log_errors.extend([2*bounded]*len(corr))
    g = sp.Matrix(nominal)

    # Propose a common right inverse and an entropy-increasing coordinate vector.
    raw = np.array(r,dtype=float)
    _,_,v = np.linalg.svd(raw,full_matrices=True)
    kernel = v[12:]
    proposed = kernel.T@(kernel@np.array(g,dtype=float).ravel())
    d = dyadic((proposed/np.linalg.norm(proposed))[:,None])
    k = dyadic(np.linalg.pinv(raw))
    d -= k*(r*d)   # exact rational nominal residual correction
    knorm = frobenius(k)
    defect = frobenius(sp.eye(12)-r*k)+rerror*knorm
    assert defect.b < c.ri(".000001")
    dnorm = frobenius(d)
    # True R K is invertible by Neumann. D=d-K(RK)^(-1)R d lies in ker R exactly.
    adjust = knorm*(frobenius(r*d)+rerror*dnorm)/(1-defect)
    assert adjust.b < c.ri(".0000004")
    true_norm = dnorm+adjust
    assert true_norm.b < c.ri("1.0001")

    # The HS log error of either conditional state is at most radius/.18 by the
    # resolvent identity; no sqrt(dimension) cost is needed on an HS perturbation.
    # The scalar logs of the two fixed parity weights cancel on every tangent.
    gerr = frobenius(sp.Matrix(log_errors)) + c.iv.sqrt(2)*c.ri(radius)/c.ri(".18")
    score = c.ri(g.dot(d))
    uncertainty = frobenius(g)*adjust+gerr*true_norm
    slope = score+c.iv.mpf([-1,1])*uncertainty
    assert slope.a > c.ri(".00001"), str(slope)

    # Certify native lottery slack at every p in the separately certified volume bracket.
    p=c.iv.mpf([".134",".1346"])
    moments=[np.array([c.ri_band(2*c.inner(state,op),2*rad) for op in ops],dtype=object)
             for state,rad in centers]
    slack=[(1-2*p)-sum(abs(m) for m in (moments[0]*(1-p)-moments[1]*p)),
           2*p-sum(abs(m) for m in (moments[1]*(2*p)))]
    assert all(s.a>c.ri(".20") for s in slack)
    # A displacement up to 0.01 changes at most nine correlation coefficients of
    # each unnormalized Pauli frame; their l1 changes are bounded by 12*t*||D||_2.
    assert c.ri(".20")-12*c.ri(".01")*true_norm.b > 0
    step=c.ri(".0000001")
    # Each normalized block has min eigenvalue >= .19-radius >.18. Including the
    # smaller parity weight, the whole chord has min eigenvalue above .024.
    assert c.ri(".134")*c.ri(".18")-step*true_norm.b > c.ri(".024")
    gain=step*slope.a-step**2*true_norm.b**2/(2*c.ri(".024"))
    assert gain.a > c.ri("7e-13")
    # YY is the fifth of the nine correlation products (XX,XY,XZ,YX,YY,...);
    # a parity-blind product read sums both unnormalized block correlations.
    yy = 4
    pauli_y = sp.Matrix([[0,-sp.I],[sp.I,0]])
    assert 2*corr[yy] == sp.kronecker_product(pauli_y,pauli_y)
    yy_slope = c.ri(2*(d[yy]+d[9+yy]))+c.iv.mpf([-1,1])*2*c.iv.sqrt(2)*adjust
    assert yy_slope.a > c.ri("1.53")
    center_words=[]
    for state,_ in centers:
        values=[c.inner(state,op)*2**35 for op in ops]
        assert all(v.is_Integer for v in values)
        center_words.append([int(v) for v in values])
    return dict(range_rows=12,correlation_coordinates=18,exact_range_kernel_dimension=6,
                channel_trace_error_upper=float(np.nextafter(float(error),np.inf)),
                conditional_HS_error_upper=float(np.nextafter(float(radius),np.inf)),
                right_inverse_defect_upper=float(np.nextafter(float(defect.b),np.inf)),
                true_kernel_correction_HS_upper=float(np.nextafter(float(adjust.b),np.inf)),
                entropy_derivative_nats=str(slope),physical_step="1/10000000",
                parity_blind_joint_YY_derivative=str(yy_slope),
                finite_entropy_increase_nats_lower=float(np.nextafter(float(gain.a),-np.inf)),
                matched_odd_weight=[".134",".1346"],lottery_slack_lower=".20",
                fixed_marginal_tomography=tomography(r,rerror),
                fixed_state_centers=dict(denominator=2**35,even_then_odd=center_words,
                    basis="IX,IY,IZ,XI,XX,XY,XZ,YI,YX,YY,YZ,ZI,ZX,ZY,ZZ; Pauli product /2; identity I/4"),
                tangent_center_even_then_odd=[float(i) for i in d])


def run():
    return verify(measure.channels()[0])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    report=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        Path(__file__).with_suffix(".json").write_text(report)
    else:
        print(report,end="")
