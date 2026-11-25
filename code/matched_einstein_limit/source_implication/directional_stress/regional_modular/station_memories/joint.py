"""Signed entropy of the four terminal loads jointly with the two parked heads."""
import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from . import channel, total

region = total.region
c = total.c
R = sp.Rational
U = 2**55


def integer_matrix(matrix):
    parts=[x.as_real_imag() for x in matrix]
    denominator=int(sp.ilcm(*(p.q for item in parts for p in item)))
    return np.array([[int(z[k]*denominator) for z in parts] for k in (0,1)],dtype=object).reshape(2,*matrix.shape),denominator


def product(a,b):
    return np.stack((a[0]@b[0]-a[1]@b[1],a[0]@b[1]+a[1]@b[0]))


def adjoint(x):
    return np.stack((x[0].T,-x[1].T))


def norm(x,denom):
    return iv.sqrt(c.ri(R(int((x*x).sum()),int(denom)**2)))


def apply(block, state):
    effects, divisor=block
    x, scale=state
    real=(effects[0]*x[0]-effects[1]*x[1]).sum(axis=(-2,-1))
    imag=(effects[1]*x[0]+effects[0]*x[1]).sum(axis=(-2,-1))
    return np.stack((real,imag)), divisor*scale


def pairing(a,b):
    x,k=a
    y,l=b
    return c.ri(R(int((x*y).sum()),k*l))


def bound_log(state, map_error):
    x,d=state
    n=x.shape[1]
    f=x.astype(float)/d
    val,vec=np.linalg.eigh(f[0]+1j*f[1])
    v=np.stack((np.rint(vec.real*U),np.rint(vec.imag*U))).astype(np.int64).astype(object)
    lam=np.rint(val*U).astype(np.int64).astype(object)
    assert (lam>0).all(), val
    logs=np.rint(-np.log(np.asarray(lam,dtype=float)/U)*U).astype(np.int64).astype(object)
    assert (lam>0).all() and all(abs(iv.ln(c.ri(R(y,U)))+c.ri(R(z,U))).b < c.ri("3e-13")
                                   for y,z in zip(lam,logs,strict=True))
    identity=np.stack([np.eye(n,dtype=object),np.zeros((n,n),dtype=object)])
    defect=norm(product(adjoint(v),v)-identity*U**2,U**2)
    reconstructed=product(v*lam,adjoint(v))
    residual=norm(x*U**3-reconstructed*d,d*U**3)
    assert defect.b<c.ri("2e-13") and residual.b<c.ri("2e-13")
    assert max(lam)<U and max(logs)<12*U
    discrepancy=residual+3*defect+map_error
    lower=c.ri(R(min(lam),U))-discrepancy
    assert lower.a>0
    # The polar unitary differs from v by <=defect; conjugating a matrix
    # of norm <=12 incurs at most 36*defect. Log is 1/lower Lipschitz here.
    error=discrepancy/lower+36*defect+c.ri("3e-13")
    return (product(v*logs,adjoint(v)),U**3), error, lower


def run():
    iv.dps=32
    ref,rad,odds=region.measure.source()
    direction,corr=region.controls()["fixed_all_number_weights"]
    prop=region.measure.bridge.native.prior
    station=region.measure.bridge.native.iso.pairing.formation
    board=station.federation(1)
    route=station.region(board)
    rings=region.measure.bridge.beam.screens(board,route)
    weights=prop.weights()[7]
    paulis=[c.prior.sym_matrix(a) for a in c.iso.basis(False)]
    axes=[integer_matrix(v) for v in paulis]
    # A compressed four-column whole-flight operator has error <1.2e-10.
    # For any Hermitian input the induced trace-norm error of A X A* is
    # <= (2+epsilon)*epsilon*||X||_1. Averaging preserves the bound.
    map_errors=[c.ri("1e-18"),c.ri("1e-18"),c.ri("2.5e-10")]
    score=err=gradient_size=curvature=0
    conditional=[]
    for odd,(x,d) in enumerate(zip(ref,direction,strict=True)):
        maps=channel.channels(board,route,rings,bool(odd),weights)
        center,gain,dlogs=0,0,[]
        gradient=[c.ri(0) for _ in axes]
        cc=integer_matrix(x)
        change=integer_matrix(d)
        nd=region.norm(d)
        output_hs_lower=[]
        for ch,ep in zip(maps,map_errors,strict=True):
            modular, log_error, least=bound_log(apply(ch,cc),ep)
            dlogs.append(log_error)
            deriv=apply(ch,change)
            center+=pairing(modular,deriv)
            gain+=ep*(norm(*modular)+log_error)
            gradient=[g+pairing(modular,apply(ch,a)) for g,a in zip(gradient,axes,strict=True)]
            # CP trace-norm contraction for each outcome block. Each input
            # has trace norm at most twice its HS norm.
            out_error=2*c.ri(corr)+2*nd*ep
            output_hs_lower.append(max(c.ri(0),norm(*deriv).a-out_error))
        gradient_size+=sum(v*v for v in gradient)
        _, _, spec=region.logarithm(x, c.ri(0))
        # A traceless 4x4 error has op norm <=sqrt(3/4) times HS norm.
        relative=iv.sqrt(c.ri(R(3,4)))*c.ri(rad)/spec[0].a
        # Positive maps preserve the two Loewner inequalities relative to x.
        # Log monotonicity and trace-preserving tracelessness permit subtraction
        # of the midpoint of [log(1-relative),log(1+relative)] from all blocks.
        source_error=(iv.ln(1+relative)-iv.ln(1-relative))/2
        total_error=2*(nd+c.ri(corr))*(source_error+max(v.b for v in dlogs)+gain)
        score+=center
        err+=total_error
        # For the physical block of weight w, eigenvalues are <=w at every
        # feasible point; d^2 S <= -||d output||_HS^2/w.
        weight=odds if odd else 1-odds
        curvature+=sum(v.a**2 for v in output_hs_lower)/weight.b
        conditional.append(dict(odd=bool(odd), nominal_slope=str(center), conditional_source_error=str(source_error),
                                maximum_log_error=str(max(v.b for v in dlogs))))
    err+=iv.sqrt(gradient_size)*c.ri(corr)
    slope=score+total.SYMMETRIC*err
    assert slope.a>c.ri("1.37e-6"), str(slope)
    t=c.ri(R(1,100))
    finite=(t*slope-t*t*curvature/2).b
    assert finite<0
    return dict(control="fixed_all_number_weights", charge_erased=True, retained_ports=list(route),
                retained_heads=2, fair_charge_weights=["1/2","1/2"], reference_slope=str(slope),
                nominal_slope=str(score), entropy_change_at_0p01_upper=str(finite),
                quadratic_entropy_curvature_lower=str(curvature.a),
                centered_tangent_error_gain=str(iv.sqrt(gradient_size).b), conditionals=conditional)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    result=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        (total.HERE/"joint.json").write_text(result)
    else:
        print(result,end="")
