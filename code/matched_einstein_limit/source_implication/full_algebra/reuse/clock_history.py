"""Information retained in pre-timeout clock bins of the same four-station survey."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import sympy as sp

from .. import certificate as c
from .. import native_precision as native
from .. import station
from . import invisible, measure

# Flattening index 12 * directed route index + (return tick-1); all chosen ticks precede 12.
BINS = (43,141,10,105,137,4,44,106,45,104,41,93,9,40,46,58,56,101)
ADDITIONAL = (93,105,22,136,104,42)
CHARGE_ROWS = (16,23,0,11,20,12,3,8,4,7,19,1,10,15,6,5,22,9)
F = native.prior.FRAC


def word(i):
    return [*measure.ROUTES[i//12],i%12+1]


def timing_error(last_tick=11):
    """Maximum bin-effect error, and the error on the usual capped time mean."""
    if last_tick not in (11, 12):
        raise ValueError("last_tick must be 11 or 12")
    rnd = Fraction(8,F)
    sampled = Fraction(1,10**19)+rnd
    first = [Fraction(0)]
    for t in range(1,12):
        first.append(sampled+sum(sampled+(1+sampled)*first[j]+rnd for j in range(1,t)))
    errors = []
    for t in range(1,last_tick+1):
        total = 0
        for a in range(1,t):
            e = first[a]+first[t-a]+first[a]*first[t-a]+rnd
            total += 2*e + e*e + rnd
        errors.append(total)
    return max(errors), sum(Fraction(12-t,4)*errors[t-1] for t in range(1,12))


def first_hits(samples,i,j):
    hits = [np.zeros((2,4,4),dtype=object)]
    for t in range(1,12):
        f = samples[:,t-1,j,i].copy()
        for a in range(1,t):
            f -= native.prior.multiply(samples[:,t-a-1,j,j],hits[a])
        hits.append(f)
    return hits


def native_events(last_tick=11):
    """Yield completed-return effects through the requested tick on the 12-tick survey."""
    if last_tick not in (11, 12):
        raise ValueError("last_tick must be 11 or 12")
    board = station.iso.pairing.formation.federation(1)
    ports = station.iso.pairing.formation.region(board)
    ids = np.array([4*p+j for p in ports for j in range(4)])
    coefficients = native.prior.weights()
    for qi,q in enumerate((-1,1)):
        for parity in (0,1):
            h = native.generator(board,ports,q,bool(parity))
            k = native.kernel_panels(h,ids,coefficients)
            hits = {r:first_hits(k,*r) for r in measure.ROUTES}
            for i,route in enumerate(measure.ROUTES):
                forward, back = hits[route], hits[route[::-1]]
                for t in range(2,last_tick+1):
                    effect = np.zeros((2,4,4),dtype=object)
                    for a in range(1,t):
                        amplitude = native.prior.multiply(back[t-a],forward[a])
                        effect += native.prior.multiply(native.prior.adjoint(amplitude),amplitude)
                    yield qi,parity,i,t,effect


def calibrated_panels():
    """Return integer coefficient numerators for pooled bins and per-charge means.

    Pooled bin rows have denominator 4F; individual charge mean rows have denominator 8F.
    """
    basis = 2*station.iso.basis(True)
    re,im = basis.real.astype(int),basis.imag.astype(int)
    assert np.array_equal(re+1j*im,basis)
    by_charge = np.zeros((2,12,11,18),dtype=object)
    for qi,parity,i,t,effect in native_events():
        by_charge[qi,i,t-1,9*parity:9*parity+9]=(re*effect[0]+im*effect[1]).sum(axis=(1,2))
    weights = np.arange(1,12)-12
    signed_means = (by_charge*weights[None,None,:,None]).sum(axis=2).reshape(24,18)
    return by_charge.sum(axis=0),signed_means


def square_certificate(matrix,op_error):
    """Bound the exact matrix's inverse in Euclidean coordinate norms."""
    # The center is rational. Its exact inverse eliminates BLAS proposal
    # variability and the corresponding center residual simultaneously.
    inverse = matrix.inv(method="DM")
    norm = invisible.frobenius(inverse)
    defect = op_error * norm
    assert defect.b < c.ri(".00001")
    bound = norm/(1-defect)
    return {"inverse_operator_norm_upper":float(np.nextafter(float(bound.b),np.inf)),
            "inverse_residual_operator_norm_upper":float(np.nextafter(float(defect.b),np.inf)),
            "smallest_singular_value_lower":float(np.nextafter(float((1/bound).a),-np.inf))}


def run():
    c.iv.dps=30
    bins, signed_means = calibrated_panels()
    eps_t, eps_m = timing_error()
    def get_rows(indices):
        return sp.Matrix([bins[i//12,i%12].tolist() for i in indices])/(4*F)
    square = get_rows(BINS)
    direct = square_certificate(square,c.ri(36*sp.Rational(eps_t)))
    assert direct["inverse_operator_norm_upper"] < 2e7
    # The original pooled means are half the sum of the two per-charge effects.
    means = sp.Matrix(signed_means)/(8*F)
    pooled = (means[:12,:]+means[12:,:])/2
    combined = square_certificate(pooled.col_join(get_rows(ADDITIONAL)),
                                  c.ri(36*sp.Rational(max(eps_t,eps_m))))
    split = square_certificate(means[list(CHARGE_ROWS),:],c.ri(36*sp.Rational(eps_m)))

    # Held source witness: the archived rational center is rounded to binary64;
    # 1e-14 HS covers all 18 coordinate conversions.
    record = json.loads(Path(invisible.__file__).with_suffix(".json").read_text())
    d = sp.Matrix([sp.Rational(x) for x in record["tangent_center_even_then_odd"]])
    d_error = c.ri(sp.Rational(record["true_kernel_correction_HS_upper"])+sp.Rational(1,10**14))
    tick = 12*measure.ROUTES.index((2,3))+9
    row = get_rows([tick]).T
    # Every one of its 18 Pauli/2 coefficients has native error <= 2 eps_t.
    row_error = c.ri(9*sp.Rational(eps_t))
    sigma = invisible.frobenius(row)*d_error+row_error*(invisible.frobenius(d)+d_error)
    slope = c.ri(row.dot(d))+c.iv.mpf([-1,1])*sigma
    assert slope.b < c.ri("-.000003")
    raw = [[str(int(j)) for j in bins[i//12,i%12]] for i in BINS]
    return dict(clock_tick_size="1/2",cap_ticks=12,kept_clock_bins=[word(i) for i in BINS],
                column_order="even then odd, each XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ, Pauli product/2",
                pooled_timestamp_matrix=dict(denominator=str(4*F),numerators=raw,
                    entry_error_upper=float(np.nextafter(float(2*eps_t),np.inf))),
                charge_hidden_pre_cap_clock_tomography=direct,
                original_means_plus_six_bins=dict(additional_bins=[word(i) for i in ADDITIONAL],**combined),
                per_charge_mean_tomography=dict(
                    settings=[[-1 if i//12==0 else 1,*measure.ROUTES[i%12]] for i in CHARGE_ROWS],**split),
                witness_event=dict(route_and_return_tick=word(tick),probability_slope=str(slope)))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    result=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        Path(__file__).with_suffix(".json").write_text(result)
    else:
        print(result,end="")
