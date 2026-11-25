"""Complete capped radar head channels with integer projection and propagation."""
import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
import json
from math import factorial, isqrt, log, sqrt
from pathlib import Path

import numpy as np

from .. import native_precision, station
from ...material_clocks.measure import density_center

ROUTES = tuple(permutations(range(4),2))
TICKS, SUBSTEPS, DEGREE = 12, 8, 22
UNIT = 2**48
LIMB = 2**24


def channel_error(d, row):
    z = Fraction(row,32)
    assert (1-z/32)**(-32) < 4
    rounding = Fraction(DEGREE*4*(1+isqrt(2*d)), UNIT)
    tail = 2*z**(DEGREE+1)/factorial(DEGREE+1)
    column = TICKS*SUBSTEPS*(rounding+tail)
    assert column < Fraction(1,10**6)
    leg = 2*column  # four source columns in each leg
    composite = 2*leg + leg**2
    branches = 1 + TICKS + TICKS*(TICKS-1)//2
    return branches*(2*composite+composite**2) + Fraction(1,10**13)


def substep(h, state):
    """Taylor for exp(-i h/32), with a shared 2**48 denominator on both real parts."""
    power = state
    total = state.copy()
    for k in range(1,DEGREE+1):
        power = (h@power)//(32*k)
        if k%2:
            sign = 1 if k%4==1 else -1
            total[:,0::2] += sign*power[:,1::2]
            total[:,1::2] -= sign*power[:,0::2]
        elif k%4==2:
            total -= power
        else:
            total += power
    return total


def moment(amplitudes):
    """Exact superoperator numerator from all position Kraus matrices; denominator UNIT**2."""
    # shape (positions, output coin, input coin, real or imaginary)
    re, im = (amplitudes[...,i].reshape(-1,16) for i in (0,1))
    x = np.hstack((re,im))
    assert np.max(abs(x)) < 2*UNIT and len(x) <= 960
    hi, lo = x//LIMB, x%LIMB
    hh, hl, ll = hi.T@hi, hi.T@lo, lo.T@lo
    gram = hh.astype(object)*LIMB**2 + (hl+hl.T).astype(object)*LIMB + ll.astype(object)
    choi = (gram[:16,:16]+gram[16:,16:], gram[16:,:16]-gram[:16,16:])
    return np.stack([c.reshape(4,4,4,4).transpose(0,2,1,3).reshape(16,16) for c in choi])


def compose(a, b):
    return np.stack((a[0]@b[0]-a[1]@b[1], a[0]@b[1]+a[1]@b[0]))


def panel(h, ports):
    d = h.shape[0]
    row = int(abs(h).sum(axis=1).max())
    assert d==3840 and row<=43 and h.dtype==np.int64 and not (h-h.T).nnz
    row_error = channel_error(d,row)
    # Exact column norm <2 during the recursion. Taylor power/sum magnitudes
    # are <9*UNIT; these integer dot products and sums stay within int64.
    assert 43*9*UNIT < 2**62
    state = np.zeros((d,len(ROUTES),4,2),dtype=np.int64)
    for k,(i,_) in enumerate(ROUTES):
        state[4*ports[i]:4*ports[i]+4,k,:,0] = UNIT*np.eye(4,dtype=np.int64)
    complete = np.zeros((2,16,16),dtype=object)
    complete[0] = UNIT**2*np.eye(16,dtype=object)
    hits = {route:[] for route in ROUTES}
    all_ends = {route:[complete] for route in ROUTES}
    cumulative = {route:np.zeros_like(complete) for route in ROUTES}
    for tick in range(1,TICKS+1):
        flat = state.reshape(d,-1)
        for _ in range(SUBSTEPS):
            flat = substep(h,flat)
        state = flat.reshape(d,len(ROUTES),4,2)
        for k,route in enumerate(ROUTES):
            end = 4*ports[route[1]]
            at = state[end:end+4,k]
            arrival = moment(at[None])
            hits[route].append(arrival)
            cumulative[route] += arrival
            at.fill(0)
            survivor = moment(state[:,k].reshape(-1,4,4,2))
            all_ends[route].append(cumulative[route]+survivor)
    channels = []
    for route in ROUTES:
        out = (all_ends[route][TICKS]-cumulative[route])*UNIT**2
        for a,first in enumerate(hits[route],1):
            out += compose(all_ends[route[::-1]][TICKS-a],first)
        assert np.max(abs(out)) <= 2*UNIT**4
        converted = np.ldexp(out.astype(float),-192)
        channels.append(converted[0]+1j*converted[1])
    return np.stack(channels),row_error


@lru_cache(maxsize=1)
def channels():
    """Return [charge(-,+), parity(e,o), directed route, flat output, flat input], uniform error."""
    board = station.iso.pairing.formation.federation(1)
    ports = station.iso.pairing.formation.region(board)
    banks, errors = [], []
    for charge in (-1,1):
        by_parity = []
        for odd in (False,True):
            matrices, error = panel(native_precision.generator(board,ports,charge,odd),ports)
            by_parity.append(matrices)
            errors.append(error)
        banks.append(by_parity)
    return np.array(banks), float(np.nextafter(float(max(errors)),np.inf))


def run():
    bank, error = channels()
    center, _, prior_x, radius = density_center()
    obs = np.kron([[0,1],[1,0]],np.eye(2))
    # The stored center is the exactly normalized rational center rounded to doubles.
    # With 4x4 contractions and |superoperator entries|<=2, 1e-11 covers its
    # binary64 state conversion, matrix sums, observable trace and decimal deposit.
    offset = error + 2*np.sqrt(8)*float(radius)+1e-11
    changes = []
    for charge in bank:
        read = [(ch@x.ravel()).reshape(-1,4,4) for ch,x in zip(charge,center)]
        final = np.trace((read[0]+read[1])@obs,axis1=1,axis2=2).real
        changes.append(final-float(prior_x))
    changes = np.array(changes)
    assert changes.max()+offset < -.00051
    trace = np.eye(4).ravel()
    tp = np.max(abs(trace@bank-trace))
    assert tp < error
    starts = 2**31
    rng = np.random.default_rng(260024)
    p0 = (1+float(prior_x))/2
    sham = int(rng.binomial(starts,p0))
    counts = rng.binomial(starts,p0+changes/2)
    rate_radius = sqrt(log(2*(1+changes.size)/.01)/(2*starts))
    assert np.all((counts-sham)/starts + 2*rate_radius < 0)
    return dict(directed_routes=[list(r) for r in ROUTES],charges=[-1,1],tick_cap=TICKS,
                initial_left_X=float(prior_x),left_X_change_center=changes.tolist(),
                all_settings_left_X_change_upper=float(np.nextafter(changes.max()+offset,np.inf)),
                all_settings_left_X_change_lower=float(np.nextafter(changes.min()-offset,-np.inf)),
                per_setting_expectation_error_upper=float(np.nextafter(offset,np.inf)),
                channel_hermitian_trace_norm_error_upper=error,
                largest_trace_preservation_entry_residual=float(tp),
                simulated_left_X=dict(starts_per_setting=starts,seed=260024,
                    sham_plus=sham,pulse_plus=counts.tolist(),other_outcome="starts minus plus",
                    simultaneous_rate_radius_99_percent=rate_radius))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    result=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        Path(__file__).with_name("results.json").write_text(result)
    else:
        print(result,end="")
