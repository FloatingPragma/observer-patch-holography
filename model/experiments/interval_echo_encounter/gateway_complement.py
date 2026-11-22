"""Resolve the finite-time complementary-band contribution to each receiving aperture."""
import argparse
from decimal import Decimal, localcontext
import hashlib
from math import factorial, log

import numpy as np
from scipy import sparse
from scipy.special import jv

from .gateway import RECEIPT, close, load_json, read, registered_receipt, save
from .gateway_source import (HERE, PairSector, candidate, federation, geometry, high_projection,
                             reference, require, settings, source)

PROTOCOL = HERE/"gateway_complement_protocol.json"
OUTPUT = HERE/"gateway_complement.json"


def coefficients():
    # Jacobi-Anger series for exp(-i*14*pi*x), |x| <= 1.
    with localcontext() as ctx:
        ctx.prec = 90
        z = 14*Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406")
        half = z/2
        values = []
        for k in range(101):
            term = half**k/Decimal(factorial(k))
            value = term
            for j in range(1,150):
                term *= -(half*half)/(j*(k+j))
                value += term
            require(abs(term)<Decimal("1e-60"), "Bessel series tail too large")
            values.append(float(value))
    require(np.max(abs(values-jv(np.arange(101),14*np.pi)))<1e-13, "Bessel weights disagree")
    w = np.asarray(values)*(-1j)**np.arange(101)
    w[1:] *= 2
    require(sum(abs(w))<10 and 2*22**101/factorial(101)/(1-22/102)<1e-21,
            "Chebyshev tail or coefficient budget exceeded")
    return w


def step(x, state, weights):
    prev = state
    curr = x @ prev
    total = weights[0]*prev + weights[1]*curr
    for coeff in weights[2:]:
        prev,curr = curr, 2*(x@curr)-prev
        total += coeff*curr
    return total


def apertures(state, receiver):
    local = abs(state[receiver,:])
    return np.vstack((local, np.linalg.norm(local,axis=0)))


def replay(record=False):
    cfg = settings()
    receipt = registered_receipt(load_json(RECEIPT))
    sector = PairSector(federation(cfg["level"]))
    layout = geometry(sector.car)
    h,admit,prepared = source(sector,candidate.KAPPA,layout)
    inds = [sector.index(*e) for e in layout[-1]]
    high_source,error = high_projection(h,prepared)
    initial = apertures(high_source,inds)
    # Evolve both complete launches independently of the sparse exponential
    # that supplied the counts. Subtract amplitudes before taking any norms.
    state = np.column_stack((high_source, prepared))
    norms = np.linalg.norm(state,axis=0)
    x = ((h+2*sparse.eye(sector.P,format="csr"))/14).tocsr()
    require(max(np.asarray(abs(x).sum(axis=1)).flat)<1+1e-14 and max(np.diff(x.indptr))<=13,
            "propagation exceeds arithmetic budget")
    weights = coefficients()
    L = log(2*(32+31)*81/cfg["alpha"])
    results = {}
    bookings = {tick:word for word,tick in reference(cfg["appointments"]).items()}
    for tick in range(max(bookings)+1):
        if tick:
            state = step(x,state,weights)
        if tick not in bookings:
            continue
        word = bookings[tick]
        require(max(abs(np.linalg.norm(state,axis=0)-norms))<1e-9, "independent evolution lost norm")
        upper = apertures(state[:, :2],inds) + error + tick*1e-8
        total = apertures(state[:, 2:],inds)
        direct_low = apertures(state[:, 2:] - state[:, :2],inds)
        lo,hi = [],[]
        for j,sign in enumerate(("plus","minus")):
            key = f"native_{sign}_{word}"
            p = np.array(receipt["weights"][key])/10**13
            require(abs(sum(p[:-1])-admit)<3.2e-12,"source enrollment differs")
            require(np.sum(abs(admit*total[:-1,j]**2 - p[:30])) < 3.2e-12,
                    "independent complete propagation disagrees with the sampling law")
            lower,higher = read(receipt["counts"][key],cfg["shots"],p,L)
            lo.append(lower)
            hi.append(higher)
        low = np.maximum(0,np.sqrt(np.asarray(lo).T)-upper)
        high = np.minimum(1,np.sqrt(np.asarray(hi).T)+upper)
        movement = np.maximum(0,low-initial-error)
        contrast = np.maximum(0,np.maximum(low[:,0]-high[:,1],low[:,1]-high[:,0]))
        direct_error = error + 2*tick*1e-8
        require(np.all((low <= direct_low + direct_error) & (direct_low <= high + direct_error))
                and np.all(contrast <= abs(direct_low[:, 0] - direct_low[:, 1]) + 2*direct_error),
                "local bounds exclude the directly evolved bound band")
        results[str(word)] = {"upper_amplitude_limit":upper.tolist(), "motion_lower":movement.tolist(),
                              "phase_norm_lower":contrast.tolist()}
    files = [PROTOCOL, HERE/"gateway_complement.py", RECEIPT]
    out = {"sources":{p.name:hashlib.sha256(p.read_bytes().replace(b"\r\n",b"\n")).hexdigest()
                      for p in files}, "bookings":results}
    if record:
        save(OUTPUT,out)
    else:
        close(out,load_json(OUTPUT))
    print("Verified complementary amplitudes, motion and phase bounds in 31 apertures at nine bookings.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    replay(parser.parse_args().record)
