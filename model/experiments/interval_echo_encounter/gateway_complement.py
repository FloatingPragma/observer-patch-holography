"""Resolve the finite-time complementary-band contribution to each receiving aperture."""
import argparse
from decimal import Decimal, localcontext
import hashlib
from math import factorial, log

import numpy as np
from scipy import sparse
from scipy.special import jv

from .gateway import RECEIPT, close, load_json, pins, read, save
from .gateway_source import (HERE, PairSector, candidate, federation, geometry, high_projection,
                             require, settings, source)

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
    receipt = load_json(RECEIPT)
    require(receipt["sources"] == pins(), "complement witness needs the registered launches")
    sector = PairSector(federation(cfg["level"]))
    layout = geometry(sector.car)
    h,admit,prepared = source(sector,candidate.KAPPA,layout)
    inds = [sector.index(*e) for e in layout[-1]]
    state,error = high_projection(h,prepared)
    initial = apertures(state,inds)
    norms = np.linalg.norm(state,axis=0)
    x = ((h+2*sparse.eye(sector.P,format="csr"))/14).tocsr()
    require(max(np.asarray(abs(x).sum(axis=1)).flat)<1+1e-14 and max(np.diff(x.indptr))<=13,
            "propagation exceeds arithmetic budget")
    weights = coefficients()
    L = log(2*(32+31)*81/cfg["alpha"])
    results = {}
    for word in range(cfg["appointments"][-1]+1):
        if word:
            state = step(x,state,weights)
        if word not in cfg["appointments"]:
            continue
        require(max(abs(np.linalg.norm(state,axis=0)-norms))<1e-9, "upper-band evolution lost norm")
        upper = apertures(state,inds) + error + word*1e-8
        lo,hi = [],[]
        for sign in ("plus","minus"):
            key = f"native_{sign}_{word}"
            p = np.array(receipt["weights"][key])/10**13
            require(abs(sum(p[:-1])-admit)<3.2e-12,"source enrollment differs")
            lower,higher = read(receipt["counts"][key],cfg["shots"],p,L)
            lo.append(lower)
            hi.append(higher)
        low = np.maximum(0,np.sqrt(np.asarray(lo).T)-upper)
        high = np.minimum(1,np.sqrt(np.asarray(hi).T)+upper)
        movement = np.maximum(0,low-initial-error)
        contrast = np.maximum(0,np.maximum(low[:,0]-high[:,1],low[:,1]-high[:,0]))
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
