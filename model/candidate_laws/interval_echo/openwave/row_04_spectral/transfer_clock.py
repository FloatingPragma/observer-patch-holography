"""Drive the donor's original 32 analyzers from a separately measured local hopping clock."""

import argparse
from hashlib import sha256
import json
from math import asin, log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from tools.strict_json import load
from .experiment import HERE, compare
from .graph import DONOR, History, donor_selection, groups
from .local import check_local, require
from .source import install, triangles


SHOTS = 2**39
T = np.exp2(np.arange(12))/256
TOL = sqrt(log(2*len(T)*17/.01)/(2*SHOTS)) + 5.1e-14


def infer(counts):
    """No strength, state, fit to a release, or target spectrum enters the reader."""
    if len(counts)!=len(T) or any(type(c) is not int or not 0<=c<=SHOTS for c in counts):
        raise ValueError("A complete clock bank is required")
    lo,hi = 0.,128.
    for t,c in zip(T,counts,strict=True):
        if t*hi/2>pi/2:
            return None  # More than one fringe can explain this appointment.
        p = c/SHOTS
        inv = lambda z: 2*asin(sqrt(max(0.,min(1.,z))))/t
        lo,hi = max(lo,inv(p-TOL)),min(hi,inv(p+TOL))
        if lo>hi:
            return None
        if p>1/8:
            return float(min(hi,max(lo,inv(p)))),float(lo),float(hi),float(t)
    return None


def targets(tree):
    vals,v = eigh(tree.lap)
    columns = []
    for band in groups(vals):
        q = v[:,band] @ v[:,band].T
        basis = []
        for w in q.T:
            y = w.copy()
            for old in basis:
                y -= old*(old @ y)
            if np.linalg.norm(y)>1e-7:
                basis.append(y/np.linalg.norm(y))
        require(len(basis)==len(band),"Donor analyzer dimensions differ")
        columns.extend(basis)
    return np.column_stack(columns)


def run_word(tree,angles,step,ratio):
    n = len(tree.paths)
    phases = np.exp(-1j*ratio*angles)
    initial = np.zeros((n,n),complex)
    initial[0,:] = 1
    U = expm(-.5j*step*8*ratio*tree.twice_scaled)
    for i in range(phases.shape[1]):
        initial = (U @ initial)*phases[:,i,:].T
    np.testing.assert_allclose(np.sum(abs(initial)**2,axis=0),1,rtol=0,atol=4e-12)
    return initial


def calibration_law(k):
    # Read the opposite port's occupation after an isolated exposure, on the
    # actual sign-even written edge. Only the two eligible columns are occupied.
    law = candidate.native_law
    H0,H1 = law.local_term_parts()
    modes = np.zeros((law.DIM,2))
    for j,(coin,w) in enumerate(((1,(1,0)),(0,(0,1)))):
        for d in (-1,1):
            modes[(1+law.LABELS.index((coin,d)))*law.N_SRC+law.SRC.index(w),j] = 1/np.sqrt(2)
    h = (H0+k*H1) @ modes
    sub = modes.T @ h
    np.testing.assert_allclose(h,modes @ sub,rtol=0,atol=2e-13)
    prob = np.array([abs(expm(-1j*t*sub)[1,0])**2 for t in T])
    np.testing.assert_allclose(prob,np.sin(k*T/2)**2,atol=3e-13,rtol=0)
    # The simulated sampler uses 13-place probabilities. The deterministic
    # rounding allowance is included in the clock interval above.
    return np.round(np.sin(k*T/2)**2,13)


def run():
    check_local()
    b = federation(0)
    tree = History(b,install(b,triangles(b)[0])[0])
    select = donor_selection(tree)["choices"]
    old = load(DONOR/"analyzers.json")
    masks = np.array([tree.sites==p for p in old["ports"]],dtype=int)
    angles = np.array([s["angles"] for s in old["pulses"]]) @ masks
    ideal = targets(tree)
    baseline = run_word(tree,angles,old["step"],1)
    def fidelity(state):
        outside = state-ideal*(ideal.conj()*state).sum(axis=0)
        return 1.-np.sum(abs(outside)**2,axis=0)
    # Bound non-target amplitude rather than subtracting two nearly equal
    # probabilities. Common scalar phase is irrelevant at negative strength.
    resid = np.linalg.norm(baseline-ideal*(ideal.conj()*baseline).sum(axis=0),axis=0)
    costs = old["layers"]*8*old["step"]*max(abs(np.linalg.eigvalsh(tree.twice_scaled)))/2 + np.max(angles,axis=2).sum(axis=1)
    strengths = [s*k for k in (.5,1,2,4,8,16,64,128) for s in (1,-1)] + [0]
    rows = []
    for i,k in enumerate(strengths):
        count = [int(x) for x in np.random.default_rng([44130,i]).binomial(SHOTS,calibration_law(k))]
        read = infer(count)
        row = {"strength_in_simulator":k,"counts":count,"calibration":read}
        if read is None:
            require(k==0,"Nonzero clock not resolved")
        else:
            value,lo,hi,_ = read
            require(lo<=abs(k)<=hi,"Clock interval missed the planted strength")
            error = max(abs(lo/value-1),abs(hi/value-1))*costs
            bounds = 1-np.minimum(1.,resid+error+1e-8)**2
            actual = run_word(tree,angles,old["step"],k/value)
            observed = fidelity(actual)
            require(np.min(bounds)>.98 and np.all(observed+1e-11>=bounds),
                    "A transferred spatial analyzer failed the registered fidelity")
            later = expm(-.5j*37*k*tree.twice_scaled) @ actual
            drift = max(np.max(abs(np.sum(abs(actual)**2,axis=0)-1)),
                        np.max(abs(np.sum(abs(later)**2,axis=0)-1)))
            require(drift<5e-12, "Phase propagation exceeds its norm allowance")
            row.update(fidelity=observed.round(10).tolist(),maximum_state_error=float(max(error)),
                       minimum_bound=float(min(bounds)),maximum_column_norm_squared_drift=float(drift),returns_37=[
                           float(1.-np.sum(abs(later[~np.isin(tree.sites,[p[-1] for p in c["paths"]]),c["analyzer"]])**2))
                           for c in select])
        rows.append(row)
    return {"schema":"oph.interval_echo.spectrometer_hop_clock.v1","shots_per_appointment":SHOTS,
            "protocol_sha256":sha256((HERE/"CLOCK_TEST.md").read_bytes().replace(b"\r\n",b"\n")).hexdigest(),
            "confidence_radius":TOL,"appointments":T.tolist(),"runs":rows}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--record",type=Path)
    args=ap.parse_args()
    result=json.loads(json.dumps(run(),allow_nan=False))
    if args.record:
        with args.record.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,separators=(",",":"),allow_nan=False)
            f.write("\n")
    else:
        compare(result,load(HERE/"clock.json"))
    print("16 nonzero strengths, 32 original words each; lowest simultaneous fidelity bound:",
          min(r["minimum_bound"] for r in result["runs"] if r["calibration"] is not None))
