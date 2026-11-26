"""Simulate complete terminal records for the frozen pentagon experiment."""
import argparse
import json
from pathlib import Path

import numpy as np

from . import apparatus as lab

HERE = Path(__file__).resolve().parent
SHOTS = 131072
SEED = 271216
GRID = 2**36


def sample(rng, probs):
    p = np.asarray(probs)
    if min(p) < -1e-12 or abs(p.sum()-1.) > 1e-11:
        raise ArithmeticError("Terminal outcome probabilities do not sum to one")
    cuts = np.rint(np.cumsum(p[:-1])*GRID).astype(np.uint64)
    values = rng.bit_generator.random_raw(SHOTS) >> np.uint64(28)
    cumulative = [0, *(int(np.count_nonzero(values < edge)) for edge in cuts), SHOTS]
    return np.diff(cumulative).tolist()


def population_probs(state, owner, port):
    # Mutually exclusive sharp ports. The other-cell outcome is retained.
    ring, anti, _ = lab.geometry()
    local = state[12*owner:12*(owner+1)]
    launch = ring[port]
    back = [anti[i] for i in launch]
    ps = abs(local)**2
    front, opposite = ps[launch].sum(), ps[back].sum()
    return [front,opposite,ps[port]+ps[anti[port]], max(0.,1.-ps.sum())]


def acquire():
    rng = np.random.default_rng(SEED)
    lab.sharp_calibration()
    edges = lab.cell().seams
    calibration = []
    for a,b in edges:
        for sign in (-1,1):
            target = np.zeros(12,complex)
            target[[a,b]] = [1/np.sqrt(2), sign*1j/np.sqrt(2)]
            wave,_ = lab.prepare(lab.SimpleNamespace(car=lab.cell(),active=edges),target,root=a)
            # The first detector is the half-exchange current quadrature. The
            # second keeps only this same physical bond active for a short dwell.
            # At +π/6 its right-minus-left occupancy is sin(π/6) times the signed input.
            read = lab.endpoint_probabilities(-1,0,wave,(a,b),False)
            t = np.pi/6
            _,h = lab.local_source(1,lab.candidate.KAPPA)
            from scipy.linalg import expm
            later = abs(expm(-1j*t*h) @ wave[[a,b]])**2
            calibration.append([a,b,sign,sample(rng,read),sample(rng,[*later,0.])])
    ring_rows=[]
    populations=[]
    def panel(level,owner,port,m,t,sorted_parity=True):
        wave = lab.evolve(level,owner,port,m,t)
        local = wave[12*owner:12*(owner+1)]
        entries=[sample(rng,lab.endpoint_probabilities(level,owner,local,e,sorted_parity)) for e in edges]
        ring_rows.append(dict(level=level,owner=owner,source=port,winding=m,time=t,
                              parity=sorted_parity,edges=entries))
        if level >= 0:
            populations.append(dict(level=level,owner=owner,source=port,winding=m,time=t,
                                    rate=1,counts=sample(rng,population_probs(wave,owner,port))))
    for port in range(12):
        for m in (1,2):
            panel(-1,0,port,m,1)
    for m in (-1,-2,0):
        panel(-1,0,0,m,1)
    for m in (1,2):
        panel(-1,0,0,m,0,False)
        panel(-1,0,0,m,3)
    for level in (0,1):
        for owner in (0,13):
            for m in (1,2):
                for t in (1,3):
                    panel(level,owner,0,m,t)
    for rate in (1,0.5):
        for m in (1,2,-1,-2):
            for t in range(7):
                wave=lab.evolve(-1,0,0,m,t,rate)
                populations.append(dict(level=-1,owner=0,source=0,winding=m,time=t,
                                        rate=rate,counts=sample(rng,population_probs(wave,0,0))))
    return dict(seed=SEED,shots=SHOTS,bond_labels=[list(e) for e in edges],
                cycles=lab.geometry()[0],antipodes=list(lab.geometry()[1]),
                quadrature_and_short_flight=calibration,circulation=ring_rows,populations=populations)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=acquire()
    path=HERE/"records.json"
    if args.write:
        path.write_text(json.dumps(result, separators=(",",":"))+"\n")
    elif result != json.loads(path.read_text()):
        raise ValueError("Port current records differ")
    from .readback import decode
    print(json.dumps(decode(result),indent=2))


if __name__ == "__main__":
    main()
