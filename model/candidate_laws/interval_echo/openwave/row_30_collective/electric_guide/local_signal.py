"""Local dimer-to-dimer intervention on the static guarded guide."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from . import measure


HOLDS = (0, .25, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64)
SHOTS = 2**26
CONTROLS = ("native", "cut", "no_flip")


def response(n, source, j, h, arm):
    c = 2*(h*np.eye(n)-j*np.eye(n,k=-1))
    if arm=="cut":
        c[source,source-1] = c[source+1,source] = 0
    elif arm=="no_flip":
        np.fill_diagonal(c,0)
    elif arm!="native":
        raise ValueError("Unknown physical intervention")
    u, e, v_transpose = np.linalg.svd(c)
    v = v_transpose.T
    result = []
    for t in HOLDS:
        co, si = np.cos(e*t), np.sin(e*t)
        aa = u @ (co*u[source])
        bb = v @ (co*v[source])
        ab = -u @ (si*v[source])
        ba = v @ (si*u[source])
        result.append(aa*bb-ab*ba)
    return np.array(result)


def acquire():
    board = measure.federation(2)
    _, path, _ = measure.electric.guide(board)
    n = 2*((len(path)-2)//4)
    rng = np.random.default_rng(302022)
    j, h = measure.electric.couplings(8)
    banks = {arm:measure.signs(response(n,n//2,j,h,arm),SHOTS,rng).tolist() for arm in CONTROLS}
    return dict(sha256=hashlib.sha256((measure.HERE/"LOCAL.md").read_bytes()).hexdigest(),
                shots=SHOTS, route=path[:2*n+2], holds=list(HOLDS), counts=banks)


def read(data):
    route = measure.electric.guide(measure.federation(2))[1]
    if (data["holds"]!=list(HOLDS) or data["route"]!=route[:90] or data["shots"]!=SHOTS or
            set(data["counts"])!=set(CONTROLS) or
            data["sha256"]!=hashlib.sha256((measure.HERE/"LOCAL.md").read_bytes()).hexdigest()):
        raise ValueError("Wrong local enrollment")
    n = 44
    radius = float(np.sqrt(2*np.log(2*3*n*len(HOLDS)/.01)/SHOTS))
    rows = []
    for arm, bank in data["counts"].items():
        readings = measure.read_sign(bank,(len(HOLDS),n),SHOTS)
        for offset in (-10,-6,6,10):
            signal = readings[:, n//2+offset]
            latest = 16 if abs(offset)==6 else 24
            allowed = [i for i,t in enumerate(HOLDS) if 1<=t<=latest]
            peak = max(allowed,key=lambda i: abs(signal[i]))
            rows.append(dict(material=arm, offset=offset, early=float(signal[1]),
                             peak_hold=HOLDS[peak], peak_contrast=float(signal[peak]),
                             early_consistent_zero=abs(float(signal[1]))<=radius,
                             resolved=abs(float(signal[peak]))>radius,
                             all_holds_consistent_zero=bool(max(abs(signal))<=radius)))
    return dict(joint_99_radius=radius, receivers=rows)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--acquire",type=Path)
    args = cli.parse_args()
    raw = acquire()
    result = read(raw)
    if args.acquire:
        with args.acquire.open("x") as out:
            json.dump(raw,out,separators=(",",":"),allow_nan=False)
            out.write("\n")
    elif raw!=json.loads((measure.HERE/"local_counts.json").read_text()):
        raise ArithmeticError("Local sign receipt differs")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
