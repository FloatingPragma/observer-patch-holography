"""Freely propagating loads and co-located collar hands on the three registered screens."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.experiments.seam_program.gate import archive, receipt
from .shell_source import Board, check, local_source

HERE = Path(__file__).resolve().parent
SOURCE_FILES = {"shell_source.py","shells.py","shell_spectroscopy.py","shell_protocol.json",
                "shell_spectroscopy_protocol.json","../../executable_model/carrier.py",
                "../../candidate_laws/contrast_memory_spherical_v1/law.py",
                "../seam_program/gate.py",
                *("../../executable_model/carriers/"+name+".npz" for name in
                  ("CELL","FEDERATION_LEVEL0","FEDERATION_LEVEL1","FEDERATION_LEVEL2"))}


def protocol(name):
    pins = json.loads((HERE/"shell_sources.json").read_text(encoding="utf8"))
    check(set(pins)==SOURCE_FILES and name in ("shell_protocol.json","shell_spectroscopy_protocol.json"),
          "incomplete source manifest or unknown experiment")
    for file, digest in pins.items():
        raw = (HERE/file).read_bytes()
        actual = sha256(raw if file.endswith(".npz") else raw.replace(b"\r\n",b"\n")).hexdigest()
        check(actual==digest,f"changed experimental source: {file}")
    return json.loads((HERE/name).read_text(encoding="utf8"))


def hands(sign_term, times):
    """Retain all four joint sharp sign outcomes. Columns are ++, +-, -+, --."""
    even, odd = np.array([1,1])/np.sqrt(2), np.array([1,-1])/np.sqrt(2)
    start = np.column_stack(((even+odd)/np.sqrt(2),(even+1j*odd)/np.sqrt(2)))
    eigen, vectors = np.linalg.eigh(sign_term)
    evolved = np.einsum("ij,tj,jk->tik",vectors,np.exp(-1j*np.outer(times,eigen-eigen[0])),vectors.T@start)
    p = abs(evolved)**2
    return (p[:,:,0,None]*p[:,None,:,1]).reshape(-1,4)


def read_hand(counts):
    check(counts.shape[-1]==4 and np.isfinite(counts).all() and counts.min()>=0
          and np.all(counts.sum(axis=-1)>0),"missing clock outcomes")
    p = counts/counts.sum(axis=-1,keepdims=True)
    x = p[...,0]+p[...,1]-p[...,2]-p[...,3]
    y = p[...,0]-p[...,1]+p[...,2]-p[...,3]
    check(np.min(np.hypot(x,y))>.5,"unresolved phase hand")
    return np.unwrap(np.arctan2(y,x),axis=-1)


def draws(rng, joint, shots):
    check(np.isfinite(joint).all() and joint.min()>=0 and
          abs(joint.sum()-1)<1e-10,"incomplete receiver distribution")
    # Quantize the cumulative probabilities before sampling for reproducibility across BLAS builds.
    cdf = np.rint(np.cumsum((joint/joint.sum()).ravel())*2**36).astype(np.int64)
    p = np.diff(cdf,prepend=0)/2**36
    return rng.multinomial(shots,p).astype(np.int32).reshape(joint.shape)


def run(spec):
    from .shell_source import prepare
    rng = np.random.default_rng(spec["seed"])
    shots = spec["shots"]
    panels = len(spec["launches"])*len(spec["laws"])*sum(4*4**l+1 for l in spec["levels"])
    # Union of the two-sided subset tests on each twenty-region distribution.
    tv_budget = np.sqrt((20*np.log(2)+np.log(panels/.01))/(2*shots))
    result = dict(angular_pair_99pct_budget=2*tv_budget,boards={})
    records = {}
    for level in spec["levels"]:
        board = Board(level,spec["matching"])
        ns = np.arange(4*4**level+1)
        native,_,shift = board.generator("native")
        initial, resources = [], []
        for name in spec["launches"]:
            state, cost = prepare(board,board.launch(name))
            initial.append(state)
            resources.append(cost)
        initial = np.column_stack(initial)
        energy = np.real(np.sum(initial.conj()*(native@initial-shift*initial),axis=0))
        # Every smooth angular launch is Lipschitz on the face graph and has unit plateau on an entire macroface.
        check(np.all(energy[:2]<=5*board.h**2/16+1e-13),"smooth launch exceeds angular gradient budget")
        sets = {}
        for member in spec["laws"]:
            H,rate,offset = board.generator(member)
            contrast = -8 if member=="half_rate_contrast_reversed" else 8
            hand = hands(local_source(rate,contrast)[0],ns/rate)
            # Remove the common phase for position probabilities only.
            step = -1j*(H-offset*sparse.eye(board.car.n)).tocsr()
            path = expm_multiply(step,initial,start=0,stop=ns[-1]/rate,num=len(ns))
            exact = np.empty((len(initial.T),len(ns),20,2))
            for i,name in enumerate(spec["launches"]):
                counts = []
                for n,state in enumerate(path[:,:,i]):
                    region = np.bincount(board.bins,weights=abs(state)**2,minlength=40).reshape(20,2)
                    exact[i,n] = region
                    counts.append(draws(rng,region[:,:,None]*hand[n],shots))
                    if member=="native":
                        uniform = state.reshape(-1,12).mean(axis=1,keepdims=True)
                        eta = np.linalg.norm(state.reshape(-1,12)-uniform)**2
                        check(eta <= energy[i]/board.gap+1e-12,"internal spectral gap bound failed")
                        imbalance = abs(region[:,0]-region[:,1]).sum()
                        check(imbalance <= 2*np.sqrt(max(eta*(1-eta),0))+eta+1e-12,
                              "depth contrast violates projection bound")
                        check(np.linalg.norm(state-initial[:,i])**2 <= 2*n*energy[i]+1e-11,
                              "ballistic-time displacement bound failed")
                records[f"{level}_{name}_{member}"] = np.asarray(counts)
            sets[member] = exact
        check(np.max(abs(sets["native"]-sets["half_rate_contrast_reversed"]))<1e-11,
              "same-clock equivalence failed")
        out = {}
        for i,name in enumerate(spec["launches"]):
            def angular(member):
                return records[f"{level}_{name}_{member}"].sum(axis=(2,3))/shots
            base = angular("native")
            peak = {m:float(np.max(np.abs(base-angular(m)).sum(axis=1)/2)) for m in spec["laws"][1:]}
            check(peak["half_rate_contrast_reversed"]<2*tv_budget,"registered equal-law count screen failed")
            err = max(np.max(abs(read_hand(records[f"{level}_{name}_{m}"].sum(axis=(1,2)))-ns))
                      for m in spec["laws"])
            depth = abs(sets["native"][i,:,:,0]-sets["native"][i,:,:,1]).sum(axis=1)/2
            moved = abs(sets["native"][i].sum(axis=2)-sets["native"][i,0].sum(axis=1)).sum(axis=1)/2
            out[name] = dict(kinetic_energy=float(energy[i]),source=resources[i],peak_control_angular_tv=peak,
                 control_separated={m:bool(v>2*tv_budget) for m,v in peak.items()},max_hand_error=float(err),
                 ballistic_region_tv=[float(moved[t*2**level]) for t in (1,2,4)],
                 slow_time_region_tv=[float(moved[t*4**level]) for t in (1,2,4)],
                 peak_depth_tv=float(max(depth)))
        result["boards"][str(level)] = dict(cells=board.car.cells,internal_gap=board.gap,
            signed_seams=6*board.car.cells,active_contacts=len(board.active),launches=out)
    return records,result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    data,summary = run(protocol("shell_protocol.json"))
    save("shell",data,summary,args.output)


def save(name,data,summary,output):
    folder = output or HERE
    if output:
        folder.mkdir(exist_ok=True,parents=True)
    archive(folder/(name+"_counts.npz"),data,output is not None)
    receipt(folder/(name+"_result.json"),summary,output is not None)
    print(f"{name}: {sum(int(a.sum()) for a in data.values())} retained simulated trials")


if __name__=="__main__":
    main()
