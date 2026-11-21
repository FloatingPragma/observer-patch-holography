"""Local spectral lower bounds after an agreed positive encounter on Interval Echo."""
import argparse
from collections import deque
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm, toeplitz
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import presence
from model.candidate_laws.contrast_memory_spherical_v1 import law as seam_law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation, seam_digest

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "return_counts.npz"
SUMMARY = HERE / "return_lines.json"


def check(ok, message):
    if not ok:
        raise ValueError(message)


def distances(car, start):
    depth = [-1]*car.n
    depth[start] = 0
    queue = deque([start])
    while queue:
        vertex = queue.popleft()
        for other in car.neighbours[vertex]:
            if depth[other] < 0:
                depth[other] = depth[vertex]+1
                queue.append(other)
    check(min(depth) == 0, "disconnected carrier")
    return depth


def line_bound(car, seam):
    """Every exchange moves exactly one token along a seam. Keep a distant witness."""
    a, b = (distances(car, i) for i in seam)
    distance, target = max((min(a[u]+b[v],a[v]+b[u]), (u,v)) for u,v in car.seams)
    return dict(seam_sha256=seam_digest(car), target=list(target), moves_at_least=distance,
                visible_energies_at_least=distance+1, positive_return_gaps_at_least=distance)


def counter(last):
    """Read the travelling double on a separate equipped native seam after each pi gate."""
    check(type(last) is int and 0<=last<32, "booking exceeds five-bit local counter")
    law = candidate.native_law
    h0,h1 = law.local_term_parts()
    swap = expm(-1j*np.pi*(h0+candidate.KAPPA*h1))
    state = np.eye(law.DIM, dtype=complex)[:,law.SRC.index((2,0))]
    echo = law.echo()
    agreed = np.arange(law.DIM)%len(law.SRC) == law.SRC.index((1,1))
    sources = [law.SRC.index((0,2)),law.SRC.index((2,0))]
    inputs = np.eye(law.DIM)[:,[law.SRC.index((0,0)),law.SRC.index((2,0))]]
    written = echo @ inputs
    effect = written[agreed].conj().T @ written[agreed]
    check(np.max(abs(effect-np.diag([0,1])))<2e-13 and np.linalg.norm(echo@written-inputs)<2e-13,
          "native detector does not read the loaded endpoint")
    words = [0]
    for index in range(1,last+1):
        state = swap @ state
        weights = abs(state[sources])**2
        check(abs(sum(weights)-1)<2e-12, "reference left the two endpoint configurations")
        # Each endpoint probes its load against blank. The opposite endpoint is
        # orthogonal on the two source rays, so its marginal is diagonal.
        powers = np.array([weights @ effect.diagonal(), weights[::-1] @ effect.diagonal()]).real
        flags = np.rint(powers).astype(int)
        check(np.max(abs(flags-powers))<2e-12,"ambiguous reference occupation")
        check(np.array_equal(flags, [1,0] if index%2==0 else [0,1]), "missed reference return")
        words.append((words[-1]+int(flags[index%2]))%32)
    return words


def acquire_laws(cfg):
    car = federation(cfg["acquisition_level"])
    sector = PairSector(car)
    home = sector.index(*cfg["source_seam"])
    swap_to = sector.index(0,5)
    check((1,5) in car.seams and candidate.KAPPA==8 and candidate.native_law is seam_law,
          "registered contact changed")
    starts = [sector.index(i,j) for i,j in combinations(range(12),2) if (i,j) not in car.seams]
    check(len(starts)==36, "wrong separated preparation")
    n = cfg["intervals"]
    clocks = counter(n+1)
    check(clocks==list(range(n+2)), "missing clock readings")
    # Two public port readers and their courier implement the local predicate.
    labels = presence.compile_records(np.array([[False],[True]]), ["nominal"], True)
    check(labels==[(0,(1,)),(1,(29,))], "first/last vote does not implement the sharp predicate")
    check(np.array_equal(sector.exchange().data[sector.exchange().data<0],
                         np.full(np.count_nonzero(sector.exchange().data<0),-.5)),
          "pair transitions changed sign or magnitude")
    e = np.zeros(sector.P)
    e[home] = 1
    for menu,kappa in (("native",candidate.KAPPA),("exchange",0)):
        h = sector.hamiltonian(kappa)
        check((h-h.T).nnz == 0, "backward columns require a real symmetric free law")
        # Transposition, without conjugation, supplies every launch column at the local seam.
        evolved = expm_multiply(-1j*np.pi*h,e,start=0,stop=n+1,num=n+2)
        check(np.max(abs(np.sum(abs(evolved)**2,axis=1)-1))<3e-10,"lost sector probability")
        first = evolved[clocks[1],starts]
        admit = float(np.mean(abs(first)**2))
        check(0 < admit < 1, "no conditional source")
        stay = evolved[clocks[:-1],home]
        bypass = evolved[np.ix_(clocks[1:],starts)]
        yes = admit*abs(stay)**2
        quiet = np.mean(abs(bypass-stay[:,None]*first)**2,axis=1)
        yes[0], quiet[0] = admit, 0.  # Identical sharp questions with no free service.
        yield menu, bins(admit,yes,quiet)
        if menu=="native":
            yield "park",bins(admit,np.full(n+1,admit),np.zeros(n+1))
            alt = np.zeros(sector.P)
            alt[swap_to] = 1
            alternate = expm_multiply(-1j*np.pi*h,alt)
            alt_weight = float(np.mean(abs(alternate[starts])**2))
            # At pi the isolated seam swaps occupations 1 and 5 up to a common sharp phase.
            law = candidate.native_law
            h0,h1 = law.local_term_parts()
            one_swap = expm(-1j*np.pi*(h0+candidate.KAPPA*h1))
            for x,y in ((0,1),(1,0)):
                check(abs(abs(one_swap[law.SRC.index((y,x)),law.SRC.index((x,y))])-1)<1e-12,
                      "isolated exchange control is not a full swap")
            odd = np.arange(n+1)%2==1
            yield "isolated_swap",bins(admit,np.where(odd,0,admit),np.where(odd,alt_weight,0))


def bins(admit, yes, quiet):
    p = np.column_stack((yes,admit-yes,quiet,1-admit-quiet))
    check(np.all(np.isfinite(p)) and p.min()>=-1e-12 and np.all(abs(p.sum(axis=1)-1)<1e-13),
          "invalid joint law")
    return np.maximum(p,0)


def outside_weight(eigen, error, levels):
    """Minimum source weight outside any chosen set of energy levels."""
    rank = 1 + levels*(levels-1)
    tail = eigen[:max(0,len(eigen)-rank)]
    lost = sum(np.maximum(0,tail-error))/len(eigen)
    return 1-math.sqrt(max(0,1-lost))


def read(data, shots, logtails):
    """Infer from pairs of public words only. Each time uses an independent cohort."""
    check(isinstance(data,np.ndarray) and data.dtype==np.int64 and data.ndim==2
          and data.shape[1]==4 and data.shape[0]>=4 and type(shots) is int
          and 0<shots<=2**60 and np.isfinite(logtails) and logtails>0
          and np.all((data>=0)&(data<=shots)) and np.all(data.sum(axis=1)==shots),
          "incomplete or invalid joint records")
    heralds = data[:,:2].sum(axis=1)
    check(np.all(heralds>0), "unresolved herald")
    check(data[0,1]==data[0,2]==0,"failure of immediate repeatability")
    measured = data[:,0]/heralds
    radius = np.sqrt(logtails/(2*heralds))
    # P(0)=1 by identical sharp reads. Every other Toeplitz diagonal uses a measured frequency.
    radius[0] = 0
    lower,upper = np.maximum(0,measured-radius),np.minimum(1,measured+radius)
    rank_error = 2*sum(radius[1:])
    eigen = np.linalg.eigvalsh(toeplitz(measured))
    numerical = 1e-11
    rank = max(1,int(np.count_nonzero(eigen > rank_error + numerical)))
    gaps = rank//2
    energies = math.ceil((1+math.sqrt(4*rank-3))/2)
    dmin,dmax = 1-upper,1-lower
    single_gap_excess = dmin[1]*dmin[3] - max(abs(dmin[2]-dmax[1]),abs(dmax[2]-dmin[1]))**2
    return dict(positive_gaps_at_least=gaps,energies_at_least=energies,rank_at_least=rank,
                outside_weight_lower={str(d):math.floor(1e7*outside_weight(eigen,rank_error+numerical,d))/1e7
                                      for d in (1,2,3,4)},
                rank_error=round(rank_error,8),smallest_eigenvalue=round(float(eigen[0]),8),
                first_three_excess_lower=round(float(single_gap_excess),8),
                conditional_return=measured.round(8).tolist(),heralds=heralds.tolist())


def run(record=False):
    protocol = (HERE/"return_protocol.json").read_bytes().replace(b"\r\n",b"\n")
    check(hashlib.sha256(protocol).hexdigest()=="0e29b10e2c596377f66fa96eeba520f64c0924150cf760a209766970c01fbd85", "registered return protocol changed")
    cfg = json.loads(protocol)
    shots = cfg["shots_per_appointment"]
    forecasts = dict(acquire_laws(cfg))
    check(set(forecasts)==set(cfg["menus"]), "registered menu changed")
    names = cfg["menus"]
    if record:
        rng = np.random.default_rng(cfg["seed"])
        fresh = {}
        for name in names:
            # Sequential binomials avoid platform-dependent multinomial trial limits.
            forecasts_for_menu = forecasts[name]
            a = forecasts_for_menu[:,:2].sum(axis=1)
            enrolled = rng.binomial(shots,a)
            yes = rng.binomial(enrolled,np.clip(forecasts_for_menu[:,0]/a,0,1))
            quiet = rng.binomial(shots-enrolled,np.clip(forecasts_for_menu[:,2]/(1-a),0,1))
            fresh[name] = np.column_stack((yes,enrolled-yes,quiet,shots-enrolled-quiet)).astype(np.int64)
        with ARCHIVE.open("xb") as output:
            np.savez_compressed(output,**fresh)
    budget = np.log(2*5*len(names)*(cfg["intervals"]+1)/cfg["alpha"])
    results = {}
    with np.load(ARCHIVE,allow_pickle=False) as data:
        check(sorted(data.files)==sorted(names), "lost acquisition arm")
        for name in names:
            counts,p = data[name],forecasts[name]
            check(counts.shape==p.shape, "missing clock appointment")
            results[name] = read(counts,shots,budget)
            # Simultaneous binomial Bernstein bounds on all unconditional frequencies.
            deviation = np.sqrt(2*p*(1-p)*budget/shots)+2*budget/(3*shots)+2e-12
            check(np.all(abs(counts/shots-p)<=deviation) and np.all(counts[p==0]==0),
                  f"joint public records disagree on {name}")
    result = dict(distance_witnesses={str(j):line_bound(federation(j),cfg["source_seam"])
                                      for j in cfg["carrier_levels"]},local_records=results,
                  starts=shots*len(names)*(cfg["intervals"]+1),
                  sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest())
    if record:
        with SUMMARY.open("x",encoding="utf8",newline="\n") as stream:
            stream.write(json.dumps(result,indent=1,sort_keys=True)+"\n")
    check(json.loads(SUMMARY.read_text(encoding="utf8"))==result,"changed spectral certificate")
    return result


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",action="store_true")
    result = run(parser.parse_args().record)
    print("Minimum measured ranks", {key:v["rank_at_least"] for key,v in result["local_records"].items()})
