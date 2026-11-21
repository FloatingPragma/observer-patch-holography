"""Resolve the four-cycle epoch of a sharp write with one material receiver.

The age-128 return records were examined before choosing this replication.
Four independent enrollment words label the mature lap classes. The blind
decoder receives just those calibration words and one terminal count word.
"""
import argparse
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
from mpmath.ctx_iv import MPIntervalContext

from .. import formation, law
from . import clock_cycle as cycle
from . import clock_factory as factory
from . import clock_inference as reader


BASE, L = 128, formation.AGE_PERIOD
CALIBRATION = tuple(BASE+h*L for h in range(4))
VISITS = (2*L+224, BASE, 3*L+8190, L+128, 4098, 2*L+8064, 3*L+302, 4*L+128)
MENUS = ("native","saturated","dephased")
N, GRID, SEED = 2**38,2**40,54110921
BANK = len(MENUS)*(len(CALIBRATION)+len(VISITS))
ARCHIVE = Path(__file__).with_name("clock_epoch_counts.npz")


def sufficient_budget():
    """Uniform single-station separation on all four mature even age classes."""
    iv = MPIntervalContext()
    iv.dps = 50
    phi = sum(iv.pi*(1-iv.cos(iv.pi*(j+1)/128))/4 for j in range(35,128))
    p = [(1+iv.cos(2*phi+h*iv.pi/2))/8 for h in range(4)]
    if not all(x>iv.mpf(".0157") for x in p) or not all(
        abs(x-y)>iv.mpf(".048") for i,x in enumerate(p) for y in p[:i]
    ):
        raise ArithmeticError("the four natal-root incidences lost separation")
    alpha = iv.mpf("0.0000005")
    r = iv.sqrt(iv.log(2*BANK/alpha)/(2*N))+iv.mpf(4)/GRID+iv.mpf("1e-13")
    k = 2_000_000_000
    log_tail = iv.log(6*BANK/alpha)
    axis_lower = iv.mpf(k)/N-iv.sqrt(log_tail/(2*N))
    q = iv.sqrt(2*log_tail/k)+4/(GRID*axis_lower)+iv.mpf("4e-11")
    chord = 2*iv.sqrt(2)*q
    angle = iv.atan2(chord,iv.sqrt(1-chord**2))
    mean = N*(iv.mpf(".0157")/2-iv.mpf(2)/GRID)
    outage = 2*BANK*iv.exp(-mean*(1-k/mean)**2/2)
    if not (4*r<iv.mpf(".048") and 2*angle+iv.mpf("1e-9")<iv.pi/2048
            and mean>k and outage<iv.mpf("1e-100")):
        raise ArithmeticError("lap or within-lap separation exceeds the booked budget")
    return {"mature_even_epochs":16128,"native_axis_outage_upper":1e-100}


def likelihoods():
    shape, root = (5,5,3,3),(2,2,1,1)
    wave = np.zeros((2,*shape,8),complex)
    wave[(slice(None),*root)] = law.EQUAL/np.sqrt(2)
    gates = factory.uniform_schedule()
    for age in range(BASE):
        wave = factory.advance(wave,age,gates)
    bank = np.empty((len(MENUS),len(CALIBRATION)+len(VISITS),6))
    for m,menu in enumerate(MENUS):
        for j,t in enumerate((*CALIBRATION,*VISITS)):
            local = cycle.transfer(wave,BASE,t-BASE,saturated=(menu=="saturated"))[(slice(None),*root)]
            # Four root effects, the remaining marked field and the refused write.
            ports = [np.linalg.norm(local[0]+p*local[1])**2/8 for p in (1,1j,-1,-1j)]
            total = sum(ports)
            if menu=="dephased":
                ports = [total/4]*4
            bank[m,j] = [*ports,.5-total,.5]
    return bank


def quantize(bank):
    numerators = np.rint(bank*GRID).astype(np.int64)
    numerators[...,4] = GRID//2-numerators[...,:4].sum(-1)
    if np.min(numerators)<0 or np.any(numerators.sum(-1)!=GRID):
        raise ArithmeticError("incomplete original-start word")
    return numerators


def decode(calibration, word):
    """Count-only lap classes and even native ages; all ties and misses remain."""
    calibration, word = np.asarray(calibration),np.asarray(word)
    if (calibration.shape!=(4,6) or word.shape!=(6,) or calibration.dtype!=np.int64
            or word.dtype!=np.int64):
        raise ValueError("four enrolled six-count calibrations and one terminal count word required")
    rows = np.vstack([calibration,word])
    if np.any(rows<0) or np.any(rows>N) or np.any(rows.sum(-1)!=N):
        raise ValueError("each original-start count word must have its booked total")
    observed = rows[:,:4].sum(-1)/N
    radius = math.sqrt(math.log(2*BANK/5e-7)/(2*N))+4/GRID+1e-13
    laps = np.flatnonzero(abs(observed[:4]-observed[4])<=2*radius).tolist()
    hand = reader.decode_elapsed(calibration[:1],word[None],rate=-math.pi/4096,
                                 cap=L-BASE-2,alpha=5e-7,family_size=BANK,
                                 cell_error=1/GRID,interval="hoeffding")
    runs = [[BASE+lo,BASE+hi] for lo,hi in hand["even_runs_step_two"]]
    return {"lap_classes":laps,"age_runs_step_two":runs,
            "candidate_count":len(laps)*hand["eligible_even_count"]}


def experiment(path=ARCHIVE, *, acquire=False, seed_replay=False):
    path = Path(path)
    exact = quantize(likelihoods())
    def draws(weights):
        return np.array([[np.random.default_rng(np.random.SeedSequence([SEED,m,j]))
                          .multinomial(N,w/GRID) for j,w in enumerate(row)]
                         for m,row in enumerate(weights)])
    if acquire:
        with path.open("xb") as handle:
            np.savez_compressed(handle,weights=exact,counts=draws(exact))
    blob = path.read_bytes()
    if not acquire and hashlib.sha256(blob).hexdigest()!=ARCHIVE.with_suffix(".sha256").read_text().split()[0]:
        raise ValueError("original epoch acquisition changed")
    with np.load(io.BytesIO(blob),allow_pickle=False) as archive:
        if sorted(archive.files)!=["counts","weights"]:
            raise ValueError("unexpected epoch acquisition keys")
        weights, counts = archive["weights"],archive["counts"]
    if (weights.shape!=exact.shape or weights.dtype!=np.int64 or np.any(weights<0)
            or np.any(weights>GRID) or np.any(weights.sum(-1)!=GRID)
            or np.any(weights[...,5]!=GRID//2) or np.max(abs(weights-exact))>2
            or np.max(abs(weights[...,:4]/GRID-likelihoods()[...,:4]))>1/GRID):
        raise ValueError("epoch source changed")
    if (counts.shape!=weights.shape or counts.dtype!=np.int64
            or np.any((weights==0)&(counts!=0)) or (seed_replay and not np.array_equal(counts,draws(weights)))):
        raise ValueError("invalid acquisition")
    records = {menu:[decode(counts[m,:4],w) for w in counts[m,4:]] for m,menu in enumerate(MENUS)}
    failed = []
    for menu,words in records.items():
        for j,(word,t) in enumerate(zip(words,VISITS,strict=True)):
            laps = list(range(4)) if menu=="saturated" else [(t//L)%4]
            ages = [[BASE,L-2]] if menu=="dephased" else [[t%L,t%L]]
            if word["lap_classes"]!=laps or word["age_runs_step_two"]!=ages:
                failed.append([menu,j])
    return {"words":records,"launches":int(counts.sum()),"coverage_failure_probability":1e-6,
            "uniform_budget":sufficient_budget(),
            "held_epochs":VISITS,"failed_gates":failed}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--acquire",action="store_true")
    parser.add_argument("--seed-replay",action="store_true")
    args = parser.parse_args()
    result = experiment(acquire=args.acquire,seed_replay=args.seed_replay)
    print(json.dumps(result,indent=2))
    raise SystemExit(int(bool(result["failed_gates"])))
