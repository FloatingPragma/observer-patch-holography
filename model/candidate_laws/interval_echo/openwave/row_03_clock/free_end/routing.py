"""Coherent two-path energy comparison on the same herald, with an explicit apparatus."""

import argparse
import json

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares

from ...row_08_annihilation import experiment
from . import measure as clock

N, TIME, SEED = 512, np.arange(512)/16, 310395
# Each bank specifies source, reference energy label, memory flight, and switch dephasing.
PANELS = dict(ground=(1,0,False,False), sham=(0,0,False,False),
              excited_reference=(1,2,False,False), matched_flight=(1,0,True,False),
              dephased=(1,0,False,True))


def phase_bank(levels, weights, ref, matched, dephase):
    delta = levels - levels[ref]
    if matched:
        delta = np.zeros_like(delta)
    phasors = np.exp(1j*np.outer(TIME,delta))
    labels = np.array([False,True,True,True,False,False,False,False])
    p = np.zeros((2,N,2,2,2))
    for tag in (0,1):
        occupied = weights*(labels == tag)
        mass = occupied.sum()
        coherent = phasors @ occupied
        if dephase:
            coherent = np.zeros_like(coherent)
        for axis, mean in enumerate((coherent.real,coherent.imag)):
            for i,b in enumerate(clock.SIGNS):
                # Both arm unitaries commute with the same target number flag after unrouting.
                p[axis,:,tag,i,tag] = (mass+b*mean)/2
    return p


def acquire():
    rng = np.random.default_rng(SEED)
    pin, waves = clock.collision()
    archive = {}
    for name,(source,ref,matched,dephased) in PANELS.items():
        banks, energy = [], []
        for phase in (0,1):
            e,v = eigh(experiment.wing(8,3,pin if phase == 0 else 2-pin))
            state = v.T @ waves[phase,:,:,source]
            weights = np.sum(abs(state)**2,axis=1)
            energy.extend(weights)
            panel = phase_bank(e, weights, ref, matched, dephased)
            shifted = phase_bank(e+7, weights, ref, matched, dephased)
            clock.check(np.max(abs(panel-shifted))<1e-12, "common energy shift registered")
            banks.append(panel)
        miss = 1-sum(energy)
        paired = np.stack(banks,axis=2).reshape(2,N,-1)
        archive[name] = np.stack([clock.sample(np.c_[p,np.full(N,miss)],rng) for p in paired])
        if name in ("ground","sham"):
            archive[name+"_energy"] = clock.sample(np.array([[*energy,miss]]),rng)
    return archive


def extract(bank, tag):
    cells = bank[:,:,:-1].reshape(2,N,2,2,2,2).sum(axis=2)
    wrong = cells[:,:,0,:,1].sum()+cells[:,:,1,:,0].sum()
    if tag is None:
        cells = cells.sum(axis=2)
    else:
        cells = cells[:,:,tag]
    signs = cells.sum(axis=-1)
    eligible = signs.sum(axis=-1)
    clock.check(np.min(eligible)>0, "no admitted sign")
    mean = (signs[:,:,1]-signs[:,:,0])/eligible
    return mean[0]+1j*mean[1], int(wrong), int(eligible.sum())


def fit(series):
    """Blind line pursuit on complex pointer records from first half only."""
    t = TIME[:N//2]
    train = series[:N//2]
    residual = train.copy()
    seeds = []
    for _ in range(3):
        f = 2*np.pi*np.fft.fftfreq(8192,1/16)
        power = abs(np.fft.fft(residual*np.hanning(len(t)),len(f)))
        valid = abs(f) <= 45
        for nu in seeds:
            valid &= abs(f-nu)>.32
        seeds.append(float(f[np.argmax(np.where(valid,power,-1))]))
        terms = np.exp(1j*np.outer(t,seeds))
        residual = train - terms @ np.linalg.lstsq(terms,train,rcond=None)[0]
    seeds = np.sort(seeds)
    def predict(x,t):
        return np.exp(1j*np.outer(t,x[:3])) @ x[3:]
    def errors(x):
        d = predict(x,t)-train
        return np.r_[d.real,d.imag]
    x = least_squares(errors,np.r_[seeds,[.3,.4,.3]],
                      bounds=(np.r_[seeds-.3,[0]*3],np.r_[seeds+.3,[1]*3])).x
    delta = predict(x,TIME[N//2:])-series[N//2:]
    return dict(lines=x[:3].tolist(),weights=x[3:].tolist(),held_complex_rmse=float(np.mean(abs(delta)**2)**.5))


def decode(archive):
    clock.check(set(archive) == set(PANELS)|{"ground_energy","sham_energy"}, "missing routing bank")
    results = {}
    for name,bank in archive.items():
        expected = (1,17) if name.endswith("energy") else (2,N,17)
        clock.check(bank.shape == expected and np.issubdtype(bank.dtype,np.integer)
                    and np.min(bank)>=0 and np.all(bank.sum(axis=-1)==clock.SHOTS), "incomplete counts")
        if name.endswith("energy"):
            counts = bank[0,:-1].reshape(2,8).sum(axis=0)
            active = counts[1:4]
            results[name] = dict(ordered=int(counts.sum()),one=int(active.sum()),
                                 one_counts=active.tolist(),miss=int(bank[0,-1]))
            continue
        z, wrong, received = extract(bank, None if name=="sham" else 1)
        result = dict(herald_reversals=wrong,received=received,starts=2*N*clock.SHOTS)
        if name in ("ground","excited_reference"):
            result.update(fit(z))
        else:
            expected = 0 if name=="dephased" else 1
            result["max_distance_from_expected"]=float(max(abs(z-expected)))
            result["mean_real"]=float(z.real.mean())
            result["mean_imag"]=float(z.imag.mean())
        results[name]=result
    return results


def validate(r):
    reference_path = clock.HERE.parents[1]/"row_08_annihilation"/"daughter_spectrum"/"results.json"
    e = np.array(json.loads(reference_path.read_text())["target_frequencies"])
    counts = np.array(r["ground_energy"]["one_counts"])
    actual = counts/counts.sum()
    for key,lines in [("ground",e),("excited_reference",e-e[1])]:
        panel = r[key]
        clock.check(max(abs(lines-panel["lines"]))<.01 and max(abs(actual-panel["weights"]))<.03
                    and panel["held_complex_rmse"]<.025, "relational energy transfer unresolved")
    for k in PANELS:
        clock.check(r[k]["herald_reversals"]==0, "herald lost")
    clock.check(r["matched_flight"]["max_distance_from_expected"]<.06 and
                r["sham"]["max_distance_from_expected"]<.015 and
                r["dephased"]["max_distance_from_expected"]<.075, "reference control failed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--reacquire", action="store_true")
    args = parser.parse_args()
    path, summary=clock.HERE/"routing_counts.npz",clock.HERE/"routing_results.json"
    if args.write:
        banks = acquire()
        r = decode(banks)
        validate(r)
        np.savez_compressed(path,**banks)
        summary.write_text(json.dumps(r,indent=2)+"\n")
    else:
        with np.load(path,allow_pickle=False) as stored:
            r = decode(stored)
            validate(r)
            expected = json.loads(summary.read_text())
            clock.check(r.keys()==expected.keys(),"summary panels differ")
            for key,row in r.items():
                clock.check(row.keys()==expected[key].keys() and all(
                    np.allclose(v,expected[key][f],atol=3e-7,rtol=3e-7) for f,v in row.items()),
                    "summary drift")
            if args.reacquire:
                fresh = acquire()
                clock.check(all(np.array_equal(stored[k],v) for k,v in fresh.items()),"sampling changed")
    print("Relational energy records verified.")


if __name__ == "__main__":
    main()
