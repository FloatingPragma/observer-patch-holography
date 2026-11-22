"""Common-clock return spectroscopy of two conserved single-load sectors."""
import argparse
from pathlib import Path

import numpy as np

from .shell_source import Board, check, local_source, prepare
from .shells import check_shots, draws, hands, protocol, read_hand, save


def first_line(values, vectors):
    check(abs(values[0])<1e-11 and values[1]>1e-9, "disconnected kinetic graph")
    cluster = vectors[:,abs(values-values[1])<=1e-9]
    diagonal = (cluster**2).sum(axis=1)
    port = int(np.flatnonzero(diagonal>1e-10)[0])
    wave = cluster@cluster[port]/np.sqrt(diagonal[port])
    return wave,cluster.shape[1],port


def analyzer(phases, weights, fuel_phase, phase):
    # An inverse splitter after W†. Complementary ports go to the third outcome.
    arm = phases@weights
    ref = fuel_phase*np.exp(-1j*phase)
    missed = (abs(phases-arm[:,None])**2)@weights/2
    return np.column_stack((abs(ref+arm)**2/4,abs(ref-arm)**2/4,missed))


def times_interval(a,b):
    ends = [x*y for x in a for y in b]
    return min(ends),max(ends)


def response(counts, shots, time_count=None):
    """Read the two quadratures only after validating every joint detector bin."""
    check_shots(shots)
    counts = np.asarray(counts)
    check(counts.ndim==4 and counts.shape[0]==2 and counts.shape[2:]==(3,4)
          and counts.shape[1]>=2 and (time_count is None or counts.shape[1]==time_count)
          and counts.dtype.kind in "iu" and np.all(counts>=0) and np.all(counts<=shots),
          "incomplete or invalid joint return counts")
    totals = counts.sum(axis=-1,dtype=np.int64)
    check(np.all(totals.sum(axis=2)==shots),"incomplete return measurement")
    quadratures = (totals[:,:,0]-totals[:,:,1])/shots
    return quadratures[0]+1j*quadratures[1]


def extract(data, hand, eta, shots, time_count):
    responses = []
    arms = [response(counts,shots,time_count) for counts in data]
    check(len(arms)==2,"missing rest or moving mode")
    for z in arms:
        if min(abs(z)) <= 1-np.sin(eta)-1e-9:
            return {"status":"unresolved quadrature"},None
        continuous = np.unwrap(np.angle(z))
        responses.append(continuous[-1]-continuous[0])
    rest, line = responses
    delta = line-rest
    if hand is None or hand[-1]-hand[0]<=2*eta:
        return {"status":"unresolved hand duration"},None
    duration = hand[-1]-hand[0]
    # Work in phase advances until the final common positive time conversion.
    first = (delta-4*eta,delta+4*eta)
    second = (line+rest-4*eta,line+rest+4*eta)
    return dict(rest=rest/duration,kinetic_increment=delta/duration,
                squared_shift=(line**2-rest**2)/duration**2),times_interval(first,second)


def run(spec,base, *, observe=None):
    board = Board(spec["level"],base["matching"])
    kinetic = (board.inside+board.outside).toarray()
    values, vectors = np.linalg.eigh(kinetic)
    wave, rank, selected = first_line(values,vectors)
    targets = (np.ones(board.car.n)/np.sqrt(board.car.n),wave)
    theta = np.arange(spec["clock_steps"]+1)/128
    shots = spec["shots"]
    # Two analyzer and two clock quadratures, for all settings and times.
    marginals = 4*len(theta)*len(targets)*len(spec["loads"])*len(spec["law_members"])*len(spec["reference_terms"])
    delta = np.sqrt(2*np.log(2*marginals/.01)/shots)
    eta = float(np.arcsin(np.sqrt(2)*(delta+1e-8)))
    rng = np.random.default_rng(spec["seed"])
    if observe is None:
        def observe(key, index, probabilities):
            return draws(rng, probabilities, shots)
    results = dict(multiplicity=rank,projector_port=selected,phase_error_99=eta,
        kinetic_spectrum=(np.round(values,10)+0.).tolist(),members={})
    records = {}
    for member in spec["law_members"]:
        contrast = -8 if member=="half_rate_contrast_reversed" else 8
        for reference in spec["reference_terms"]:
            arm_data = []
            exact = []
            for amount in spec["loads"]:
                H,rate,offset = board.generator(member,amount)
                check(np.max(abs(H.toarray()-offset*np.eye(len(values))-rate*kinetic))<1e-12,
                      "two loads see different kinetic operators")
                times = theta/rate
                ev = rate*values+offset
                phases = np.exp(-1j*np.outer(times,ev))
                clock = hands(local_source(rate,contrast,amount)[0],times)
                fuel = np.exp(-1j*contrast*reference*amount**2*times)
                modes = []
                expected = (np.array([0,rate*values[1]])+offset-contrast*reference*amount**2)/rate
                check(np.max(abs(expected))*np.diff(theta).max()+2*eta<np.pi,"frequency can alias")
                exact.append(expected)
                for name,target,frequency in zip(spec["modes"],targets,expected,strict=True):
                    psi,_ = prepare(board,target,amount,root=0)
                    weights = abs(vectors.conj().T@psi)**2
                    check(abs(sum(weights)-1)<1e-12,"off-spectrum arm")
                    # Bound preparation/eigensolver error throughout the observed interval.
                    check(np.max(abs(phases@weights*fuel.conj()-np.exp(-1j*frequency*theta)))<1e-10,
                          "eigenline visibility or phase lost")
                    key = f"{member}_{reference}_{amount}_{name}"
                    trials = []
                    for analyzer_index, phi in enumerate((0., np.pi/2)):
                        probabilities = analyzer(phases,weights,fuel,phi)
                        trials.append([observe(key, (analyzer_index, n), p[:, None]*hand)
                                       for n, (p, hand) in enumerate(zip(probabilities, clock, strict=True))])
                    rows = np.asarray(trials)
                    records[key] = rows
                    modes.append(rows)
                arm_data.append(modes)
            hand = read_hand(sum(c.sum(axis=(0,2),dtype=np.int64) for modes in arm_data for c in modes),
                             allow_unresolved=True)
            duration = None if hand is None else hand[-1]-hand[0]
            observed, bounds = zip(*(extract(modes,hand,eta,shots,len(theta)) for modes in arm_data))
            interval = None
            if all(bound is not None for bound in bounds):
                low, high = bounds[1][0]-bounds[0][1],bounds[1][1]-bounds[0][0]
                interval = list(times_interval((low,high),((duration+2*eta)**-2,(duration-2*eta)**-2)))
            truth = np.array([w[1]**2-w[0]**2 for w in exact])
            results["members"][f"{member}_{reference}"] = dict(loads=list(observed),
                difference=float(truth[1]-truth[0]),measured_difference_interval=interval,
                common_squared_shift="excluded" if interval and (interval[0]>0 or interval[1]<0) else "unresolved",
                hand_duration=duration)
    return records,results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    counts,result = run(protocol("shell_spectroscopy_protocol.json"),protocol("shell_protocol.json"))
    save("shell_spectroscopy",counts,result,args.output)


if __name__=="__main__":
    main()
