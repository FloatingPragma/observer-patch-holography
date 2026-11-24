"""Independent material clock histograms and counting precision of the signed area test."""
import argparse
from itertools import combinations
import json

import numpy as np

from . import measure
from ..collar_pairing.measure import compare

SHOT_LEVELS = (2**30, 2**60)
EDGES = list(combinations(range(4), 2))


def clock_bank():
    banks = [{}, {}]
    times = []
    def record(first, second, legs, free, hand):
        k = measure.radar.TICKS
        appointment = np.arange(1, k+1)
        arrivals = (hand[np.add.outer(appointment, appointment)]/2, hand[appointment])
        joint = measure.radar.round_trip_law(legs)
        if not times:
            times.extend([a.ravel() for a in arrivals])
        for i, j in EDGES:
            returns = joint[i,j,:-1,:-1].ravel()
            pulse = (free[i,j]+free[j,i])/(2*k)
            for bank, hits in zip(banks,(returns,pulse),strict=True):
                bank.setdefault((first,second),[]).append(np.r_[hits,1-hits.sum()])
    measure.read_word_laws(0, True, record_independent=record)
    return [np.stack([bank[k] for k in sorted(bank)]) for bank in banks], times


def shell_area(distances):
    d = np.zeros((4,4))
    for (i,j), length in zip(EDGES,distances,strict=True):
        d[i,j] = d[j,i] = length
    zero = np.zeros_like(d)
    faces = list(combinations(range(4),3))
    area = sum(measure.simplex(d,zero,f)[0] for f in faces)
    gradient = []
    for i,j in EDGES:
        unit = np.zeros_like(d)
        unit[i,j] = unit[j,i] = 1
        gradient.append(sum(a*da for a,da,_ in (measure.simplex(d,unit,f) for f in faces)))
    return area,np.array(gradient)


def read(distribution, hand):
    hits = distribution[:,:-1]
    success = hits.sum(axis=1)
    means = hits@hand/success
    # Misses have zero influence on the ratio of summed times to number of hits.
    variance = (hits*(hand-means[:,None])**2).sum(axis=1)/success**2
    area, gradient = shell_area(means)
    return area, (gradient**2)@variance


def clock_counts(rng, shots, probabilities):
    # Some return-time bin probabilities are below machine epsilon. A naively
    # rounded 1-p loses genuine hits on these very large simulated cohorts.
    batch = min(shots, 2**50)
    left = np.full(shots//batch, batch, dtype=np.int64)
    tail = np.cumsum(probabilities[::-1])[::-1]
    counts = []
    for prob, remaining in zip(probabilities[:-1], tail[:-1], strict=True):
        p = prob/remaining
        if p == 0:
            counts.append(0)
            continue
        if p < 1e-8 and left.max()*p < 24:
            uniform = rng.random(len(left))
            mass = np.exp(left*np.log1p(-p))
            cdf = mass.copy()
            draws = np.zeros_like(left)
            select = uniform > cdf
            while np.any(select):
                k = draws[select]
                mass[select] *= (left[select]-k)*p/((k+1)*(1-p))
                cdf[select] += mass[select]
                draws[select] += 1
                select = uniform > cdf
        else:
            draws = rng.binomial(left,p)
        counts.append(int(draws.sum()))
        left -= draws
    counts.append(int(left.sum()))
    return counts


def run():
    witness = json.loads((measure.HERE/'signed_area.json').read_text())['boards'][0]
    stored = json.loads((measure.HERE/'results.json').read_text())['screens'][0]
    reference = np.array(stored['word_moments']).reshape(9,len(measure.PAIRS),2,2)
    banks,times = clock_bank()
    rng = np.random.default_rng(291147)
    result = {}
    for menu,(bank,hand) in enumerate(zip(banks,times,strict=True)):
        moments = np.stack((bank[:,:,:-1].sum(axis=-1),bank[:,:,:-1]@hand),axis=-1)
        indices = [measure.PAIRS.index(e) for e in EDGES]
        np.testing.assert_allclose(moments,reference[:,indices,menu,:],atol=5e-13,rtol=0)
        rows = {}
        for label,test in witness['witnesses'].items():
            populations = [np.einsum('i,iej->ej', w, bank) for w in test['finite']['probabilities']]
            calculated = [read(probs, hand) for probs in populations]
            change = calculated[1][0]-calculated[0][0]
            # The independent dense propagation also has to resolve the signed contrast.
            check = test['finite']['area_secants'][menu]/16
            np.testing.assert_allclose(change, check, rtol=5e-5, atol=3e-13)
            per_release_variance = sum(v for _,v in calculated)
            acquisitions = []
            for shots in SHOT_LEVELS:
                records, observed = [], []
                for probs in populations:
                    counts = np.asarray([clock_counts(rng, shots, p) for p in probs])
                    hits = counts[:,:-1].sum(axis=1)
                    expected = shots*probs[:,:-1].sum(axis=1)
                    if np.any(abs(hits-expected) > 9*np.sqrt(expected)+12):
                        raise ValueError("Clock success counts disagree with the independent flight law")
                    if np.any(counts.sum(axis=1) != shots):
                        raise ValueError("Clock bank lost releases")
                    records.append(counts.tolist())
                    if np.any(counts[:,:-1].sum(axis=1) == 0):
                        observed.append(None)
                    else:
                        try:
                            observed.append(float(read(counts/shots, hand)[0]))
                        except ArithmeticError:
                            observed.append(None)
                acquisitions.append(dict(shots_per_edge_arm=shots, counts=records, observed_areas=observed,
                                         observed_area_change=(observed[1]-observed[0]
                                                               if None not in observed else None),
                                         asymptotic_counting_se=np.sqrt(per_release_variance/shots)))
            rows[label] = dict(population_area_change=change, acquisitions=acquisitions,
                               projected_shots_per_edge_arm_for_five_se=25*per_release_variance/change**2)
        result[measure.NAMES[menu]]=rows
    return dict(cells=20, edges=EDGES, clock_times=[t.tolist() for t in times],
                menus=result)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    result=run()
    dest=measure.HERE/'signed_counts.json'
    # Read back through JSON to normalize tuple-valued edge labels before comparison.
    result=json.loads(json.dumps(result,allow_nan=False))
    if parser.parse_args().write:
        dest.write_text(json.dumps(result,indent=2)+'\n')
    else:
        compare(result,json.loads(dest.read_text()))
    for menu,rows in result['menus'].items():
        for label,row in rows.items():
            print(menu,label,row['population_area_change'],row['projected_shots_per_edge_arm_for_five_se'],
                  [(a['observed_area_change'], a['asymptotic_counting_se']) for a in row['acquisitions']])
