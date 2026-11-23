"""Compare exterior log multiplicity with a counted tag transmission on subdivided patches."""

import argparse
import hashlib
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.special import gammaln, logsumexp

from model.candidate_laws.interval_echo import candidate
from model.executable_model.refinement import load_refinement
from model.executable_model.carrier import CELL_SEAMS
from ..shape_stationarity.geometry import selected_shapes


HERE = Path(__file__).resolve().parent
SEED = 761932
AREA_SHOTS = 2**34
INFO_SHOTS = 2**24
FILLINGS = ('two', 'quarter', 'half', 'three_quarters')


def tag_probability():
    law = candidate.native_law
    h = law.seam_term(candidate.KAPPA)
    start, stop = [law.SRC.index(q) for q in ((1, 0), (0, 1))]
    pulse = expm(-.5j * np.pi * h)[:, start]
    p = abs(pulse[stop])**2
    np.testing.assert_allclose([p, abs(pulse[start])**2], [.5, .5], atol=5e-14)
    return float(p)


def class_law(v, m, n):
    k = np.arange(max(0, n+m-v), min(m, n)+1, dtype=int)
    low = k[:-1]
    steps = np.log((m-low).astype(float)*(n-low)/((low+1.)*(v-m-n+low+1)))
    unnormalized = np.r_[0., np.cumsum(steps)]
    probs = np.exp(unnormalized-logsumexp(unnormalized))
    probs /= probs.sum()
    log_choices = gammaln(m+1.)-gammaln(k+1.)-gammaln(m-k+1.)
    return k, probs, log_choices


def edge_examples(k, probability, cut, rng):
    return [sorted(rng.choice(cut, size=int(count), replace=False).tolist())
            for count in rng.choice(k, p=probability, size=2)]


def weighted_score(counts, scores):
    probs = counts / counts.sum()
    mean = float(probs @ scores)
    return mean, math.sqrt(float(probs @ (scores-mean)**2)/counts.sum())


def tangent(v, m, n, k, prob, choices, strength, rng):
    center = m//2
    width = max(1, math.isqrt(m)//2)
    idx = np.searchsorted(k, [center-width, center, center+width])
    step = min(prob[idx])/4
    direction = np.zeros(len(k))
    direction[idx] = [.5, -1., .5]
    arms = []
    for sign in (-1, 1):
        law = prob + sign*step*direction
        counts_bulk = rng.multinomial(INFO_SHOTS, law)
        counts_edge = rng.multinomial(INFO_SHOTS, law)
        # Constant log binom(V,N) drops out of all normalized bulk variations.
        bs, bse = weighted_score(counts_bulk, -choices)
        cs, cse = weighted_score(counts_edge, choices)
        chance = strength*(m - float(law @ k))/v
        transmitted = int(rng.binomial(AREA_SHOTS, chance))
        arms.append(dict(bulk_counts=counts_bulk.tolist(), exterior_counts=counts_edge.tolist(),
                         tag_hits=transmitted, interior_score=bs, interior_se=bse,
                         exterior_score=cs, exterior_se=cse))
    bulk = (arms[1]['interior_score']-arms[0]['interior_score'])/(2*step)
    edge = (arms[1]['exterior_score']-arms[0]['exterior_score'])/(2*step)
    bse = math.hypot(*(a['interior_se'] for a in arms))/(2*step)
    return dict(class_transfer=k[idx].tolist(), displacement=step, arms=arms,
                exact_tag_slope=float(-direction @ k),
                exact_bulk_slope=float(-direction @ choices), measured_bulk_slope=bulk,
                measured_bulk_se=bse, measured_edge_slope=edge,
                measured_tag_slope=(arms[1]['tag_hits']-arms[0]['tag_hits'])*v/(AREA_SHOTS*strength*2*step))


def enroll(board, parent_cells, depth):
    divisions = 4**depth
    members = np.repeat(np.isin(np.arange(board.cells//divisions), list(parent_cells)), divisions)
    faces = board.external//12
    sides = members[faces]
    cross = sides[:, 0] ^ sides[:, 1]
    cut = board.external[cross][sides[cross]]
    if len(np.unique(cut)) != len(cut):
        raise ValueError('Two contact channels share an inner endpoint')
    return 12*int(members.sum()), np.sort(cut)


def energy_meter(v, cut):
    """Class-conditional mean of the unchanged parked-region free Hamiltonian."""
    m = len(cut)
    boundary = set(map(int, cut))
    cells_with_cut = {port//12 for port in boundary}
    cut_cut = sum(12*cell+a in boundary and 12*cell+b in boundary
                  for cell in cells_with_cut for a,b in CELL_SEAMS)
    cross = 5*m - 2*cut_cut
    contact_count = 30*(v//12) + (3*(v//12)-m)//2
    return dict(cut_cut=cut_cut, cut_interior=cross, interior_interior=contact_count-cut_cut-cross,
                degree_five_interior=9*(v//12))


def conditional_mean_energy(v, m, n, k, meter):
    """A quadratic in k, computed exactly on uniform sharp words in that class."""
    def unlike(size, occupied):
        return Fraction(2*occupied*(size-occupied), size*(size-1))
    u, z = Fraction(k,m), Fraction(n-k,v-m)
    unequal = (meter['cut_cut']*unlike(m,k) + meter['interior_interior']*unlike(v-m,n-k)
               + meter['cut_interior']*(u+z-2*u*z))
    contrast = 4*Fraction(candidate.KAPPA)
    return (contrast+Fraction(1,2))*unequal + contrast*meter['degree_five_interior']*z


def acquire(v, cut, host, n, filling, strength, seed):
    m = len(cut)
    k, law, choices = class_law(v, m, n)
    rng = np.random.default_rng(seed)
    volume_hits = int(rng.binomial(AREA_SHOTS, v/host))
    p_tag = strength * (m - law@k) / v
    transmitted = int(rng.binomial(AREA_SHOTS, p_tag))
    volume = host*volume_hits/AREA_SHOTS
    area = volume*transmitted/(AREA_SHOTS*strength)
    cap_counts = rng.multinomial(INFO_SHOTS, law)
    observed_capacity, cap_se = weighted_score(cap_counts, choices)
    f = n/v
    thermal = -f*math.log(f)-(1-f)*math.log1p(-f)
    capacity = float(law @ choices)
    if not -1e-9 <= m*thermal-capacity <= math.log(m+1)+1+1e-9:
        raise ArithmeticError('Conditional-word entropy bound failed')
    result = dict(filling=filling, v=v, m=m, n=n, volume_hits=volume_hits, tag_hits=transmitted,
                  volume_read=volume, available_transmission_read=area, population_transmission=m*(1-f),
                  capacity_read_nats=observed_capacity, capacity_se_nats=cap_se,
                  population_capacity_nats=capacity, first_class=int(k[0]), class_counts=cap_counts.tolist(),
                  boundary_word_examples=edge_examples(k, law, cut, rng))
    if filling == 'half':
        result['energy_meter'] = energy_meter(v,cut)
        result['matched_variation'] = tangent(v, m, n, k, law, choices, strength, rng)
    return result


def run(directory):
    hosts = {i: load_refinement(directory/f'FEDERATION_LEVEL{i}.npz', level=i) for i in range(8)}
    strength = tag_probability()
    panels = []
    settings = sorted({(0,i) for i in range(8)} | {(j,j+1) for j in range(4)})
    for parent, depth in settings:
        host = hosts[parent+depth]
        parent_shapes = selected_shapes(hosts[parent])
        results = []
        for idx, (name, cells) in enumerate(parent_shapes):
            v, cut = enroll(host, cells, depth)
            rows = []
            for fi,(filling,n) in enumerate(zip(FILLINGS, [2, v//4, v//2, 3*v//4], strict=True)):
                data = acquire(v, cut, 12*host.cells, n, filling, strength, [SEED,parent,depth,idx,fi])
                rows.append(data)
            results.append(dict(name=name, parent_cells=sorted(map(int,cells)), sources=rows))
        comparisons = []
        for cycle, path in zip(results[::2], results[1::2], strict=True):
            for a,b in zip(cycle['sources'], path['sources'], strict=True):
                if a['v']!=b['v'] or b['m']<=a['m']:
                    raise ValueError('Shape comparison failed equal-volume qualification')
                def contrast(label):
                    return b[label]-a[label]
                measured = contrast('capacity_read_nats')/contrast('available_transmission_read')
                exact = contrast('population_capacity_nats')/contrast('population_transmission')
                comparisons.append(dict(region=cycle['name'], filling=a['filling'],
                                        population_nats_per_available_contact=exact,
                                        observed_nats_per_available_contact=measured,
                                        frozen_half_fill_relative_error=abs(measured/(2*math.log(2))-1)))
        panels.append(dict(parent_level=parent, microscopic_subdivisions=depth,
                           host_cells=host.cells, geometries=results, shape_transfers=comparisons))
    return dict(seed=SEED, protocol_sha256=hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest(),
                carrier_hashes={i:b.sha256 for i,b in hosts.items()}, address_or_tag_trials=AREA_SHOTS,
                information_trials=INFO_SHOTS, isolated_tag_probability=strength, panels=panels)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--carriers', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    encoded = json.dumps(run(args.carriers), ensure_ascii=False, separators=(',', ':'), allow_nan=False)
    with args.output.open('x') as stream:
        stream.write(encoded+'\n')
