"""Clock every shortest protected tag route, and test Euclidean volume and area on the same source."""
import argparse
from itertools import combinations, islice
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.spatial import ConvexHull

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ..autonomous_selection.measure import stationary_pair_tangent
from ..equilibrium_refinement.measure import apparatus, matched_energy_and_load
from ..shape_stationarity.geometry import selected_shapes
from ..shape_stationarity.measure import native_sector, partition

HERE = Path(__file__).resolve().parent
SHOTS = 2**38
SEED = 351729


def pulse_calibration():
    law = candidate.native_law
    endpoint = np.eye(law.N_SRC)[:,law.SRC.index((1,0))]
    receiver = law.SRC.index((0,1))
    h = law.seam_term(candidate.KAPPA)
    forward = expm(-1j*np.pi*h) @ endpoint
    round_trip = expm(-1j*np.pi*h) @ forward
    np.testing.assert_allclose(abs(forward[receiver])**2, 1, atol=1e-13)
    np.testing.assert_allclose(abs(round_trip @ endpoint)**2, 1, atol=1e-13)
    # A second isolated contact with a single load acts as the material hand.
    # On fresh copies a 3π/2 analyzer delay turns the cosine into positive sine.
    times = np.linspace(0,np.pi,5)
    source = law.SRC.index((1,0))
    prob = np.array([[abs((expm(-1j*(t+delay)*h) @ endpoint)[source])**2 for t in times]
                     for delay in (0,1.5*np.pi)])
    hand = np.unwrap(np.arctan2(2*prob[1]-1,2*prob[0]-1))
    duration = abs(hand[-1]-hand[0])
    np.testing.assert_allclose(duration, np.pi, atol=2e-13)
    return dict(forward_probability=float(abs(forward[receiver])**2),
                return_probability=float(abs(round_trip @ endpoint)**2),
                quadrature_probabilities=prob.tolist(), decoded_single_step=float(duration))


def ambient_reduction(car, ports, borders):
    """Keep every excursion outside the region, eliminating only unoccupied external vertices."""
    A = car.adjacency().tocsr()
    outer = np.setdiff1d(np.arange(car.n), ports)
    local = A[ports][:,ports].toarray()
    points = sorted(borders)
    mask = np.array([np.searchsorted(ports,a) for a in points])
    outside_gate = [np.searchsorted(outer,borders[p]) for p in points]
    between = shortest_path(A[outer][:,outer], indices=outside_gate, directed=False,
                            unweighted=True, method='D')[:,outside_gate] + 2
    for i,j in combinations(range(len(points)),2):
        a,b = mask[i],mask[j]
        if np.isfinite(between[i,j]) and (local[a,b] == 0 or between[i,j] < local[a,b]):
            local[a,b] = local[b,a] = between[i,j]
    return local


def distances(local, resident):
    retained = np.setdiff1d(np.arange(len(local)),resident)
    table = shortest_path(local[np.ix_(retained,retained)], directed=False, method='FW')
    return retained, table


def measured_hull(d):
    """Convex-hull reconstruction in a three-coordinate classical distance scaling."""
    n = len(d)
    squared = d*d
    mean = squared.mean(axis=1)
    gram = -(squared - mean[:,None] - mean[None,:] + mean.mean())/2
    val, vec = np.linalg.eigh(gram)
    coords = vec[:,-3:]*np.sqrt(np.maximum(val[-3:],0))
    delta = coords[:,None,:]-coords[None,:,:]
    ranges = np.sqrt((delta**2).sum(axis=2))
    error = np.linalg.norm(ranges-d)/np.linalg.norm(d) if n>1 else 0.
    if n < 4 or val[-3] < 1e-9:
        return 0.,0.,float(error)
    hull = ConvexHull(coords)
    return float(hull.volume),float(hull.area),float(error)


def quartet_bound(d, indices):
    """Any four Euclidean points obey the quadrilateral inequality for every ordering."""
    a,b,c,e = indices.T
    paired = np.column_stack([d[a,b]**2+d[c,e]**2,
                              d[a,c]**2+d[b,e]**2, d[a,e]**2+d[b,c]**2])
    ratio = np.sqrt(paired.max(axis=1)/(paired.sum(axis=1)-paired.max(axis=1)))
    eps = np.maximum(0, (ratio-1)/(ratio+1))
    idx = eps.argmax()
    return float(eps[idx]), indices[idx].tolist()


def certificates(d, ports):
    result = [(0., []) for _ in range(2)]
    sets = combinations(range(len(d)),4)
    while batch := list(islice(sets,65536)):
        quads = np.array(batch)
        full = quartet_bound(d,quads)
        result[0] = max(result[0],full)
        owners = np.sort(ports[quads]//12,axis=1)
        distinct = np.all(np.diff(owners,axis=1)>0,axis=1)
        if np.any(distinct):
            result[1] = max(result[1],quartet_bound(d,quads[distinct]))
    return [dict(min_relative_distortion=w, ports=ports[q].tolist()) for w,q in result]


def features(d):
    components = []
    seen = set()
    for row in range(len(d)):
        if row in seen:
            continue
        ids = np.flatnonzero(np.isfinite(d[row]))
        seen.update(ids)
        components.append(measured_hull(d[np.ix_(ids,ids)]))
    a = np.array(components)
    return [float(a[:,0].sum()),float(a[:,1].sum()),float(a[:,2].max()),len(a)]


def distribution_derivatives(words, cut, energy):
    interior_ids = partition(words,cut)[0]
    k = np.array([sum(v in cut for v in word) for word in words])
    matched = matched_energy_and_load(energy,k,interior_ids)
    base = np.array([(k==j)/np.count_nonzero(k==j) for j in range(3)])
    strict = np.array(matched['sector_transfers']) @ base
    strict[matched['high_word']] += matched['compensation']
    strict[matched['low_word']] -= matched['compensation']
    pert = np.stack([base[1]-base[0],base[2]-base[0],strict])
    np.testing.assert_allclose(pert.sum(axis=1), 0, atol=1e-13)
    return pert


def sample_direction(direction, measures, rng):
    w = np.full(len(direction), 1/len(direction))
    amplitude = .7 / (len(w)*np.abs(direction).max())
    groups, labels = np.unique(measures[:,2],return_inverse=True)
    values = []
    variances = []
    counts = []
    for sign in (-1,1):
        p = w + sign*amplitude*direction
        # Separate fresh copies feed the range survey and the cut-count reader.
        geometry = rng.multinomial(SHOTS,p)
        edge = rng.multinomial(SHOTS,np.bincount(labels,weights=p))
        geo_mean = geometry @ measures[:,:2] / SHOTS
        edge_mean = edge @ groups / SHOTS
        values.append(np.r_[geo_mean,edge_mean])
        variances.append(np.r_[geometry @ ((measures[:,:2]-geo_mean)**2) / SHOTS,
                               edge @ ((groups-edge_mean)**2) / SHOTS])
        counts.append(dict(geometry=geometry.tolist(),cut_classes=edge.tolist()))
    change = (values[1]-values[0])/(2*amplitude)
    se = np.sqrt(np.sum(variances,axis=0)/SHOTS)/(2*amplitude)
    return dict(mixture_displacement=float(amplitude), population=(direction @ measures).tolist(),
                sampled=change.tolist(), standard_errors=se.tolist(), cohort_counts=counts)



def measure_case(car, level, name, cells, host, rng):
    ports, contacts, cut, words = apparatus(car,cells)
    ports = np.array(ports)
    local = ambient_reduction(car,ports,cut)
    _, blank = distances(local,[])
    np.testing.assert_array_equal(blank,host[np.ix_(ports,ports)])
    lower, coarse = certificates(blank,ports)
    table = []
    word_indices = np.searchsorted(ports,np.array(words))
    # Each cell has at least one intact icosahedron for two residents. Its
    # internal two-step certificate holds in the ambient metric: no external
    # seam joins two ports of the same cell, and every port has <=1 external seam.
    same_cell = np.array(list(combinations(range(12),4)))
    regional_min = []
    for owner in sorted(cells):
        ids = np.flatnonzero(ports//12 == owner)
        regional_min.append(quartet_bound(blank[np.ix_(ids,ids)],same_cell)[0])
    guaranteed = min(regional_min)
    for resident in word_indices:
        retained, d = distances(local,resident)
        table.append(features(d))
    for i in (0, len(words)//2,len(words)-1):
        resident = np.array(words[i])
        active = np.setdiff1d(np.arange(car.n),resident)
        kept = np.setdiff1d(ports,resident)
        region = np.searchsorted(active,kept)
        full = car.adjacency().tocsr()[active][:,active]
        direct = shortest_path(full, method='D', directed=False, indices=region, unweighted=True)[:,region]
        _, reduced = distances(local,word_indices[i])
        np.testing.assert_array_equal(direct,reduced)
    energy, _ = native_sector(car,contacts,words)
    k = np.array([sum(v in cut for v in w) for w in words])
    capacity = np.log([comb(len(cut),int(kk)) for kk in k])
    pert = distribution_derivatives(words,cut,energy)
    table = np.asarray(table)
    measures = np.column_stack([table[:,:2],capacity])
    basis_grad = pert @ measures
    v = basis_grad[:,0]
    a = basis_grad[:,1]
    vhat = v / np.linalg.norm(v)
    fixed_volume = a - vhat * (vhat @ a)
    fixed_both = np.cross(v,a)
    dirs = {}
    for label_,c in [('area_at_fixed_volume',fixed_volume),('both_geometry_fixed',fixed_both)]:
        t = c @ pert
        if np.max(abs(t)) < 1e-11:
            dirs[label_] = None
            continue
        t /= np.maximum(t,0).sum()
        if t[np.argmax(abs(t))]<0:
            t *= -1
        dirs[label_] = sample_direction(t,measures,rng)
    free,_,_ = stationary_pair_tangent(car,cells,words)
    return dict(level=level,shape=name,cells=sorted(cells),port_order=ports.tolist(),
                component_shapes_of_uniform_word_survey=list(map(int, sorted(set(table[:,3])))),
                blank_diameter=int(blank.max()), fraction_of_host_diameter=float(blank.max()/host.max()),
                blank_vol_area_range_error=features(blank)[:3],
                blank_quartet=lower,blank_four_cell_quartet=coarse,
                all_resident_words_distortion_lower_bound_from_intact_cell=guaranteed,
                min_word_relative_range_error=float(table[:,2].min()),
                max_word_relative_range_error=float(table[:,2].max()),
                uniform_mean_volume_area_capacity=measures.mean(axis=0).tolist(),
                columns=['volume','area','exterior_capacity_nats'],
                basis_slopes=basis_grad.tolist(), fixed_volume_experiments=dirs,
                normalized_free_stationary_transfer=(free @ measures).tolist(),
                per_word_volume_area_rms=table[:,:3].tolist())


def run():
    pulse = pulse_calibration()
    rng = np.random.default_rng(SEED)
    boards = []
    cases=[]
    for level in range(3):
        car = federation(level)
        host = shortest_path(car.adjacency(),method='D',directed=False,unweighted=True)
        diam = int(host.max())
        boards.append(dict(level=level,host_range_diameter=diam,
                           mean_ball_count=[float(np.mean(np.sum(host<=i,axis=1)))
                                            for i in range(diam+1)]))
        for name,cells in selected_shapes(car):
            cases.append(measure_case(car,level,name,cells,host,rng))
    first = cases[0]['fixed_volume_experiments']['area_at_fixed_volume']['population']
    eta = first[2]/first[1]
    for c in cases:
        for conf in c['fixed_volume_experiments'].values():
            if conf:
                conf['capacity_minus_calibrated_area'] = conf['population'][2] - eta*conf['population'][1]
    return dict(seed=SEED,shots_per_signed_cohort=SHOTS,pulse=pulse,boards=boards,
                training_nats_per_hull_area=eta,cases=cases)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    target=parser.parse_args().output
    result=run()
    with target.open('x') as out:
        json.dump(result,out,separators=(',',':'),allow_nan=False)
        out.write('\n')
