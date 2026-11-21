"""Finite success-and-failure preparation followed by transport through an OPH federation."""

import numpy as np
from scipy import sparse as sp
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.scalar_channel.carrier import require
from model.experiments.cage import construction as cell
from .source import chain
from .instruments import acquire, clock_amplitude, clock_read, count_partition, fuel_law


def run(spec, rng):
    menu = spec['flight']
    car = cell.federation(menu['federation'])
    edge = np.array(car['edges'])
    E, W = car['C'].shape
    a, b = next((a,b) for a,b in car['edges'] if a//12 != b//12)
    sender, receiver = a//12, b//12
    site = next(i for i,(s,t) in enumerate(car['edges']) if a in (s,t) and s//12 == t//12)
    # An intercell seam has its own outcome. A wedge belongs to its centre's cell.
    owners = np.r_[np.where(edge[:,0]//12 == edge[:,1]//12, edge[:,0]//12, -1),
                   [x[1]//12 for x in car['wedges']]]
    labels = np.where(owners == sender, 0, np.where(owners == receiver, 1, np.where(owners == -1, 2, 3)))
    B = sp.csr_matrix(car['C'], dtype=float)
    D = sp.bmat([[None,B], [B.T,None]], format='csr')  # generator in the oscillator's units
    reference_steps = np.arange(4*max(menu['waits'])+1)/4
    clock = acquire(rng, fuel_law(clock_amplitude(reference_steps), (0,2,1,3)), spec['shots'])
    receiver_times = clock_read(clock, spec['shots'])[np.array(menu['waits'])*4]
    rows = []
    for length in range(1, spec['depth']+2):
        H = chain(spec['depth'], length)/cell.HOP
        K = (np.eye(len(H))+expm(-1j*menu['filter_clock_interval']*H))/2
        psi = np.eye(len(H))[:,0]
        failures = []
        for _ in range(menu['filter_rounds']):
            out = K @ psi
            failures.append(float(np.vdot(psi, psi).real - np.vdot(out, out).real))
            psi = out
        values, modes = np.linalg.eigh(H)
        light = abs(values) <= 1+1e-12
        success = float(np.vdot(psi, psi).real)
        purity = float(np.linalg.norm(modes[:,light].T @ psi)**2 / success)
        source = np.zeros((len(H), E+W), complex)
        source[:,site] = psi
        grade = sp.diags(np.r_[np.ones(length), -np.ones(length)])
        generator = sp.kron(sp.csr_matrix(H), sp.eye(E+W), format='csr') + sp.kron(grade, D, format='csr')
        waits = menu['waits']
        laws, counts, stopped = [], [], []
        for t in waits:
            state = expm_multiply(-1j*t*generator, source.ravel(), traceA=0).reshape(source.shape)
            regions = np.bincount(labels, weights=(abs(state)**2).sum(axis=0), minlength=4)
            law = np.r_[failures, regions]
            require(law.min() >= -1e-13 and abs(law.sum()-1)<1e-11, 'incomplete flight readout')
            law = np.maximum(law, 0)
            laws.append(law.tolist())
            counts.append(count_partition(rng, law, spec['shots']))
            # Remove just the exterior hopping. The complete internal evolution remains.
            rest = expm(-1j*t*H) @ psi
            off = np.r_[failures, np.vdot(rest,rest).real, 0., 0., 0.]
            stopped.append(count_partition(rng, off, spec['shots']))
        rows.append(dict(length=length, success=success, light_weight_given_pass=purity,
                         laws=laws, counts=counts, stopped_outer_counts=stopped))
    return dict(sender_cell=int(sender), receiver_cell=int(receiver), seam_seed=list(map(int,edge[site])),
                source_modes=(E+W)*2*((4**(spec['depth']+1)-1)//3), clock_counts=clock.tolist(),
                receiver_times=receiver_times.tolist(), histories=rows)


def assess(spec, transit):
    menu = spec['flight']
    require([row['length'] for row in transit['histories']] == list(range(1, spec['depth']+2)),
            'incomplete or repeated transfer histories')
    for row in transit['histories']:
        require(0 < row['success'] <= 1 and
                menu['minimum_light_weight_given_pass'] <= row['light_weight_given_pass'] <= 1+1e-12,
                'filter admits too much bulk weight')
        laws = np.asarray(row['laws'])
        require(laws.shape == (len(menu['waits']), menu['filter_rounds']+4)
                and np.all(laws >= 0.) and np.allclose(laws.sum(axis=1),1.,atol=1e-11,rtol=0)
                and np.allclose(laws[:,-4:].sum(axis=1),row['success'],atol=1e-11,rtol=0),
                'invalid transfer probabilities')
        for key in ('counts','stopped_outer_counts'):
            counts = np.asarray(row[key])
            require(counts.shape == (len(menu['waits']), menu['filter_rounds']+4) and counts.dtype.kind in 'iu'
                    and np.all(counts >= 0) and np.all(counts.sum(axis=1)==spec['shots']), 'lost flight outcomes')
        neighbor = np.asarray(row['counts'])[:,-3]/spec['shots']
        require(neighbor[0] == 0 and neighbor.max() > menu['minimum_unconditional_neighbor_probability']
                and not np.any(np.asarray(row['stopped_outer_counts'])[:,-3:]), 'flight control failed')
