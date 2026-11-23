"""Charge flight and the regional quantum entropy derivative of a written-ring vacuum."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.electron_field.flux_string import search
from model.experiments.interval_echo_einstein_station.written_field_records.measure import ring_laplacian
from .distance import distances

HERE = Path(__file__).resolve().parent
TIMES = np.arange(1, 4) / 16
SIGNS = (-1, 1)
SHOTS = 2**32
SEED = 359611


def shapes(patch):
    graph = nx.Graph(patch.edges)
    chosen = []
    for n in (3, 4, 5):
        cuts = [(len(list(nx.edge_boundary(graph, v))), v)
                for v in combinations(range(12), n) if 0 in v and nx.is_connected(graph.subgraph(v))]
        for seam_count in (min(cuts)[0], max(cuts)[0]):
            chosen.append(min(v for m, v in cuts if m == seam_count))
    return chosen


def charge_sector(patch, words, anchor, sign):
    """Extend the anchored charge fiber by every allowed ring flip and signed transfer."""
    basis = [(anchor, int(w)) for w in words]
    index = {state: i for i, state in enumerate(basis)}
    rows, cols, rates = [], [], []
    for i, (pos, word) in enumerate(basis):
        moves = [((pos, other), .5) for _, other in patch.ring_moves(word)]
        moves.extend((state, 4.) for state in patch.charge_moves(pos, word, sign))
        for state, rate in moves:
            if state not in index:
                index[state] = len(basis)
                basis.append(state)
            rows.append(i)
            cols.append(index[state])
            rates.append(rate)
    adj = coo_matrix((rates, (rows, cols)), shape=(len(basis),)*2).tocsr()
    if (adj-adj.T).nnz:
        raise ArithmeticError("Local signed moves lost Hermiticity")
    positions = np.array([p for p, _ in basis])
    degree = search.federation(0).degree[:12]
    h = diags(np.asarray(adj.sum(axis=1)).ravel() + 8*(6-degree[positions]))-adj
    return positions, h


def launched_tables(patch, words, field, derivatives):
    """For every signed anchor, retain the twelve terminal-port outcomes."""
    columns = field.shape[1]
    outcomes = np.zeros((len(SIGNS), 12, len(TIMES), 12, columns))
    sizes = np.zeros((2,12), dtype=int)
    for si, sign in enumerate(SIGNS):
        for anchor in range(12):
            patch.compensator(anchor, sign)
            ports, h = charge_sector(patch, words, anchor, sign)
            sizes[si, anchor] = len(ports)
            starts = np.zeros((len(ports), columns), dtype=complex)
            starts[:len(words)] = field
            for t, states in enumerate(expm_multiply(-1j*h, starts, start=0, stop=TIMES[-1],
                                                      num=4, traceA=-1j*h.diagonal().sum())):
                if not t:
                    continue
                amplitudes = abs(states)**2
                if derivatives:
                    amplitudes[:,1:] = 2*np.real(states[:,:1].conj()*states[:,1:])
                for endpoint in range(12):
                    outcomes[si, anchor, t-1, endpoint] = amplitudes[ports==endpoint].sum(axis=0)
    norm = outcomes.sum(axis=3)
    target = np.ones(columns)
    if derivatives:
        target[1:] = 0
    np.testing.assert_allclose(norm, np.broadcast_to(target, norm.shape), atol=1e-12, rtol=0)
    return outcomes, sizes


def crossing_laws(outcomes, vertices):
    inside = np.isin(np.arange(12), vertices)
    rows = []
    for launched, received in ((inside, ~inside), (~inside, inside)):
        rows.append(outcomes[:,launched][:,:,:,received].sum(axis=3).mean(axis=1))
    return np.asarray(rows).reshape(12, outcomes.shape[-1])


def contact_survey(patch, words, vertices):
    inside = set(vertices)
    results = []
    for sign in SIGNS:
        row = []
        for origins in (inside, set(range(12))-inside):
            crossed = np.array([sum((patch.move(p, int(w), edge, sign)[0] not in origins)
                                    for p in origins for edge in patch.edges if p in edge) for w in words])
            if np.ptp(crossed):
                raise ArithmeticError("Gate census depends on a legal field word")
            row.append(int(crossed[0]))
        results.append(row)
    counts = np.array(results)
    p_out = counts[:,0].mean()/(5*len(inside))
    p_in = counts[:,1].mean()/(5*(12-len(inside)))
    return dict(enabled_out_in_counts_by_sign=results, accessible_ports=12*p_in/(p_in+p_out),
                calibrated_area_contacts=int(counts.sum(axis=0)[0]))


def entropy(p):
    p = np.asarray(p)
    return float(-p[p>0] @ np.log(p[p>0]))


class Region:
    def __init__(self, patch, words, vertices):
        internal = [e for e in patch.edges if set(e).issubset(vertices)]
        mask = sum(patch.bit[e] for e in internal)
        _, self.rows = np.unique(words & mask, return_inverse=True)
        _, self.cols = np.unique(words & ~mask, return_inverse=True)
        flux = np.zeros((len(words), 12), dtype=int)
        for a,b in internal:
            head = ((words & patch.bit[a,b]) != 0).astype(int)
            flux[:,a] += head
            flux[:,b] -= head
        self.flux, self.center = np.unique(flux, axis=0, return_inverse=True)
        self.center_counts = np.bincount(self.center)
        self.center_info = -np.log(self.center_counts / len(words))
        self.interior_info = -np.log(np.bincount(self.rows) / len(words))
        self.completions = np.log(np.bincount(self.rows))
        self.complete = all(len(set(self.rows[self.center==k]))*len(set(self.cols[self.center==k]))==n
                            for k,n in enumerate(self.center_counts))
        self.exterior_reads_center = all(len(set(self.center[self.cols==k]))==1
                                        for k in range(self.cols.max()+1))
        if not self.exterior_reads_center:
            raise ValueError("The exterior heads fail to identify this flux sector")
        self.exterior_to_center = np.zeros(self.cols.max()+1,dtype=int)
        self.exterior_to_center[self.cols] = self.center
        g = np.ones(len(words))/np.sqrt(len(words))
        density = self.density(g)
        values, basis = np.linalg.eigh(density)
        active = values > 1e-13
        self.schmidt = values[active]
        self.eigenbasis = basis[:, active]
        self.modular = -np.log(self.schmidt)
        self.modal_scores = np.r_[self.modular, 0.]
        self.vacuum_entropy = entropy(self.schmidt)

    def matrix(self, wave):
        coefficients = np.zeros((self.rows.max()+1,self.cols.max()+1), dtype=complex)
        coefficients[self.rows,self.cols] = wave
        return coefficients

    def density(self, wave):
        amplitudes = self.matrix(wave)
        return amplitudes @ amplitudes.conj().T

    def modular_record(self, wave):
        projected = self.eigenbasis.conj().T @ self.matrix(wave)
        probs = (abs(projected)**2).sum(axis=1)
        probs = np.r_[probs,max(0,1-sum(probs))]
        return probs/probs.sum()

    def gradients(self, base, variations):
        weights = 2*np.real(base[:,None].conj()*variations)
        coeff = self.matrix(base)
        scored = self.eigenbasis @ (self.modular[:,None]*(self.eigenbasis.conj().T@coeff))
        return dict(modular=2*np.real(scored[self.rows,self.cols].conj()@variations),
                    center=self.center_info[self.center]@weights,
                    interior_word=self.interior_info[self.rows]@weights,
                    exterior_completion=self.completions[self.rows]@weights)


def blind_direction(crossings, bulk):
    _, sigma, vt = np.linalg.svd(crossings, full_matrices=False)
    rank = int(sum(sigma > max(1e-13,1e-11*sigma[0])))
    projection = bulk - vt[:rank].T @ (vt[:rank] @ bulk)
    if np.linalg.norm(projection) < 1e-12:
        return np.zeros_like(projection), rank, sigma
    return projection/np.linalg.norm(projection), rank, sigma


def signed_read(counts, score, theta):
    counts = np.array(counts)
    p = counts/counts.sum(axis=1)[:,None]
    estimates = p@score
    variances = p@(score**2) - estimates**2
    return dict(measured_secant=float(np.diff(estimates)[0]/(2*theta)),
                standard_error=float(np.sqrt(np.maximum(variances,0).sum()/SHOTS)/(2*theta)))


def sample_entropy(region, waves, theta, rng):
    center_laws, modal_laws, entropies = [], [], []
    for wave in waves.T:
        exterior_law = np.bincount(region.cols,abs(wave)**2)
        center_laws.append(np.bincount(region.exterior_to_center,exterior_law))
        modal_laws.append(region.modular_record(wave))
        values = np.linalg.eigvalsh(region.density(wave))
        entropies.append(entropy(values[values>1e-14]))
    result = {}
    for name, laws, scores in [('modular',modal_laws,region.modal_scores),
                                ('center',center_laws,region.center_info)]:
        laws = np.array(laws)
        counts = [rng.multinomial(SHOTS,x/x.sum()).tolist() for x in laws]
        result[name] = dict(counts=counts, expected_secant=float(np.diff(laws@scores)[0]/(2*theta)),
                            **signed_read(counts,scores,theta))
    result.update(von_neumann_entropies=entropies, von_neumann_secant=float(np.diff(entropies)[0]/(2*theta)),
                  null_support_probabilities=np.array(modal_laws)[:,-1].tolist())
    return result


def endpoint_calibration(gradient, tangent, entropy_gradient):
    response = gradient @ tangent
    tv = np.sum(abs(response), axis=3)/2
    setting = np.unravel_index(np.argmax(tv), tv.shape)
    si, origin, time = map(int,setting)
    decoder = (response[setting] > 1e-13).astype(float)
    forecasts = {}
    shared = gradient.mean(axis=0)
    reciprocal = (shared + np.swapaxes(shared,0,2))/2
    for name, matrix in (('complete_signed',gradient), ('common_reciprocal',reciprocal)):
        matrix = matrix.reshape(-1, len(tangent))
        u, singular, vt = np.linalg.svd(matrix, full_matrices=False)
        rank = int(sum(singular > max(1e-13,1e-11*singular[0])))
        error = entropy_gradient - vt[:rank].T @ (vt[:rank] @ entropy_gradient)
        forecasts[name] = dict(rank=rank, singular_values=singular.tolist(),
                               unresolved_fraction=float(np.linalg.norm(error)/np.linalg.norm(entropy_gradient)),
                               largest_entry_slope=float(np.max(abs(matrix@tangent))))
    return dict(sign=SIGNS[si],origin=origin,exposure=float(TIMES[time]),
                decoder_ports=np.flatnonzero(decoder).tolist(),
                terminal_tv_derivative=float(tv[setting]),
                decoder_slope=float(response[setting]@decoder), forecasts=forecasts),setting,decoder


def held_endpoints(tables, theta, pair, setting, decoder, rng):
    counts = np.zeros(tables.shape[:-1]+(2,), dtype=np.int64)
    for index in np.ndindex(tables.shape[:3]):
        for branch, probabilities in enumerate(tables[index][:,pair].T):
            counts[index+(slice(None),branch)] = rng.multinomial(SHOTS,probabilities/probabilities.sum())
    laws = tables[setting][:,pair].T
    read = signed_read(counts[setting].T,decoder,theta)
    return dict(full_endpoint_counts=counts.tolist(),
                selected_held_tv=float(np.abs(np.diff(laws,axis=0)).sum()/2),
                selected_expected_secant=float(np.diff(laws@decoder)[0]/(2*theta)),
                **read)


def run():
    patch = search.Patch()
    words = np.array(patch.field_component(), dtype=np.int64)
    g = np.ones(len(words))/np.sqrt(len(words))
    h = ring_laplacian(patch, words)
    np.testing.assert_allclose(h@g,0,atol=1e-13)
    z = np.column_stack([np.where((words&patch.bit[e])!=0,-1,1) for e in patch.edges])
    tangent = expm_multiply(-1j*h,1j*g[:,None]*z)
    tables, sizes = launched_tables(patch,words,np.column_stack([g,tangent]), True)
    original = json.loads((Path(search.__file__).parent/'results.json').read_text())
    anchor = patch.faces[0][0]
    for si,sign in enumerate(SIGNS):
        held = next(x for x in original['release'] if x['dimension']==sizes[si,anchor])
        at = next(x for x in held['records'] if x['time']==.125)
        np.testing.assert_allclose(tables[si,anchor,1,:,0],at['uniform']['charge_port_probabilities'],
                                   atol=2e-12,rtol=0)
    reports = []
    splits = []
    for vertices in shapes(patch):
        region = Region(patch,words,vertices)
        slopes = region.gradients(g,tangent)
        geometry = crossing_laws(tables,vertices)
        direction, rank, sigma = blind_direction(geometry[:,1:],slopes['modular'])
        reports.append(dict(vertices=vertices,contact=contact_survey(patch,words,vertices),
                            rectangular_flux_sectors=region.complete,
                            exterior_reads_center=region.exterior_reads_center,
                            sector_populations=region.center_counts.tolist(),
                            vacuum_schmidt_entropy=region.vacuum_entropy,
                            vacuum_center_entropy=entropy(region.center_counts/len(words)),
                            vacuum_schmidt=region.schmidt.tolist(),
                            derivative_by_head={k:v.tolist() for k,v in slopes.items()},
                            free_crossing_vacuum_and_derivatives=geometry.tolist(),
                            singular_values=sigma.tolist(),rank=rank,
                            fitted_phase_direction=direction.tolist(),
                            blind_entropy_slope=float(slopes['modular']@direction),
                            max_blind_crossing_slope=float(max(abs(geometry[:,1:]@direction)))))
        splits.append(region)
    winner=int(np.argmax([r['blind_entropy_slope'] for r in reports]))
    selected=reports[winner]
    direction=np.array(selected['fitted_phase_direction'])
    vacuum_ports=tables[:,:,:,:,0].copy()
    endpoint, endpoint_setting, decoder=endpoint_calibration(
        tables[:,:,:,:,1:], direction, np.array(selected['derivative_by_head']['center']))
    phase=z@direction
    angles=[1/64,1/128]
    signed_angles=np.ravel([[-t,t] for t in angles])
    starts=g[:,None]*np.exp(1j*np.outer(phase,signed_angles))
    field=expm_multiply(-1j*h,starts)
    tables,_=launched_tables(patch,words,field,False)
    crossing=crossing_laws(tables,selected['vertices'])
    rng=np.random.default_rng(SEED)
    held=[]
    for k,theta in enumerate(angles):
        pair=slice(2*k,2*k+2)
        protocols=[]
        for probs in crossing[:,pair]:
            clicks=rng.binomial(SHOTS,probs)
            trials=np.column_stack([SHOTS-clicks,clicks])
            protocols.append(dict(crossings=clicks.tolist(),probabilities=probs.tolist(),
                                  expected_secant=float(np.diff(probs)[0]/(2*theta)),
                                  **signed_read(trials,np.array([0.,1.]),theta)))
        held.append(dict(angle=theta, geometry=protocols,
                         endpoints=held_endpoints(tables,theta,pair,endpoint_setting,decoder,
                                                  np.random.default_rng([SEED,2,k])),
                         energy=np.sum(field[:,pair].conj()*(h@field[:,pair]),axis=0).real.tolist(),
                         entropy=sample_entropy(splits[winner],field[:,pair],theta,rng)))
    # A sharp word and its freely propagated copies retain their ground overlap;
    # ground energetic selection by a positive graph Laplacian is independent of attraction.
    sharp=np.zeros(len(words)); sharp[words.tolist().index(patch.seed)]=1.
    launches=[]
    for time in (0,1,4,16):
        prepared=expm_multiply(-1j*time*h,sharp) if time else sharp
        local=splits[winner].modular_record(prepared)
        counts=rng.multinomial(SHOTS,local)
        sharp_law=np.bincount(splits[winner].rows,abs(prepared)**2)
        launches.append(dict(time=time,ground_probability=float(abs(g@prepared)**2),
                             expected_energy=float(np.real(np.vdot(prepared,h@prepared))),
                             modular_score=float(local@splits[winner].modal_scores),
                             modular_counts=counts.tolist(),
                             regional_null_weight=float(local[-1]),
                             sharp_counts=rng.multinomial(SHOTS,sharp_law/sharp_law.sum()).tolist()))
    return dict(protocol_sha256=hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest(),
                endpoint_protocol_sha256=hashlib.sha256((HERE/'ENDPOINT_PROTOCOL.md').read_bytes()).hexdigest(),
                distance_protocol_sha256=hashlib.sha256((HERE/'DISTANCE_PROTOCOL.md').read_bytes()).hexdigest(),
                distance_charts=distances(tables,vacuum_ports,held),
                endpoint_calibration=endpoint, field_dimension=len(words),seed=SEED,shots=SHOTS,
                observable_phase_rank=int(np.linalg.matrix_rank(z-z.mean(axis=0),tol=1e-10)),
                sign_order=SIGNS,crossing_order='out then in, vacancy then excess, three increasing exposures',
                exposures=TIMES.tolist(), calibration_crossing_probabilities=np.sin(4*TIMES)**2,
                contact_clock_counts=rng.binomial(SHOTS,np.sin(4*TIMES)**2).tolist(),
                compensated_sector_dimensions=sizes.tolist(),head_order=patch.edges,
                regions=reports,selected_region=winner,held=held,free_sharp=launches,
                vacuum_sharp_counts=rng.multinomial(SHOTS,np.bincount(splits[winner].rows)/len(words)).tolist(),
                vacuum_modular_counts=rng.multinomial(SHOTS,splits[winner].modular_record(g)).tolist())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    record=run()
    with args.output.open('x') as file:
        json.dump(record,file,indent=2,allow_nan=False,default=lambda a:a.tolist())
        file.write('\n')
