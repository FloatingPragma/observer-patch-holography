"""Ground-energy and reference-arm records for one freely hopping dressing load."""

import argparse
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier

from . import measure, spectrum

SEED = 220401
SHOTS = 2**24
ARMS = {'ABC': (True, True, True), 'A': (True, False, False),
        'BC': (False, True, True), 'empty': (False, False, False)}


def geometry(board, graph, ports, owners):
    live = np.array([p for p in range(board.n) if p // 12 not in owners])
    mobile = board.laplacian()[live][:, live] / 2
    index = {int(p): i for i, p in enumerate(live)}
    bridges = []
    for s, owner in enumerate(owners):
        for other in sorted(set(graph[owner]) - set(owners)):
            bridges.append((s, ports[owner, other] % 12, index[ports[other, owner]]))
    return mobile, live, bridges


def ground(h):
    # A fixed nonsymmetric start excites each permitted sector.
    rng = np.random.default_rng(1220)
    evals, vecs = eigsh(h, k=2, which='SA', v0=rng.normal(size=h.shape[0]),
                        ncv=45, tol=5e-12, maxiter=18000)
    v = vecs[:, 0]
    residual = float(np.linalg.norm(h @ v - evals[0]*v))
    if residual > 1e-8:
        raise ArithmeticError('Unresolved mobile minimum')
    return float(evals[0]), v, float(evals[1]-evals[0]), float(1e-10*np.ceil(residual*1e10))


def phase_read(energy, rng):
    phase = np.exp(-1j*energy)
    prob = np.clip([(1+phase.real)/2, (1+phase.imag)/2], 0, 1)
    counts = rng.binomial(SHOTS, prob.round(13))
    x, y = 2*counts/SHOTS-1
    angle = -np.arctan2(y, x)
    # On this receipt the shifted mobile energies lie in (-pi,pi).
    if abs(energy) >= np.pi:
        raise ArithmeticError('Clock requires unwrapping')
    se = np.sqrt((y*y*(1-x*x) + x*x*(1-y*y))/SHOTS) / (x*x+y*y)
    return dict(x_y_counts=list(map(int, counts)), energy=angle, se=se)


def source_embedding(label):
    e = np.eye(3)[:, :1]
    wedge = spectrum.antitriplet_embedding()[:9, :3]
    return np.kron(np.eye(3) if ARMS[label][0] else e,
                   wedge if ARMS[label][1] else np.kron(e, e))


def full_source_h(graph, ports, owners, frames, mobile, bridges, k, label):
    embed = source_embedding(label)
    base, contact = spectrum.fixed_source_matrices(graph, ports, owners, frames, k, ARMS[label])
    hs = embed.T @ (base+contact) @ embed
    hs -= np.eye(len(hs))*np.linalg.eigvalsh(hs)[0]
    hs[abs(hs)<1e-13] = 0
    sites = mobile.shape[0]
    h = sparse.kron(sparse.eye(len(hs)), mobile, format='csr')
    h += sparse.kron(hs, sparse.eye(sites), format='csr')
    for i, p, n in bridges:
        if ARMS[label][i]:
            ops = [np.eye(3)]*3
            ops[i] = np.outer(frames[i][p], frames[i][p])
            charge = embed.T @ spectrum.tensor(ops) @ embed
            charge[abs(charge)<1e-13] = 0
            site = sparse.csr_matrix(([1.], ([n], [n])), shape=(sites, sites))
            h -= (1+2*k)*sparse.kron(charge, site, format='csr')
    return h


def trial(graph, ports, owners, frames, geom, k, mode, rng):
    mobile, live, bridges = geom
    arms = {}
    for label, active in ARMS.items():
        if mode == 'singlet_clamp':
            energy_on_port = np.zeros(len(live))
            for i, p, n in bridges:
                if active[i]:
                    energy_on_port[n] -= (1+2*k)*np.linalg.norm(frames[i][p])**2 / 3
            h = mobile + sparse.diags(energy_on_port)
        else:
            h = full_source_h(graph, ports, owners, frames, mobile, bridges, k, label)
        en, v, gap, residual = ground(h)
        weights = (v.reshape(-1, len(live))**2).sum(axis=0)
        boundary = [float(sum(weights[n] for j, _, n in bridges if i == j)) for i in range(3)]
        arms[label] = dict(ground=en, next_gap=gap, residual_upper_bound=residual, boundary_A_B_C=boundary,
                           phase_read=phase_read(en, rng))
    def calc(key):
        v = {label: arm['ground'] if key == 'exact' else arm['phase_read']['energy']
             for label, arm in arms.items()}
        # The isolation reference is chosen by the independent isolated-source spectra.
        best = min(('A', 'BC'), key=lambda label: arms[label]['ground'])
        return dict(connected=v['ABC']-v['A']-v['BC']+v['empty'],
                    best_source=best, beyond_best=v['ABC']-v[best])
    obs = calc('observed')
    errors = {label: arm['phase_read']['se'] for label, arm in arms.items()}
    obs['connected_se'] = np.sqrt(sum(x*x for x in errors.values()))
    obs['beyond_best_se'] = np.hypot(errors['ABC'], errors[obs['best_source']])
    return dict(arms=arms, exact=calc('exact'), observer=obs)


def slopes(rows):
    summaries = {}
    for name in rows[0]['trials']:
        records = [row['trials'][name]['observer'] for row in rows]
        first, second = records[:2]
        slope = second['beyond_best']-first['beyond_best']
        se = np.hypot(first['beyond_best_se'], second['beyond_best_se'])
        held = []
        for d, read in enumerate(records[2:], 3):
            prediction = first['beyond_best'] + (d-1)*slope
            err = np.sqrt(read['beyond_best_se']**2 + ((d-1)*second['beyond_best_se'])**2
                          + ((2-d)*first['beyond_best_se'])**2)
            held.append(dict(distance=d, forecast=prediction, observed=read['beyond_best'],
                             residual_z=(read['beyond_best']-prediction)/err))
        summary = dict(train_slope=slope, train_se=se, held=held)
        if len(rows) >= 10:
            # Descriptive after acquisition; the held contact-line test above was preregistered.
            late = records[5:]
            design = np.column_stack((np.ones(len(late)), np.arange(6, len(rows)+1)))
            w = np.diag([1/read['beyond_best_se']**2 for read in late])
            cov = np.linalg.inv(design.T @ w @ design)
            beta = cov @ design.T @ w @ [read['beyond_best'] for read in late]
            summary['posthoc_tail'] = dict(distances=list(range(6, len(rows)+1)),
                                           slope=beta[1], se=np.sqrt(cov[1,1]))
        summaries[name] = summary
    return summaries


def run():
    band = spectrum.physical_frame(currents._service(carrier.cell(), False, 1.)[4][0])
    previous = json.loads((measure.HERE / 'records.json').read_text())['levels']
    data = {}
    rng = np.random.default_rng(SEED)
    for lev in (0, 1):
        board = carrier.federation(lev)
        graph, ports = measure.skeleton(board)
        layout = previous[str(lev)]['colour']
        b, c = layout['antitriplet_owners']
        first = layout['rows'][0]['owner']
        frames = measure.aligned(band, ports[first, b] % 12), measure.aligned(band, ports[b, first] % 12), band
        rows = []
        for placement in layout['rows']:
            a, d = placement['owner'], placement['distance']
            owners = (a, b, c)
            geom = geometry(board, graph, ports, owners)
            modes = {}
            for k in ((int(candidate.KAPPA), 0, -int(candidate.KAPPA)) if lev == 0 else (int(candidate.KAPPA),)):
                modes[f'singlet_{k}'] = trial(graph, ports, owners, frames, geom, k, 'singlet_clamp', rng)
            if lev == 0:
                modes['relaxed_8'] = trial(graph, ports, owners, frames, geom, 8, 'relaxed', rng)
            rows.append(dict(distance=d, owners=owners, mobile_ports=len(geom[1]), trials=modes))
        data[str(lev)] = rows
    return currents._portable(dict(seed=SEED, shots=SHOTS, levels=data,
                                  slopes={level: slopes(rows) for level, rows in data.items()}))


def replay(write=False):
    path = measure.HERE / 'dressing.json'
    acquired = run()
    if write:
        path.write_text(json.dumps(acquired, indent=1, allow_nan=False) + '\n')
    elif acquired != json.loads(path.read_text()):
        raise ValueError('OW-22 dressing replay differs')
    print('OW-22 freely mobile dressing records verified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    replay(parser.parse_args().write)
