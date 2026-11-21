"""Free flight with a retained seam pointer on the unmodified spherical federation."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.executable_model.carrier import seam_digest


HERE = Path(__file__).resolve().parent
LAW = candidate.native_law
KAPPA = candidate.KAPPA
TIMES = (0, .25, .5, 1, 2, 4, 8)
SHOTS = 16384


def local_indices(labels, occupancies):
    return [9 * (1 + LAW.LABELS.index(label)) + LAW.SRC.index(pair)
            for label in labels for pair in occupancies]


def exact_collar():
    """Eight times the nonblank Hamiltonian in the registered integer write basis."""
    write = np.rint(np.sqrt(2) * LAW.write()).astype(np.int64)
    g0, g1 = LAW.generator_parts()
    g2 = np.rint(2 * (g0 + KAPPA * g1)[9:, 9:]).astype(np.int64)
    image2 = write @ write.T
    complement2 = 2 * np.eye(90, dtype=np.int64) - image2
    collar8 = complement2 @ g2 @ complement2
    assert not np.any(g2 @ image2 - image2 @ g2)
    assert not np.any(complement2 @ g2 @ write)  # blank / marked coupling
    h0, h1 = LAW.local_term_parts()
    native = h0 + KAPPA * h1
    assert np.max(abs(8 * native[9:, 9:] - collar8)) < 1e-11
    assert np.max(abs(native[:9, 9:])) < 1e-12

    abs_contrast = np.repeat([abs(d) for _, d in LAW.LABELS], 9)
    delta = abs_contrast[:, None] - abs_contrast
    assert not np.any(collar8 * delta)
    retained_double = [i - 9 for i in local_indices(LAW.LABELS, [(0, 0), (2, 0), (0, 2)])]
    diff = collar8 - 4 * g2
    assert not np.any(diff[:, retained_double])
    retained_single = [i - 9 for i in local_indices(LAW.LABELS, [(1, 0), (0, 1)])]
    single_residual = max(abs(diff[:, retained_single]).flat)
    assert single_residual > 0
    return native, dict(native_states=99, marked_double_or_vacuum_columns=len(retained_double),
                        double_integer_residual=0, marked_single_columns=len(retained_single),
                        single_max_entry_difference=float(single_residual / 8),
                        fixed_absolute_contrast=True, fixed_blank_status=True)


def uniform_background(native):
    """Recover the coin matrix units from the two directed single-unit transfers."""
    labels = [(c, d) for c in (0, 1) for d in (-1, 1)]
    a = local_indices(labels, [(1, 0)])
    b = local_indices(labels, [(0, 1)])
    # The collar test certifies eight times every marked entry by integer arithmetic.
    eight_h = np.rint(8 * native).astype(np.int64)
    out = []
    for parity in (1, -1):
        embed = np.kron(np.eye(2, dtype=int), [[1], [parity]])
        right = eight_h[np.ix_(b, a)]
        left = eight_h[np.ix_(a, b)]
        g = int(KAPPA) + (1 - parity) // 2
        unit = np.array([[0, 1], [0, 0]])
        # Using unnormalized sign eigenvectors avoids square roots.
        assert np.array_equal(right @ embed, -4 * g * embed @ unit)
        assert np.array_equal(left @ embed, -4 * g * embed @ unit.T)
        units = [unit, unit.T, unit @ unit.T, unit.T @ unit]
        products = np.array([np.kron(x, y).ravel() for x in units for y in units])
        assert np.array_equal(products @ products.T, np.eye(16, dtype=int))
        out.append(dict(sign_parity=parity, forward_coin_0_from_1=-g / 2,
                        reverse_coin_1_from_0=-g / 2, one_coin_algebra_dimension=4,
                        two_coin_algebra_dimension=16))
    vac = local_indices(labels, [(0, 0)])
    endpoints = [(1, 0), (0, 1)]
    occ = local_indices(labels, endpoints)
    blank = [LAW.SRC.index(pair) for pair in endpoints]
    difference = (eight_h[np.ix_(occ, occ)] - np.kron(eight_h[np.ix_(vac, vac)], np.eye(2, dtype=int))
                  - np.kron(np.eye(4, dtype=int), eight_h[np.ix_(blank, blank)]))
    one_edge_rank = int(np.linalg.matrix_rank(difference))
    return dict(parities=out, one_retained_edge_endpoint_dimension=8,
                one_retained_edge_numerical_perturbation_rank=one_edge_rank,
                uniform_edges=[len(federation(level).seams) for level in (0, 1, 2)],
                lower_bound_common_internal_dimension="2^number_of_marked_contrast_one_seams")


def generator(car, native, load, labels, edge, coupled=True, marked=True):
    """Sum all free terms and replace one blank input block by its marked block."""
    endpoints = [(load, 0), (0, load)]
    blank_indices = [LAW.SRC.index(pair) for pair in endpoints]
    blank = native[np.ix_(blank_indices, blank_indices)]
    lift = sparse.csr_matrix((np.ones(2), (edge, [0, 1])), shape=(car.n, 2))
    # All blank seams and degree-reference terms, with a constant 6*kappa*n^2 removed.
    full = car.laplacian() * .5
    assert np.allclose(blank, [[KAPPA * load**2 + .5, -.5], [-.5, KAPPA * load**2 + .5]])
    if not marked:
        return full
    blank_except_gateway = full - lift @ sparse.csr_matrix(blank) @ lift.T
    vac = local_indices(labels, [(0, 0)])
    occ = local_indices(labels, endpoints)
    vacuum = native[np.ix_(vac, vac)]
    occupied = native[np.ix_(occ, occ)]
    for block in (vac, occ):
        omitted = np.setdiff1d(np.arange(LAW.DIM), block)
        assert np.max(abs(native[np.ix_(omitted, block)])) < 1e-11
    count = len(labels)
    h = sparse.kron(sparse.eye(count), blank_except_gateway) + sparse.kron(vacuum, sparse.eye(car.n))
    if coupled:
        on_edge = sparse.kron(sparse.eye(count), lift, format='csr')
        interaction = sparse.csr_matrix(occupied - np.kron(vacuum, np.eye(2)))
        h = h + on_edge @ interaction @ on_edge.T
    return h.tocsr()


def native_receiver(load):
    """Contrast read after the echo on a blank receiver gives its occupancy union."""
    inputs = [LAW.SRC.index(pair) for pair in [(0, 0), (load, 0), (0, load)]]
    selected = local_indices([l for l in LAW.LABELS if abs(l[1]) == load], LAW.SRC)
    columns = LAW.echo()[np.ix_(selected, inputs)]
    assert np.max(abs(columns.T @ columns - np.diag([0, 1, 1]))) < 1e-12


def fly(operator, initial):
    states = expm_multiply(-1j * operator, initial, start=0, stop=8, num=33)
    assert np.max(abs(np.sum(abs(states)**2, axis=1) - 1)) < 5e-11
    return states[[round(4 * t) for t in TIMES]]


def geometry(car):
    a, b = car.external[0]
    other_a = next(p for p in car.neighbours[a] if p // 12 == a // 12)
    other_b = next(p for p in car.neighbours[b] if p // 12 == b // 12)
    adjacency = car.adjacency().tolil()
    adjacency[a, b] = adjacency[b, a] = 0
    dist = shortest_path(adjacency.tocsr(), directed=False, indices=[a, b], unweighted=True)
    assert np.isfinite(dist).all()
    # No marked link is used to transport either preparation to the remote port.
    far = int(np.argmax(dist.min(axis=0)))
    return (a, b), (b, other_b), other_a, far, [int(d[far]) for d in dist]


def single(car, native, edge, receiver, far, rng):
    labels = [(coin, sign) for coin in (0, 1) for sign in (-1, 1)]
    prepared = LAW.echo()[:, LAW.SRC.index((1, 0))]
    for coin, endpoint in [(0, (0, 1)), (1, (1, 0))]:
        ix = local_indices([(coin, 1)], [endpoint])[0]
        assert abs(abs(prepared[ix])**2 - .5) < 1e-12
    native_receiver(1)
    out = {}
    for name in ['near', 'far', 'decoupled']:
        h = generator(car, native, 1, labels, edge, coupled=name != 'decoupled')
        experiments = []
        for coin in (0, 1):
            source_port = far if name == 'far' else edge[1 - coin]
            initial = np.zeros((len(labels), car.n))
            initial[labels.index((coin, 1)), source_port] = 1
            mask_coin = np.array([label[0] != coin for label in labels])[:, None]
            mask_receiver = np.isin(np.arange(car.n), receiver)
            outcome = (2 * mask_coin + mask_receiver).ravel().astype(int)
            states = fly(h, initial.ravel())
            probabilities = np.stack([np.bincount(outcome, abs(state)**2, minlength=4) for state in states])
            experiments.append(probabilities)
        entries = []
        for step, t in enumerate(TIMES):
            n0 = int(rng.binomial(SHOTS, .5))
            arms = []
            for coin, n in enumerate((n0, SHOTS - n0)):
                p = experiments[coin][step]
                clicks = rng.multinomial(n, p / p.sum())
                arms.append(dict(initial_coin=coin, joint_probabilities=np.round(p, 10).tolist(),
                                 preparation_count=n, joint_counts=clicks.tolist()))
            entries.append(dict(exposure=t, arms=arms))
        out[name] = entries
    def largest_flip(name):
        return max(sum(a['joint_probabilities'][2:]) for row in out[name] for a in row['arms'])
    assert largest_flip('near') > .01 and largest_flip('decoupled') < 1e-6
    return dict(maximum_coin_flip={name: round(largest_flip(name), 10) for name in out}, settings=out)


def double(car, native, edge, receiver, source_port, rng):
    labels = [(0, 0), (1, 0)]
    native_receiver(2)
    vacuum_preparation = LAW.echo()[:, LAW.SRC.index((0, 0))]
    coin_indices = local_indices(labels, [(0, 0)])
    assert np.max(abs(vacuum_preparation[coin_indices] - 1 / np.sqrt(2))) < 1e-12
    cases = {}
    for marked in (False, True):
        h = generator(car, native, 2, labels, edge, marked=marked)
        initial = np.zeros((2 if marked else 1, car.n))
        initial[:, source_port] = vacuum_preparation[coin_indices] if marked else 1
        states = fly(h, initial.ravel()).reshape(len(TIMES), -1, car.n)
        probability = (abs(states[:, :, receiver])**2).sum(axis=(1, 2))
        cases['vacuum_write' if marked else 'blank'] = [
            dict(exposure=t, receiver_probability=round(float(p), 10),
                 receiver_clicks=int(rng.binomial(SHOTS, p))) for t,p in zip(TIMES, probability, strict=True)]
    return cases


def run():
    native, audit = exact_collar()
    boards = []
    for level in [0, 1]:
        car = federation(level)
        edge, receiver, near, far, distance = geometry(car)
        rng = np.random.default_rng(304050 + level)
        boards.append(dict(cells=car.cells, ports=car.n, seams=len(car.seams),
                           carrier_sha256=seam_digest(car), marked_seam=edge, receiver_seam=receiver,
                           double_source_port=near, far_single_port=far, far_distances=distance,
                           single=single(car, native, edge, receiver, far, rng),
                           double=double(car, native, edge, receiver, near, rng)))
    return dict(operator=audit, exposure_coordinate='native service exposure',
                joint_order=['coin held, quiet', 'coin held, receiver', 'coin changed, quiet', 'coin changed, receiver'],
                trials_per_setting=SHOTS, boards=boards)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--output', type=Path, help='write the free-propagation receipt to a fresh file')
    group.add_argument('--algebra-output', type=Path, help='write the uniform-background certificate to a fresh file')
    args = parser.parse_args()
    algebra = uniform_background(exact_collar()[0])
    receipt = json.loads(json.dumps(run()))
    if args.output or args.algebra_output:
        path, data = (args.output, receipt) if args.output else (args.algebra_output, algebra)
        with path.open('x') as f:
            json.dump(data, f, indent=2)
            f.write('\n')
    else:
        assert receipt == json.loads((HERE / 'results.json').read_text()), 'saved receipt differs'
        assert algebra == json.loads((HERE / 'algebra.json').read_text()), 'saved algebra certificate differs'
    print(receipt['operator'])
    print(algebra)
    for b in receipt['boards']:
        print(b['cells'], b['single']['maximum_coin_flip'])
