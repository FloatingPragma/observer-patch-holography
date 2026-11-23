"""Test autonomous endpoint selection on the native two-double-load regions."""

import argparse
from functools import lru_cache
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, diags

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.measure import angular_frame
from model.executable_model import federation
from model.experiments.interval_echo_space.shell_source import prepare
from ..equilibrium_refinement.measure import apparatus, seam
from ..hidden_entropy.measure import two_dark_modes
from ..shape_stationarity.geometry import selected_shapes
from ..shape_stationarity.measure import partition

HERE = Path(__file__).resolve().parent
SHOTS = 2**28
SEED = 23112377


@lru_cache(None)
def local_modes(gateway_addresses):
    local, frame, _, _ = angular_frame()
    target = two_dark_modes(frame, np.array(gateway_addresses))
    equipment = SimpleNamespace(car=local, active=local.seams)
    # A complete orthonormal launch menu allows transferring probability without
    # any postselection. Same-cell two-load outcomes use the sharp pair basis.
    basis = [column.copy() for column in target.T]
    for sharp in np.eye(local.n):
        for col in basis:
            sharp -= np.dot(col, sharp)*col
        if np.linalg.norm(sharp) > 1e-9:
            basis.append(sharp/np.linalg.norm(sharp))
    prepared = np.column_stack([prepare(equipment, mode, amount=2)[0] for mode in basis])
    np.testing.assert_allclose(prepared.conj().T @ prepared, np.eye(local.n), atol=1e-12)
    return prepared[:, :2]


def stationary_pair_tangent(car, cells, words):
    """Signed sharp-word direction and normalized product-wave columns for A and B."""
    three = sorted(cells)[:3]
    waves = {}
    for i in three:
        gates = tuple(np.flatnonzero(car.degree[12*i:12*i+12] == 6))
        waves[i] = local_modes(gates)
    pair = np.asarray(words)
    products = []
    for later in three[1:]:
        mask = (pair[:, 0]//12 == three[0]) & (pair[:, 1]//12 == later)
        cols = np.zeros((len(words), 4), dtype=complex)
        left = waves[three[0]][pair[mask,0] % 12]
        right = waves[later][pair[mask,1] % 12]
        cols[mask] = (left[:, :, None] * right[:, None, :]).reshape(-1, 4)
        np.testing.assert_allclose(cols.conj().T @ cols, np.eye(4), atol=1e-12)
        products.append(cols)
    a, b = products
    return ((abs(a)**2 - abs(b)**2).sum(axis=1)/4, a, b)


def pair_hamiltonian(car, contacts, words):
    """Enumerate the actual occupied-vacant exchange on every accessible two-load word."""
    block = seam()
    np.testing.assert_array_equal(np.diag(block), [0, 32.5, 32.5, 0])
    neighbors = {p: set() for w in words for p in w}
    for a, b in contacts:
        neighbors[a].add(b)
        neighbors[b].add(a)
    index = {w: i for i, w in enumerate(words)}
    rows, cols, energy = [], [], []
    for i, (a,b) in enumerate(words):
        for src, resident in ((a,b), (b,a)):
            for dest in neighbors[src] - {resident}:
                rows.append(index[tuple(sorted([dest, resident]))])
                cols.append(i)
        active = len(neighbors[a]) + len(neighbors[b]) - 2*(a in neighbors[b])
        energy.append(block[1,1]*active + 4*candidate.KAPPA*(12-int(car.degree[a])-int(car.degree[b])))
    hops = coo_matrix((np.full(len(rows), block[1,2]), (rows, cols)), shape=(len(words),len(words)))
    return (hops.tocsr() + diags(energy)).tocsr()


def tv(p,q):
    return float(np.sum(abs(p-q))/2)


def stationary_experiment(car, cells, words, h, inside, cuts, cut_words, seed):
    size = len(words)
    direction, a, b = stationary_pair_tangent(car, cells, words)
    energy = float(np.vdot(a[:,0], h @ a[:,0]).real)
    eigen_errors = [float(np.linalg.norm(h @ x - energy*x, axis=0).max()) for x in (a,b)]
    if max(eigen_errors) > 1e-11:
        raise ArithmeticError("Enrolled modes do not have one common pair energy")
    cut_grad = np.bincount(cuts, weights=direction)
    if max(abs(cut_grad)) > 2e-14:
        raise ArithmeticError("Transferred matter changes the exterior word")
    bulk_grad = np.bincount(inside, weights=direction)
    bulk_ref = np.bincount(inside)/size
    samples = []
    for j, alpha in enumerate((-.75,0.,.75)):
        word_prob = (1 + 4*alpha*direction)/size
        bulk = np.bincount(inside, weights=word_prob)
        cut = np.bincount(cuts, weights=word_prob)
        rng = np.random.default_rng([SEED,seed,0,j])
        b_counts = rng.multinomial(SHOTS, bulk/bulk.sum())
        e_counts = rng.multinomial(SHOTS, cut/cut.sum())
        b_read = b_counts/SHOTS
        slope = -float(bulk_grad @ np.log(bulk))
        measured = -float(bulk_grad @ np.log(b_read))
        se = float(np.sqrt(np.sum(bulk_grad**2 / b_read)/SHOTS))
        samples.append(dict(alpha=alpha, bulk_counts=b_counts.tolist(), cut_counts=e_counts.tolist(),
                            independent_shots_per_meter=SHOTS, signed_bulk_slope_nats=slope,
                            measured_slope_nats=measured, slope_standard_error=se,
                            bulk_tv_from_uniform_words=tv(bulk,bulk_ref)))
    return dict(common_energy=energy, eigenvector_residuals=eigen_errors,
                full_energy_distribution_equal_by_eigenspace_trace=True,
                certified_cut_word_tv=float(np.abs(cut_grad).sum()/2),
                preparations=samples)


def launch_pairs(energy):
    target = energy.mean()
    low = energy[energy < target].max()
    high = energy[energy > target].min()
    weight = (high-target)/(high-low)
    indexes = np.array([[np.flatnonzero(energy == e)[j] for e in (low,high)] for j in (0,-1)])
    if np.unique(indexes).size != 4:
        raise ArithmeticError("The equal-mean preparations require distinct inputs")
    return indexes, np.array([weight,1-weight])


def diagonal_record_laws(values, v, launches, tolerance):
    """Group unresolved energies before squaring, retaining all within-block coherence."""
    groups = np.split(np.arange(len(values)), np.flatnonzero(np.diff(values) > tolerance)+1)
    rows = v[launches]
    squared = v*v
    result = squared @ (rows*rows).T
    for idx in groups:
        if len(idx) > 1:
            block = v[:,idx] @ rows[:,idx].T
            result += block*block - squared[:,idx] @ (rows[:,idx]**2).T
    np.testing.assert_allclose(result.sum(axis=0), 1, atol=5e-12)
    return result, len(groups)


def free_records(h, inside, cut_ids, seed):
    diagonal = h.diagonal()
    indices, weights = launch_pairs(diagonal)
    energy, v = eigh(h.toarray(), driver='evd', check_finite=False)
    centered = energy-energy[0]
    # These are independently prepared sharp endpoint branches of a classical two-word mixture.
    launches = indices.ravel()
    def endpoints(t):
        contrib = np.cos(t*centered) - 1j*np.sin(t*centered)
        waves = v @ (contrib[:, None]*v[launches].T)
        return abs(waves)**2
    def pool(raw):
        return raw.reshape(len(energy),2,2) @ weights
    bulk_ref = np.bincount(inside)/len(inside)
    exterior_ref = np.bincount(cut_ids)/len(cut_ids)
    def read(raw, rng):
        result = []
        for p in pool(raw).T:
            bulk = np.bincount(inside, weights=p)
            cut = np.bincount(cut_ids, weights=p)
            bulk = np.maximum(bulk,0); bulk /= bulk.sum()
            cut = np.maximum(cut,0); cut /= cut.sum()
            b_counts = rng.multinomial(SHOTS,bulk)
            e_counts = rng.multinomial(SHOTS,cut)
            # Allow independent replay of the observation without storing every sparse high-shot bin.
            result.append(dict(bulk_tv=tv(bulk,bulk_ref), exterior_tv=tv(cut,exterior_ref),
                               sampled_bulk_tv=tv(b_counts/SHOTS,bulk_ref),
                               sampled_exterior_tv=tv(e_counts/SHOTS,exterior_ref),
                               coarse_interior_counts=b_counts.tolist(), cut_counts=e_counts.tolist()))
        return result
    finite = []
    for j, time in enumerate((1,4,16,64)):
        finite.append(dict(exposure=time, preparations=read(endpoints(time),
                            np.random.default_rng([SEED,seed,1,j]))))
    # Pool single terminal shots taken at uniformly randomized appointments.
    appointments = np.linspace(64.,128.,17)
    window = sum(endpoints(t) for t in appointments)/len(appointments)
    timed = read(window, np.random.default_rng([SEED,seed,2]))
    fine, n = diagonal_record_laws(energy, v, launches, 1e-10)
    coarse, n_coarse = diagonal_record_laws(energy,v,launches,1e-8)
    delta = pool(fine)-pool(coarse)
    sensitivity = float(np.max(np.abs(delta).sum(axis=0))/2)
    dephased = read(fine, np.random.default_rng([SEED,seed,3]))
    # The lowest levels can be unresolved at double precision. Compute a word
    # event common to every vector in the isolated cluster and bound its TV.
    cluster = np.flatnonzero(centered < 1e-8)
    reference = v[:,cluster]
    word_prob = np.sum(reference**2,axis=1)/len(cluster)
    low_bulk = np.bincount(inside,weights=word_prob)
    event = low_bulk > bulk_ref
    restricted = reference[event[inside]]
    low_event = np.linalg.eigvalsh(restricted.T @ restricted)
    ground_tv_floor = max(0., float(low_event[0] - bulk_ref[event].sum()))
    cluster_residual = float(np.max(np.linalg.norm(h @ reference - reference*energy[cluster],axis=0)))
    return dict(sharp_launch_indices=indices.tolist(), launch_weights=weights.tolist(),
                common_native_mean_energy=float(diagonal.mean()), fixed_sector_ground_energy=float(energy[0]),
                ground_eigenvector_residual=float(np.linalg.norm(h @ v[:,0] - energy[0]*v[:,0])),
                ground_gap=float(energy[1]-energy[0]), lowest_cluster_dimension=len(cluster),
                gap_to_next_cluster=float(centered[len(cluster)]-centered[cluster[-1]]),
                maximum_cluster_eigenvector_residual=cluster_residual,
                lowest_cluster_bulk_tv=float(tv(low_bulk,bulk_ref)),
                uniform_tv_lower_bound_for_any_lowest_cluster_state=ground_tv_floor,
                finite_appointments=finite, uniform_appointment_grid=appointments.tolist(),
                appointment_average=timed, spectral_average=dephased,
                spectral_groups_at_1e10_tolerance=n, spectral_groups_at_1e8_tolerance=n_coarse,
                spectral_tolerance_word_tv_sensitivity=sensitivity)


def run():
    rows = []
    for level in range(3):
        car = federation(level)
        for j,(name, cells) in enumerate(selected_shapes(car)):
            _, contacts, cut, words = apparatus(car, cells)
            h = pair_hamiltonian(car,contacts,words)
            inside, cuts, _, inside_words, cut_words = partition(words,cut)
            seed = 4*level+j
            rows.append(dict(board_cells=car.cells, shape=name, cells=sorted(cells),
                             pair_words=len(words), interior_words=inside_words, cut_words=cut_words,
                             stationary=stationary_experiment(car,cells,words,h,inside,cuts,cut_words,seed),
                             flight=free_records(h,inside,cuts,seed)))
            print(level, name, rows[-1]['stationary']['preparations'][-1]['signed_bulk_slope_nats'],
                  [p['bulk_tv'] for p in rows[-1]['flight']['spectral_average']], flush=True)
    return dict(seed=SEED, shots_per_meter=SHOTS, rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    dest = parser.parse_args().output
    result = run()
    with dest.open('x') as f:
        json.dump(result, f, separators=(',',':'), allow_nan=False)
        f.write('\n')
