"""Universal comparison bounds for finite internal transports on the written cell."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import breadth_first_order
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[6]
sys.path.insert(0, str(REPO))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.family_projector.measure import full_spectrum_bound
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.transfer.measure import extract_sources


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def spin_generators():
    # Exact realification of the two complex spinor components, in order
    # (Re z0, Re z1, Im z0, Im z1). Conjugate transpose becomes transpose.
    x = np.array([[0, 1], [1, 0]], int)
    z = np.diag([1, -1])
    zero = np.zeros((2, 2), int)
    gx, gz = (np.block([[zero, -p], [p, zero]]) for p in (x, z))
    identity = np.eye(4, dtype=int)
    for g in (gx, gz):
        np.testing.assert_array_equal(g.T @ g, identity)
        np.testing.assert_array_equal(g @ g, -identity)
    np.testing.assert_array_equal(gx@gz@gx.T@gz.T, -identity)
    return identity, gx, gz


def lifted_matrix(diagonal, row, col, weights, blocks):
    # Every holonomy used for the numerical control is a signed permutation
    # of four real coordinates, so each block needs only four sparse entries.
    count, rank, _ = blocks.shape
    pick = np.argmax(abs(blocks), axis=2)
    values = np.take_along_axis(blocks, pick[:, :, None], axis=2).reshape(count, rank)
    np.testing.assert_array_equal(np.abs(values), 1)
    source = row[:, None]*rank+np.arange(rank)[None, :]
    dest = col[:, None]*rank+pick
    nodes = len(diagonal)*rank
    h = sparse.csr_matrix(((weights[:, None]*values).ravel(),
                           (source.ravel(), dest.ravel())), shape=(nodes, nodes))
    return h + h.T + sparse.diags(np.repeat(diagonal, rank), dtype=np.int64, format='csr')


def cycle_sites(u, v, predecessor):
    def to_root(s):
        trail = [int(s)]
        while predecessor[trail[-1]] >= 0:
            trail.append(int(predecessor[trail[-1]]))
        return trail
    l, r = to_root(u), to_root(v)
    r_index = {p: i for i, p in enumerate(r)}
    lca = next(i for i in l if i in r_index)
    return l[:l.index(lca)+1] + r[:r_index[lca]][::-1] + [int(u)]


def extremum(native, which):
    start = np.random.default_rng(31817).normal(size=native.shape[0])
    energies, vectors = eigsh(native, k=1, which=which, tol=1e-11, v0=start)
    energy, v = float(energies[0]), vectors[:, 0]
    residual = float(np.linalg.norm(native @ v - energy*v))
    assert residual < 5e-9
    return energy, residual


def test_connections(doubled, states, positions, scalar_ground):
    count = len(states)
    ident, gx, gz = spin_generators()
    off = sparse.triu(doubled, k=1, format='coo')
    row, col, weights = off.row, off.col, off.data
    matter = positions[row] != positions[col]
    np.testing.assert_array_equal(weights, np.where(matter, -8, -1))
    graph = doubled.copy()
    graph.setdiag(0)
    graph.eliminate_zeros()
    _, predecessor = breadth_first_order(graph, i_start=0, directed=False)
    tree = (predecessor[col] == row) | (predecessor[row] == col)
    assert sum(tree) == count-1
    # Exact configuration-dependent local gauges: U_ij = G_i G_j inverse.
    frames = np.stack([ident, gx, gz, gx@gz])[np.arange(count)%4]
    transformed = lifted_matrix(doubled.diagonal(), row, col, weights,
                                frames[row]@frames[col].transpose(0, 2, 1))
    gauge = sparse.bsr_matrix((frames, np.arange(count), np.arange(count+1)), shape=(count*4, count*4))
    ordinary = sparse.kron(doubled, sparse.eye(4, dtype=np.int64), format='csr')
    conjugacy = transformed - gauge @ ordinary @ gauge.T
    conjugacy.eliminate_zeros()
    assert conjugacy.nnz == 0
    # A tree path transports the identity. Each non-tree edge supplies a
    # fundamental cycle based at its vertex, referred to the root along the tree.
    kinds = np.where(tree, 0, np.where(matter, 1, 2))
    blocks = np.stack([ident, gx, gz])[kinds]
    curved = lifted_matrix(doubled.diagonal(), row, col, weights, blocks).astype(float)/2
    examples = []
    for kind, name in ((1, 'material'), (2, 'ring')):
        edge = int(np.flatnonzero(kinds == kind)[0])
        sites = cycle_sites(row[edge], col[edge], predecessor)
        assert sites[0] == sites[-1] and len(set(sites)) == len(sites)-1
        for a, b in zip(sites[:-2], sites[1:-1]):
            assert predecessor[a] == b or predecessor[b] == a
        assert doubled[sites[-2], sites[-1]] == weights[edge]
        examples.append({'closing_contact': name, 'chords_of_this_kind': int(np.sum(kinds == kind)),
                         'fundamental_cycle_state_ids': sites,
                         'cycle_charge_ports': [int(states[i][0]) for i in sites],
                         'cycle_head_words': [int(states[i][1]) for i in sites]})
    scalar_max, residual = extremum(doubled.astype(float)/2, 'LA')
    low, low_err = extremum(curved, 'SA')
    high, high_err = extremum(curved, 'LA')
    assert low > scalar_ground and high <= np.asarray(abs(doubled).sum(axis=1)).max()/2+1e-10
    return {'vertices': count, 'edges': len(weights), 'matter_edges': int(sum(matter)),
            'ring_edges': int(sum(~matter)), 'tree_edges': int(sum(tree)),
            'pure_gauge_full_matrix_conjugacy_entry_errors': conjugacy.nnz,
            'pure_gauge_spectrum_complex_multiplicity_factor': 2,
            'pure_gauge_ground_and_ceiling': [scalar_ground, scalar_max],
            'scalar_ceiling_eigenpair_residual': residual, 'curved_real_matrix_dimension': curved.shape[0],
            'curved_examples': examples, 'curved_rooted_holonomy_commutator': 'minus identity',
            'curved_extrema': [low, high], 'curved_eigenpair_residuals': [low_err, high_err]}


def neutral_bounds(doubled, states, positions, port, certificate):
    rows = np.flatnonzero(positions == port)
    ring = doubled[rows][:, rows].copy()
    ring.setdiag(0)
    ring.eliminate_zeros()
    np.testing.assert_array_equal(ring.data, -1)
    native_diagonal = -np.asarray(ring.sum(axis=1)).ravel()/2
    i = int(np.argmin(native_diagonal))
    vmax = float(native_diagonal[i])
    numerator, denominator = certificate['strict_lower_native_fraction']
    lower, upper = numerator/denominator, certificate['upper_native_gershgorin']
    assert 0 <= vmax < lower
    ratio = (upper-vmax)/(lower-vmax)
    r = np.sqrt(ratio)
    assert ratio < 64
    return {'neutral_sharp_configurations': len(rows), 'charged_reference_port': int(port),
            'neutral_trial_head_word': int(states[rows[i]][1]), 'flippable_faces_at_trial': int(2*vmax),
            'all_connections_neutral_ground_range': [0, vmax],
            'charged_gap_strict_lower_at_largest_vacuum': lower-vmax,
            'charged_gap_upper_at_largest_vacuum': upper-vmax, 'all_gap_ratio_upper': float(ratio),
            'all_gap_koide_upper': float((r*r+2)/(r+2)**2)}


def run():
    prior_folder = HERE.parent/'family_projector'
    prior = json.loads((prior_folder/'results.json').read_bytes())
    assert prior['executable_sha256'] == digest(prior_folder/'measure.py')
    verified_keys = ['upper_native_gershgorin', 'strict_lower_native_fraction',
                     'min_positive_ground_witness_entry', 'native_ground_witness_scale',
                     'ground_witness_int64_le_sha256', 'minimum_integer_margin',
                     'max_to_min_ratio_bound', 'all_triplets_koide_upper_bound']
    sectors = []
    with TemporaryDirectory() as tmp:
        hashes, _ = extract_sources(Path(tmp))
        assert hashes == prior['archived_dependency_sha256']
        from search import Patch, charged_sector, lowest
        from bridge import matched_bridge
        patch = Patch()
        gateway = matched_bridge()[0]
        for s, expected in zip((-1, 1), prior['sectors'], strict=True):
            assert s == expected['charge_sign']
            doubled, states, positions = charged_sector(patch, s)
            levels, modes, _, _ = lowest(doubled)
            checked = full_spectrum_bound(doubled, modes[:, 0])
            for key in verified_keys:
                assert checked[key] == expected['whole_spectrum_bound'][key], key
            sectors.append({'charge_sign': s, 'verified_scalar_keys': {k: checked[k] for k in verified_keys},
                            'unitary_controls': test_connections(doubled, states, positions, float(levels[0])),
                            'matched_neutral': neutral_bounds(doubled, states, positions, gateway, checked)})
    # Separate neutral subtractions for opposite charge backgrounds permit
    # a wider common envelope when pooling levels from both signs.
    lowest_gap = min(s['matched_neutral']['charged_gap_strict_lower_at_largest_vacuum'] for s in sectors)
    largest_gap = max(s['verified_scalar_keys']['upper_native_gershgorin'] for s in sectors)
    r = np.sqrt(largest_gap/lowest_gap)
    pooling = {'strict_lower_gap': lowest_gap, 'upper_gap': largest_gap,
               'ratio_upper': largest_gap/lowest_gap, 'koide_upper': float((r*r+2)/(r+2)**2)}
    return {'design_sha256': digest(HERE/'DESIGN.md'), 'executable_sha256': digest(Path(__file__)),
            'scalar_certificate_sha256': digest(prior_folder/'results.json'), 'source_digests': hashes,
            'scope': 'every finite-dimensional unitary per transition with fixed positive weights and scalar port costs',
            'opposite_sign_pooling_with_independent_neutral_subtractions': pooling, 'sectors': sectors}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'results.json')
    parser.parse_args().output.write_text(json.dumps(run(), indent=2, allow_nan=False)+'\n')
