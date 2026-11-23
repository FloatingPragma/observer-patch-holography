"""Attach the exact A5 port projectors to the accessible written charge field."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from math import lcm, sqrt
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory

import numpy as np
from sympy import Matrix

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[6]
sys.path.insert(0, str(REPO))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.transfer.measure import extract_sources
from model.executable_model.carrier import CELL_SEAMS, cell, isomorphism


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def theorem_projectors():
    lean = REPO / 'Lean/Screen/A5FamilyBand.lean'
    content = lean.read_text()

    def table(name):
        match = re.search(r'def '+name+r' : List \(List Int\) := (\[.*?\]\])', content, re.DOTALL)
        if not match:
            raise ValueError(f'Cannot read coefficient table {name}')
        data = np.asarray(ast.literal_eval(match[1]), int)
        assert data.shape == (12, 12)
        return data

    a, x, y = (table(name) for name in ('At', 'Xt', 'Yt'))
    theorem_seams = [tuple(edge) for edge in np.argwhere(np.triu(a, 1))]
    mapping = isomorphism(theorem_seams, CELL_SEAMS)
    assert mapping is not None and len(set(mapping)) == 12
    order = np.argsort(mapping)
    a, x, y = (m[np.ix_(order, order)] for m in (a, x, y))
    np.testing.assert_array_equal(a, cell().adjacency().toarray())
    np.testing.assert_array_equal(a@x, 5*y)
    np.testing.assert_array_equal(a@y, x)
    np.testing.assert_array_equal(x@x+5*y@y, 40*x)
    np.testing.assert_array_equal(x@y+y@x, 40*y)
    antipode = np.eye(12, dtype=int) - x//10
    np.testing.assert_array_equal(antipode.sum(axis=1), 1)
    np.testing.assert_array_equal(antipode@x, -x)
    np.testing.assert_array_equal(antipode@y, -y)
    result = {}
    for sign in (1, -1):
        project = (x + sign*np.sqrt(5)*y)/40
        eigenvalues, vec = np.linalg.eigh(project)
        basis = vec[:, eigenvalues > .5]
        assert basis.shape == (12, 3)
        np.testing.assert_allclose(basis@basis.T, project, atol=2e-15, rtol=0)
        # Physical scalar positive Laplacian on the twelve-port cell.
        h = 4*(5*np.eye(12)-a)
        e = 4*(5-sign*np.sqrt(5))
        remainder = h @ basis - e*basis
        np.testing.assert_allclose(remainder, 0, atol=4e-14)
        result[sign] = (project, basis, e, float(np.linalg.norm(remainder, 2)))
    return result, {'lean_file': str(lean.relative_to(REPO)), 'lean_sha256': digest(lean),
                    'lean_to_physical_port': list(map(int, mapping)), 'projector_rank': 3, 'physical_antipodes': np.argmax(antipode, axis=1).tolist(),
                    'scalar_control_hopping': 4,
                    'scalar_control': [{'adjacency_sign': s, 'energy': data[2], 'max_leakage_bound': data[3]}
                                       for s, data in result.items()]}


def exact_head_contrast(doubled, ports):
    # Work in unnormalized sharp-word coordinates. Contrasting rows at the
    # same charge port cancels the entirely local diagonal potential.
    indicators = np.eye(12, dtype=np.int64)[ports]
    advanced = doubled @ indicators
    anchors = np.array([np.flatnonzero(ports == p)[0] for p in range(12)])
    contrast = advanced - advanced[anchors[ports]]
    np.testing.assert_array_equal(contrast.sum(axis=1), 0)
    # Select and then certify a nonzero exact minor. The constant vector is
    # already in the kernel. Check all extra null directions exactly below.
    prime = 65521
    work = contrast.copy() % prime
    row_ids = np.arange(len(work))
    pivots = []
    for c in range(12):
        r = len(pivots)
        choices = np.flatnonzero(work[r:, c])
        if not len(choices):
            continue
        i = r + choices[0]
        work[[r, i]] = work[[i, r]]
        row_ids[[r, i]] = row_ids[[i, r]]
        work[r] = work[r]*pow(int(work[r, c]), -1, prime) % prime
        work[r+1:] = (work[r+1:] - work[r+1:, c, None]*work[r]) % prime
        pivots.append(c)
    rows = row_ids[:len(pivots)]
    exact_det = int(Matrix(contrast[np.ix_(rows, pivots)].tolist()).det())
    assert exact_det % prime
    ker = Matrix(contrast[rows].tolist()).nullspace()
    kernel = []
    for v in ker:
        scale = lcm(*(int(x.q) for x in v))
        w = np.array([int(x * scale) for x in v])
        np.testing.assert_array_equal(contrast @ w, 0)
        kernel.append(w.tolist())
    assert kernel in ([[1]*12], [[1]+[0]*11, [0]+[1]*11])
    return {'same_port_reference_state_indices': anchors.tolist(), 'independent_contrast_rows': rows.tolist(),
            'independent_address_columns': pivots, 'nonzero_integer_minor_determinant': exact_det,
            'real_contrast_rank': len(pivots), 'integer_kernel_in_sharp_word_coordinates': kernel}


def full_spectrum_bound(doubled, ground):
    # The native Hamiltonian is doubled/2. Off-diagonals of doubled are
    # nonpositive; shifting its negative to a nonnegative matrix permits
    # Perron comparison with this explicitly positive integer witness.
    scale = 10**9
    x = np.rint(scale*ground*np.sign(ground[0])).astype(np.int64)
    assert np.min(x) > 0
    radius = int(np.asarray(abs(doubled).sum(axis=1)).max())
    # All products including the rational comparison fit the signed dtype.
    assert 1000*radius*int(x.max()) < np.iinfo(np.int64).max
    action = doubled @ x
    denominator = 1000
    numerator = int(np.min((denominator*action)//(2*x)))
    margin = denominator*action - 2*numerator*x
    assert numerator > 0 and int(np.min(margin)) > 0
    upper = radius/2
    lower = numerator/denominator
    ratio = upper/lower
    r = sqrt(ratio)
    kmax = (r*r+2)/(r+2)**2
    # An integer-only guard: ratio < 64 implies sqrt(ratio) < 8,
    # strictly inside the threshold 4 + 3*sqrt(2).
    assert radius*denominator < 128*numerator
    return {'upper_native_gershgorin': upper, 'strict_lower_native_fraction': [numerator, denominator],
            'min_positive_ground_witness_entry': int(x.min()), 'native_ground_witness_scale': scale,
            'ground_witness_int64_le_sha256': hashlib.sha256(x.astype('<i8').tobytes()).hexdigest(),
            'minimum_integer_margin': int(margin.min()), 'max_to_min_ratio_bound': ratio,
            'all_triplets_koide_upper_bound': kmax, 'koide_target': '2/3'}


def test_triplet(basis, w, h, low_states, eig_groups):
    u = w @ basis
    np.testing.assert_allclose(u.T@u, np.eye(3), atol=1e-13)
    evolved = h @ u
    port_component = w.T @ evolved
    inside = basis.T @ port_component
    port_leak = port_component-basis@inside
    head_leak = evolved-w@port_component
    full_leak = evolved-u@inside
    covariance = full_leak.T@full_leak
    np.testing.assert_allclose(covariance, port_leak.T@port_leak+head_leak.T@head_leak, atol=1e-12)
    svals = np.linalg.svd(full_leak, compute_uv=False)
    np.testing.assert_allclose(svals**2, np.linalg.eigvalsh(covariance)[::-1], atol=8e-13)
    overlap = u.T @ low_states
    group_weights = []
    for ids in eig_groups:
        block = overlap[:, ids]
        group_weights.append({'low_level_indices': [int(i) for i in ids],
                              'total_projector_weight': float(np.sum(block**2))})
    return {'leakage_singular_values': svals.tolist(), 'leakage_operator_norm': float(svals[0]),
            'smallest_initial_escape_probability_coefficient': float(svals[-1]**2),
            'projected_native_energies': np.linalg.eigvalsh(inside).tolist(),
            'leakage_into_other_uniform_port_modes_singular_values': np.linalg.svd(port_leak, compute_uv=False).tolist(),
            'leakage_into_orthogonal_head_patterns_singular_values': np.linalg.svd(head_leak, compute_uv=False).tolist(),
            'three_ground_level_weights': np.sum(overlap[:, :3]**2, axis=0).tolist(),
            'triplet_to_lowest_three_principal_singular_values': np.linalg.svd(overlap[:, :3], compute_uv=False).tolist(),
            'low_groups_projector_weights': group_weights}


def run():
    projectors, bound = theorem_projectors()
    out = []
    with TemporaryDirectory() as tmp:
        hashes, reference = extract_sources(Path(tmp))
        from search import Patch, charged_sector, lowest

        patch = Patch()
        for sign, archived in zip((-1, 1), reference['charge_sectors'], strict=True):
            doubled, states, ports = charged_sector(patch, sign)
            levels, modes, eig_groups, residual = lowest(doubled)
            assert archived['sign'] == sign
            h = doubled.astype(float)/2.
            counts = np.bincount(ports, minlength=12)
            assert np.all(counts > 0) and sum(counts) == h.shape[0]
            # A constant signed amplitude on all ports except possibly 0
            # cannot be odd on antipodes after positive word-count scaling.
            antipode = np.array(bound['physical_antipodes'])
            assert np.all(antipode != np.arange(12)) and antipode[0] != 0
            assert any(p != 0 and antipode[p] != 0 for p in range(12))
            w = np.zeros((len(states), 12))
            w[np.arange(len(states)), ports] = 1/np.sqrt(counts[ports])
            np.testing.assert_allclose(w.T@w, np.eye(12), atol=1e-13)
            assert all(len(group) == 1 for group in eig_groups[:3])
            out.append({'charge_sign': sign, 'sharp_configurations': len(states), 'reachable_word_counts_by_port': counts.tolist(),
                        'exact_head_contrast_certificate': exact_head_contrast(doubled, ports),
                        'whole_spectrum_bound': full_spectrum_bound(doubled, modes[:, 0]),
                        'native_low_energies': levels.tolist(), 'low_eigenpair_residual': residual,
                        'canonical_plus_sqrt5': test_triplet(projectors[1][1], w, h, modes, eig_groups),
                        'comparison_minus_sqrt5': test_triplet(projectors[-1][1], w, h, modes, eig_groups)})
    return {'design_sha256': digest(HERE/'DESIGN.md'), 'executable_sha256': digest(Path(__file__)),
            'archived_dependency_sha256': hashes, 'finite_port_theorem': bound,
            'embedding': 'normalized equal coherent sum of the connected head words at each charge port',
            'geometry': 'source cell of the 3–48 matched bridge on the 20-cell board',
            'native_generator': 'one half of the archived integer charged-sector matrix', 'sectors': out}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', '--out', dest='out', type=Path, default=HERE/'results.json')
    output = parser.parse_args().out
    output.write_text(json.dumps(run(), indent=2, allow_nan=False)+'\n')
