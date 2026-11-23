"""Measure the same field-head preparations with interior and perimeter observers."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.electron_field.flux_string.search import Patch, around


HERE = Path(__file__).resolve().parent
TIMES = (0, 0.5, 1, 2)
WEIGHTS = (0.625, 0.875)
SHOTS = 2**24


def ring_laplacian(patch, words):
    idx = {w: i for i, w in enumerate(words)}
    i, j = zip(*[(i, idx[v]) for i, w in enumerate(words) for _, v in patch.ring_moves(w)])
    return coo_matrix((np.r_[np.full(len(i), -.5), np.full(len(i), .5)],
                       (i+i, j+i)), shape=(len(words),)*2).tocsr()


def count_contrast(pure, rng):
    mixtures = np.array([p * pure[0] + (1-p) * pure[1] for p in WEIGHTS])
    mixtures = np.maximum(mixtures, 0)
    mixtures /= mixtures.sum(axis=1)[:, None]
    counts = np.array([rng.multinomial(SHOTS, probs) for probs in mixtures])
    return mixtures, counts


def score_slopes(scores, pure, counts):
    gradients = counts @ scores / SHOTS
    means = counts @ (scores**2) / SHOTS
    se = np.sqrt(max(0, sum(means - gradients**2)) / SHOTS) / (WEIGHTS[1]-WEIGHTS[0])
    return dict(population=float((pure[0]-pure[1]) @ scores),
                measured=float((gradients[1]-gradients[0]) / (WEIGHTS[1]-WEIGHTS[0])),
                standard_error=float(se))


def read_heads(amplitudes, labels, rng):
    pure = np.array([np.bincount(labels, abs(wave)**2, minlength=64) for wave in amplitudes.T])
    mixtures, counts = count_contrast(pure, rng)
    scored = (pure[0] - pure[1] > 1e-14).astype(float)
    binary = score_slopes(scored, pure, counts)
    ref = mixtures.mean(axis=0)
    surprise = np.zeros(len(ref))
    surprise[ref > 0] = -np.log(ref[ref > 0])
    return dict(counts=counts.tolist(), decoder_words=np.flatnonzero(scored).tolist(),
                held_law_tv=float(np.abs(mixtures[1]-mixtures[0]).sum()/2),
                significance=binary['measured']/binary['standard_error'] if binary['standard_error'] else 0.,
                entropy_slope=score_slopes(surprise, pure, counts))


def read_interior(amplitudes, rows, cols, rng):
    density = []
    for wave in amplitudes.T:
        coeff = np.zeros((rows.max()+1, cols.max()+1), dtype=complex)
        coeff[rows, cols] = wave
        density.append(coeff @ coeff.conj().T)
    ref = .75*density[0] + .25*density[1]
    values, vectors = np.linalg.eigh(ref)
    active = values > 1e-14
    values, vectors = values[active], vectors[:, active]
    pure = np.array([np.sum(vectors.conj() * (rho @ vectors), axis=0).real for rho in density])
    np.testing.assert_allclose(pure.sum(axis=1), 1, atol=2e-12, rtol=0)
    # Resolve each distinct eigenvalue only. Basis choice within a degeneracy is immaterial.
    groups = np.split(np.arange(len(values)), np.where(np.diff(values)>1e-13)[0]+1)
    scores = -np.log([values[g].mean() for g in groups])
    grouped = np.array([pure[:, g].sum(axis=1) for g in groups]).T
    _, counts = count_contrast(grouped, rng)
    return dict(entropy_slope=score_slopes(scores, grouped, counts),
                entropy_scores=scores.tolist(), counts=counts.tolist())


def forgotten_head_read(initial, head, h, t, labels, decoder, rng):
    # Keep both dephasing outcomes and both phase preparations, then trace the pointer.
    branches = np.column_stack([initial[:, sign] * (head == outcome)
                                for sign in range(2) for outcome in (-1, 1)])
    propagated = expm_multiply(-1j*t*h, branches) if t else branches
    word_laws = (abs(propagated)**2).reshape(len(head), 2, 2).sum(axis=2)
    pure = np.array([np.bincount(labels, law, minlength=64) for law in word_laws.T])
    mixtures, counts = count_contrast(pure, rng)
    score = np.zeros(64)
    score[decoder] = 1
    held = score_slopes(score, pure, counts)
    return dict(held_law_tv=float(np.abs(mixtures[1]-mixtures[0]).sum()/2),
                decoder_slope=held, counts=counts.tolist())


def run():
    patch = Patch()
    words = np.array(patch.field_component(), dtype=np.int64)
    h = ring_laplacian(patch, words)
    ground = np.ones(len(words)) / np.sqrt(len(words))
    np.testing.assert_allclose(h @ ground, 0, atol=2e-15)
    incidence = Counter(tuple(sorted(e)) for f in patch.disk for e in around(patch.faces[f]))
    inside = [e for e in patch.edges if incidence[e] == 2]
    boundary = [tuple(sorted(e)) for e in around(patch.boundary)]
    labels = sum(((words & patch.bit[e]) != 0).astype(int)*(1<<i) for i,e in enumerate(boundary))
    mask = sum(patch.bit[e] for e in inside)
    _, rows = np.unique(words & mask, return_inverse=True)
    _, cols = np.unique(words & ~mask, return_inverse=True)
    output = []
    for e_id, edge in enumerate(inside):
        head = np.where((words & patch.bit[edge]) != 0, -1., 1.)
        plus = ground * np.exp(1j * math.pi * head / 4)
        initial = np.stack([plus, plus.conj()], axis=1)
        energies = np.sum(initial.conj() * (h @ initial), axis=0).real
        np.testing.assert_allclose(energies[0], energies[1], atol=1e-14)
        for ti, t in enumerate(TIMES):
            evolved = expm_multiply(-1j*t*h, initial) if t else initial
            np.testing.assert_allclose(np.sum(abs(evolved)**2, axis=0), 1, atol=1e-12)
            inside_obs = read_interior(evolved, rows, cols, np.random.default_rng([712953,e_id,ti,0]))
            cut_obs = read_heads(evolved, labels, np.random.default_rng([712953,e_id,ti,1]))
            dephased = forgotten_head_read(initial, head, h, t, labels, cut_obs['decoder_words'],
                                         np.random.default_rng([712953,e_id,ti,2]))
            output.append(dict(forgotten_head_control=dephased, interior_seam=edge, exposure=t, mean_ring_energies=energies.tolist(),
                               interior=inside_obs, perimeter=cut_obs))
    held = [r for r in output if r['exposure']==1]
    bulk = np.array([r['interior']['entropy_slope']['population'] for r in held])
    edge = np.array([r['perimeter']['entropy_slope']['population'] for r in held])
    coefficient = -bulk[0]/edge[0] if abs(edge[0]) > 1e-10 else None
    best = -(bulk @ edge)/(edge @ edge)
    measured_bulk = np.array([r['interior']['entropy_slope']['measured'] for r in held])
    measured_edge = np.array([r['perimeter']['entropy_slope']['measured'] for r in held])
    measured_best = -(measured_bulk @ measured_edge)/(measured_edge @ measured_edge)
    return dict(empirical_optimal_coefficient=float(measured_best),
                empirical_optimal_relative_residual=float(
                    np.linalg.norm(measured_bulk + measured_best*measured_edge)/np.linalg.norm(measured_bulk)),
                source='positive written ring Laplacian, one neutral 458-word component',
                shots_per_preparation=SHOTS, mixture_weights=WEIGHTS, seed=712953,
                boundary_edges=boundary, interior_seams=inside, records=output,
                calibration_first_seam_coefficient=coefficient,
                largest_held_calibrated_residual=float(np.max(abs(bulk[1:]+coefficient*edge[1:])))
                if coefficient is not None else None, optimal_common_coefficient=float(best),
                optimal_relative_residual=float(np.linalg.norm(bulk+best*edge)/np.linalg.norm(bulk)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    result = run()
    with parser.parse_args().output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
