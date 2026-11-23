"""An entropy-independent distance chart on labeled material endpoints."""

from itertools import combinations

import numpy as np


def square_distances(features):
    differences = features[:,None,:]-features[None,:,:]
    return np.sum(differences**2,axis=2)


def gram(squares):
    mean = squares.mean(axis=0)
    return (mean[:,None]+mean[None,:]-squares-mean.mean())/2


def chart(squares):
    eigenvalues = np.linalg.eigvalsh(gram(squares))[::-1]
    values, axes = np.linalg.eigh(gram(squares[:8,:8]))
    anchor_coordinates = axes[:,-3:] * np.sqrt(np.maximum(values[-3:],0))
    held_coordinates = []
    norms = np.sum(anchor_coordinates**2,axis=1)
    for station in range(8,12):
        offsets = squares[station,:8]-norms
        target = -(offsets-offsets.mean())/2
        held_coordinates.append(np.linalg.lstsq(anchor_coordinates,target,rcond=None)[0])
    projected = square_distances(np.array(held_coordinates))
    truth = squares[8:,8:]
    pairs = np.array(list(combinations(range(4),2))).T
    expected=np.sqrt(truth[tuple(pairs)])
    estimated=np.sqrt(projected[tuple(pairs)])
    return dict(gram_spectrum=eigenvalues.tolist(),
                spectral_fraction_beyond_three=float(np.maximum(eigenvalues[3:],0).sum()/eigenvalues.sum()),
                held_expected=expected.tolist(),held_forecast=estimated.tolist(),
                largest_held_relative_error=float(np.max(abs(estimated/expected-1))))


def stationary(kernel):
    a = kernel.T - np.eye(len(kernel))
    a[-1,:] = 1
    b = np.zeros(len(kernel))
    b[-1] = 1
    weights = np.linalg.solve(a,b)
    if min(weights) <= 0 or max(abs(weights@kernel-weights))>2e-12:
        raise ArithmeticError("Vacuum terminal kernel has no positive stationary weights")
    return weights


def distances(fresh, free_vacuum, held):
    positive = np.random.default_rng(56123).normal(size=(12,3))
    control = chart(square_distances(positive))
    if control['largest_held_relative_error']>1e-11:
        raise ArithmeticError("Exact three-dimensional calibration failed")
    rows=[]
    for si, sign in enumerate((-1,1,0)):
        vacuum = free_vacuum[si] if sign else free_vacuum.mean(axis=0)
        field = fresh[si] if sign else fresh.mean(axis=0)
        for ti in range(vacuum.shape[1]):
            baseline = vacuum[:,ti]
            weights = stationary(baseline)
            scale=np.sqrt(weights)
            records=[]
            for w in range(field.shape[-1]):
                exact = field[:,ti,:,w]
                counts=np.asarray(held[w//2]['endpoints']['full_endpoint_counts'])
                observed=(counts[si,:,ti,:,w%2] if sign else counts[:,:,ti,:,w%2].sum(axis=0))
                empirical=observed/observed.sum(axis=1)[:,None]
                records.append(dict(outcome=chart(square_distances(exact/scale)),
                                    counts=chart(square_distances(empirical/scale))))
            rows.append(dict(sign=sign,appointment=ti,stationary_weights=weights.tolist(),
                             vacuum=chart(square_distances(baseline/scale)),perturbations=records))
    return dict(positive_control=control,readings=rows)
