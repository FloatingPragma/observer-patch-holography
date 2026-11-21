"""Compare the modular charge of the pair-source mixture with the cell energy meter."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import block_diag, eigh

from model.executable_model import cell, federation
from model.experiments.interval_echo_null_balance.experiment import hardware, states
from .__main__ import HERE, check

TIME = 4
SHOTS, SEED = 2**18, 2026112111
BLOCKS = (slice(0, 1), slice(1, 13), slice(13, 79))


def marginal(sector, wave):
    """Partial trace grouped by the exact exterior occupancy word."""
    a, b = sector.aa, sector.bb
    two = wave[b < 12]
    one = np.zeros((12, 24), complex)
    mixed = (a < 12) & (b >= 12)
    one[a[mixed], b[mixed] - 12] = wave[mixed]
    return block_diag(np.sum(abs(wave[a >= 12]) ** 2), one @ one.conj().T, np.outer(two, two.conj()))


def rays(rho, rng=None):
    """Return the local matrix reconstructed from independent word and two-ray readings."""
    diagonal = rho.diagonal().real
    diagonals = diagonal if rng is None else rng.multinomial(SHOTS, diagonal / diagonal.sum()) / SHOTS
    result = np.diag(diagonals).astype(complex)
    data = {'diagonal': diagonals}
    for n, s in enumerate(BLOCKS[1:], start=1):
        i, j = np.triu_indices(s.stop - s.start, k=1)
        i, j = i + s.start, j + s.start
        population = (diagonal[i] + diagonal[j]) / 2
        z = rho[i, j]
        for which, delta, factor in (('real', z.real, 1), ('imag', -z.imag, -1j)):
            p, q = population + delta, population - delta
            probabilities = np.column_stack((p, q, 1 - p - q))
            assert np.min(probabilities) > -1e-13
            if rng is None:
                difference = p - q
            else:
                probabilities = np.maximum(probabilities, 0.)
                counts = rng.multinomial(SHOTS, probabilities / probabilities.sum(axis=1, keepdims=True))
                difference = (counts[:, 0] - counts[:, 1]) / SHOTS
                data[f'block{n}_{which}'] = counts
            result[i, j] += factor * difference / 2
        result[j, i] = result[i, j].conj()
    return result, data


def physical_word_readers(laboratory, samples):
    """Compile the physically equivalent rays with the workspace's native two-load compiler."""
    import networkx as nx

    stress = HERE.parents[3] / 'LAW/current-interval-clock/theorem-calibration/search/stress'
    if not (stress / 'meter.py').is_file():
        raise FileNotFoundError('The optional native compiler replay requires the LAW research workspace.')
    sys.path.insert(0, str(stress))
    import meter as historical

    board = nx.Graph()
    for a, b in federation(0).external:
        board.add_edge(a // 12, b // 12, ports={a // 12: a, b // 12: b})
    geom = nx.Graph(cell().seams)
    detector = historical.probe_hardware(board, 0, set(laboratory), geom)
    compiler = historical.BinaryCompiler(cell())
    compiler.verify_native_macros()
    plans, errors, costs = [], [], []
    for n, s in enumerate(BLOCKS[1:], start=1):
        rows, error, cost = historical.click_rows(n, compiler, geom)
        plans.append((s.start, s.stop, rows))
        errors.append(error)
        costs.append(cost)
    err = max(np.linalg.norm(historical.tomography(x, plans) - x, 2) for x in samples)
    assert err < 1e-11 and set(detector['owner_cells'][1:]).isdisjoint(laboratory)
    return dict(reader_cells=detector['owner_cells'], flight_cells=laboratory,
                all_words=detector['sharp_occupancy_labels_checked'], settings=1+sum(len(p[-1]) for p in plans),
                projector_errors=errors, max_elementary_ops_per_setting=costs,
                max_reconstruction_error=err, hardware=detector)


def contrast(reference, samples, energies, cutoff):
    w, v = np.linalg.eigh(reference)
    mask = w > cutoff
    kept = v[:, mask]
    proj = kept @ kept.conj().T
    k = (kept * -np.log(w[mask])) @ kept.conj().T
    charges = np.array([np.trace(k @ x).real for x in samples])
    dk = charges[::2] - charges[1::2]
    de = energies[::2] - energies[1::2]
    if abs(dk[0]) < 1e-7 or min(abs(de)) < 1e-7:
        gain, predicted, residual, relative, valid = (None,) * 5
    else:
        gain = float(de[0] / dk[0])
        predicted = float(gain * dk[1])
        residual = predicted - de[1]
        relative = residual / abs(de[1])
        valid = bool(abs(residual) <= max(1e-9, 1e-5 * max(abs(de))))
    return dict(cutoff=cutoff, ranks=[int(sum(np.linalg.eigvalsh(reference[s,s]) > cutoff)) for s in BLOCKS],
                min_retained_eigenvalue=float(min(w[mask])), omitted_weight=float(np.sum(w[~mask])),
                max_outside_support_frobenius=max(float(np.linalg.norm(x - proj @ x @ proj)) for x in samples),
                charge_plus_minus=charges.tolist(), modular_contrasts=dk.tolist(), energy_contrasts=de.tolist(),
                port0_gain=gain, port9_prediction=predicted, port9_residual=residual,
                residual_over_abs_held=relative, matches_held=valid)


def measure(verify_compiler=False):
    sector, nodes, links, meter, cell_h, bonds = hardware()
    h = sum(cell_h) + sum(b.h for b in bonds)
    eigen, v = eigh(h.toarray())
    keys, samples, energies = [], [], []
    for port, sign, start in states(sector):
        wave = v @ (np.exp(-1j * TIME * eigen) * (v.T @ start))
        local = marginal(sector, wave)
        recovered, _ = rays(local)
        assert np.linalg.norm(recovered - local) < 1e-12
        keys.append([port, sign])
        samples.append(recovered)
        energies.append(meter.read(wave))
    assert keys == [[0,1], [0,-1], [9,1], [9,-1]]
    energies = np.array(energies)
    reference = np.mean(samples, axis=0)
    report = dict(time=TIME, cells=nodes, links=links, launch_keys=keys, energy=energies.tolist(),
                  reference_sector_weights=[float(np.trace(reference[s,s]).real) for s in BLOCKS],
                  cutoffs=[contrast(reference, samples, energies, tol) for tol in (1e-14, 1e-12, 1e-10)])
    # A local instrument which records only the twelve-port binary word.
    diagonal = [np.diag(np.diag(x)) for x in samples]
    report['diagonal_only'] = contrast(np.mean(diagonal, axis=0), diagonal, energies, 1e-12)
    finite, data = rays(reference, np.random.default_rng(SEED))
    vals = np.linalg.eigvalsh(finite)
    report['reference_survey'] = dict(shots_per_setting=SHOTS, seed=SEED, settings=4423,
                                      min_eigenvalue=float(vals.min()), negative_eigenvalues=int(sum(vals < -1e-12)),
                                      frobenius_error=float(np.linalg.norm(finite - reference)))
    data['diagonal'] = np.rint(SHOTS * data['diagonal']).astype(np.int64)
    if verify_compiler:
        report['compiler'] = physical_word_readers(nodes, samples)
    return report, data


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--verify-compiler', action='store_true')
    cli.add_argument('--write', action='store_true')
    args = cli.parse_args()
    report, data = measure(args.verify_compiler)
    path = HERE / 'pair.json'
    if args.write:
        if not args.verify_compiler:
            raise ValueError('Archive the native reader replay alongside the source results.')
        path.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
        np.savez_compressed(HERE / 'pair_reads.npz', **data)
    else:
        ref = json.loads(path.read_text())
        if not args.verify_compiler:
            del ref['compiler']
        check(ref, report)
        with np.load(HERE / 'pair_reads.npz') as archived:
            assert set(archived.files) == set(data)
            for name, array in data.items():
                np.testing.assert_array_equal(archived[name], array)
    for key in ('cutoffs', 'diagonal_only', 'reference_survey'):
        print(key, report[key])


if __name__ == '__main__':
    main()
