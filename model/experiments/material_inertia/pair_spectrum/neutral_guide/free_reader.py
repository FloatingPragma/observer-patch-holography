"""Read unrestricted recoil and distant signed charge on the same four mounted sources."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from ... import neutral_defect
from .. import measure as spectrum
from ..sum_rule_inertia import measure as meter
from .guide import physical_rings, source_vector


def free_sources(ring, n):
    size = len(ring)
    v, r = np.indices((size, size - 1))
    ids = neutral_defect.ordered_index(np.asarray(ring)[(v + r + 1) % size], np.asarray(ring)[v], n)
    for sign in (1, -1):
        for off in (0, 1):
            k = size // 2 + off
            source = np.zeros(n * (n - 1), complex)
            source[ids] = (np.exp(2j * np.pi * k * v / size) / np.sqrt(size)) * source_vector(size, k, sign)
            yield sign, off, source


def charge_after_flight(e, basis, g, admission, x, y, rng):
    q = (x // 12 == 0).astype(int) - (y // 12 == 0).astype(int)
    tick = 10**6
    arrived = basis[:, :2] @ (g * np.exp(-1j * (e[:2] - e[0]) * tick))
    p = np.bincount(q + 1, weights=abs(arrived)**2, minlength=3)
    records = {}
    for name, bright in [('parity', p[0] + p[2]), ('sign', (1 - (p[2] - p[0])) / 2)]:
        records[name] = meter.sample([admission * bright, admission * (1 - bright), 1 - admission], rng)
    return dict(owner=0, time=tick, true_charge_probabilities_from_minus_to_plus=p.tolist(),
                receiver_exits=['20', '02', 'initial rest-filter refused'], records=records)


def report():
    board = federation(0)
    _, rings, _ = physical_rings()
    x, y, _, _, h, cost = neutral_defect.generators(board)
    reduced = spectrum.PairReduction(h, cost)
    e, basis, residual, _ = reduced.lowest(h, 24)
    cut = (e[1] + e[2]) / 2
    if [reduced.inertia(c)['upper_lines_below'] for c in (cut, spectrum.CUT)] != [2, 24]:
        raise ArithmeticError('unrestricted spectral admission is ambiguous')
    parent_receipt = json.loads((spectrum.HERE / 'results.json').read_text())
    band_error = max(parent_receipt['cluster_projector_error'], residual / (spectrum.CUT - e[-1])) + 1e-12
    rest_error = residual / (min(cut - e[1], spectrum.CUT - e[-1]) - residual) + 1e-12
    operator_error = 8 * np.ptp(e) * (rest_error + band_error) + 4 * residual + 240 * band_error**2
    edges, all_profiles = meter.spatial(board)
    profiles = [p for p in all_profiles if p['center'] == 0 and p['shell'] in (1, 3)]
    rng = np.random.default_rng(190028)
    results = []
    for sign, off, source in free_sources(rings[20]['ports'], board.n):
        g = basis[:, :2].T @ source
        admission = float(np.vdot(g, g).real)
        if np.sqrt(admission) <= 2 * rest_error:
            results.append(dict(orientation=sign, offset=off, admission=admission,
                                unresolved='rest projection smaller than spectral error'))
            continue
        g /= np.sqrt(admission)
        psi = basis[:, :2] @ g
        coefficient_error = min(2., 2 * rest_error / np.sqrt(admission))
        errors = []
        for owner in range(board.cells):
            c = ((x // 12 == owner).astype(float) + (y // 12 == owner).astype(float)) / 2
            local = basis[:, :2].T @ (c[:, None] * basis[:, :2])
            errors.append(2 * rest_error + np.sqrt(max(np.linalg.eigvalsh(local))) * coefficient_error)
        panel = meter.build_source(e, basis, psi, g, admission, x, y, edges, profiles,
                                   np.array(errors), coefficient_error, operator_error, rng)
        panel.update(label=f'orientation_{sign}_offset_{off}', orientation=sign, offset=off,
                     signed_remote=charge_after_flight(e, basis, g, admission, x, y, rng),
                     remote_effect_probability_error=float(coefficient_error + 2 * rest_error + 10**6 * residual))
        first, third = [meter.intersect(row['records']) for row in panel['response']['native']]
        panel['within_source'] = dict(first=first, third=third, result=meter.decision(first, third))
        results.append(panel)
    verdict = None if 'unresolved' in results[0] else meter.inference([r for r in results if 'unresolved' not in r])
    return dict(seed=190028, shots_per_setting=meter.SHOTS, geometry_edges=edges.tolist(), profiles=profiles,
                band_projector_error=band_error, rest_projector_error=rest_error,
                normalized_source_projector='lowest native doublet in the unrestricted 24-line band',
                calibration_and_held=verdict, mounted_sources=results)


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--output', type=Path, required=True)
    args = cli.parse_args()
    args.output.write_text(json.dumps(report(), indent=1, allow_nan=False) + '\n')
