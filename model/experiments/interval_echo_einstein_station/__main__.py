"""Local tomography and a fixed-cut test of the G4 entropy premise."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
SEED, BLOCKS, TRIALS = 2026112110, 32, 16384
STEPS = (1 / 128, 1 / 256, 1 / 512, 1 / 1024)
SETTINGS = ("A", "B", "half_AB_read_A", "delay_BD_half_AB_read_A")


def hardware(level):
    car = federation(level)
    cut = [(a, c) for a, c in car.external if a < 12]
    assert len(cut) == 3
    boundary_ports = {a for a, _ in cut}
    a, c = cut[0]
    b = next(v for v in car.neighbours[a] if v < 12 and v not in boundary_ports)
    d = next(v for v in car.neighbours[b] if v < 12 and v not in boundary_ports and v != a)
    used = {a, b, c, d}
    probes = [next(v for v in car.neighbours[t] if v < 12 and v not in used) for t in (a, b)]
    return dict(cells=car.cells, region_cell=0, rails_ABCD=[a, b, c, d],
                boundary_seams=cut, vacant_read_neighbours_AB=probes)


def native_generators():
    # Single double on four rails ordered A, B, C, D; the other local state is (0,0).
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    h = h0 + candidate.KAPPA * h1
    local = [law.SRC.index((2, 0)), law.SRC.index((0, 2))]
    q = law.SRC.index((0, 0))
    isometry = np.eye(law.DIM)[:, local + [q]]
    restricted = isometry.T @ h @ isometry
    assert np.linalg.norm(h @ isometry - isometry @ restricted) < 1e-11
    assert np.linalg.norm(restricted[2]) < 1e-11
    pointer = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
    signal = (abs(law.echo()[pointer][:, [local[0], q]]) ** 2).sum(axis=0)
    assert np.allclose(signal, [1, 0], atol=1e-12, rtol=0)
    return restricted[:2, :2]


def contact(block, i, j, duration):
    u = np.eye(4, dtype=complex)
    u[np.ix_([i, j], [i, j])] = expm(-1j * duration * block)
    return u


def reduced(state, indices):
    """Vacuum direct sum the one-double block after tracing out the discarded rails."""
    occupied = state[np.ix_(indices, indices)]
    result = np.zeros((1 + len(indices), 1 + len(indices)), dtype=complex)
    result[0, 0] = np.trace(state) - np.trace(occupied)
    result[1:, 1:] = occupied
    return result


class Patch:
    def __init__(self):
        block = native_generators()
        self.half = contact(block, 0, 1, np.pi / 2)
        self.cross = contact(block, 0, 2, np.pi / 3)
        delay = contact(block, 1, 3, np.pi / 65)
        settings = [(np.eye(4), 0), (np.eye(4), 1), (self.half, 0), (self.half @ delay, 0)]
        self.effects = np.array([np.outer(u[t].conj(), u[t]) for u, t in settings])
        basis = [np.diag([1, 0]), np.diag([0, 1]), np.array([[0, 1], [1, 0]]),
                 np.array([[0, 1j], [-1j, 0]])]
        self.design = np.array([[np.trace(e[:2, :2] @ z).real for z in basis] for e in self.effects])
        assert np.linalg.matrix_rank(self.design) == 4

    def state(self, p, splitter=True):
        u = self.cross @ self.half if splitter else self.cross
        return (u * [0.8 * p, 0.8 * (1 - p), 0.2, 0]) @ u.conj().T

    def table(self, p):
        state = self.state(p)
        return np.einsum('ijk,kj->i', self.effects, state).real

    def invert(self, clicks):
        a, b, re, im = np.linalg.solve(self.design, clicks)
        return np.array([[1 - a - b, 0, 0], [0, a, re + 1j * im], [0, re - 1j * im, b]])

    def region(self, p):
        matrix = self.invert(self.table(p))
        assert np.linalg.norm(matrix - reduced(self.state(p), [0, 1])) < 1e-12
        return matrix


def entropy(matrix):
    eig = np.linalg.eigvalsh(matrix)
    if min(eig) <= 0:
        raise ValueError(f'tomography left the faithful support: {eig}')
    return float(-np.dot(eig, np.log(eig)))


def entropy_components(matrix):
    weights = [matrix[0, 0].real, np.trace(matrix[1:, 1:]).real]
    centre = -np.dot(weights, np.log(weights))
    conditional = weights[1] * entropy(matrix[1:, 1:] / weights[1])
    total = entropy(matrix)
    assert abs(centre + conditional - total) < 1e-13
    return dict(centre=float(centre), weighted_within_occupied=float(conditional),
                total=total, sector_weights=weights)


def derivatives(patch, p):
    ref = patch.region(p)
    eig, v = np.linalg.eigh(ref)
    k = (v * -np.log(eig)) @ v.conj().T
    tangent = patch.region(1) - patch.region(0)
    slope = float(np.trace(k @ tangent).real)
    estimates = [(entropy(patch.region(p + d)) - entropy(patch.region(p - d))) / (2 * d)
                 for d in STEPS]
    assert abs(estimates[-1] - slope) < 2e-5
    return dict(p=p, reference_entropy=entropy_components(ref), regional_modular_slope=slope,
                reference_spectrum=eig.tolist(), steps=STEPS, entropy_slopes=estimates,
                maximum_diagonal_derivative=float(max(abs(tangent.diagonal()))),
                reference_click_probabilities=patch.table(p).tolist())


def finite_reader(patch, p, rng):
    step = 1 / 16
    probabilities = np.array([patch.table(p - step), patch.table(p + step)])
    counts = rng.binomial(TRIALS, probabilities[:, None, :], size=(2, BLOCKS, len(SETTINGS)))
    def statistic(rates):
        return (entropy(patch.invert(rates[1])) - entropy(patch.invert(rates[0]))) / (2 * step)
    pooled = counts.sum(axis=1)
    slope = statistic(pooled / (TRIALS * BLOCKS))
    leaves = [statistic((pooled - counts[:, i]) / (TRIALS * (BLOCKS - 1))) for i in range(BLOCKS)]
    se = float(np.sqrt((BLOCKS - 1) / BLOCKS * np.sum((leaves - np.mean(leaves)) ** 2)))
    return dict(p=p, step=step, expected_finite_difference=statistic(probabilities), measured=slope,
                standard_error=se, four_se_excludes_zero=bool(abs(slope) > 4 * se),
                blocks_of_click_counts=counts.tolist())


def measure():
    patch = Patch()
    boards = [hardware(level) for level in range(3)]
    cut = [0, 2]
    def boundary_distance(splitter):
        d = reduced(patch.state(1, splitter), cut) - reduced(patch.state(0, splitter), cut)
        return float(sum(abs(np.linalg.eigvalsh(d))) / 2)
    boundary = boundary_distance(True)
    bypass = boundary_distance(False)
    assert boundary < 1e-12 and bypass > 0.1
    exact = [derivatives(patch, p) for p in (0.75, 0.5)]
    rng = np.random.default_rng(SEED)
    records = [finite_reader(patch, p, rng) for p in (0.75, 0.5)]
    return dict(law=candidate.NAME, kappa=candidate.KAPPA, hardware=boards,
                settings=SETTINGS, tomography_condition=float(np.linalg.cond(patch.design)),
                boundary_trace_distance=boundary, splitter_bypass_boundary_distance=bypass,
                exact=exact, seed=SEED, blocks=BLOCKS, trials_per_block_setting=TRIALS, finite=records)


def check(reference, measured):
    if isinstance(reference, dict):
        assert reference.keys() == measured.keys()
        for key, value in reference.items():
            check(value, measured[key])
    elif isinstance(reference, list):
        assert len(reference) == len(measured)
        for left, right in zip(reference, measured, strict=True):
            check(left, right)
    elif isinstance(reference, float):
        assert abs(reference - measured) < 1e-9
    else:
        assert reference == measured


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='deposit a fresh JSON record at this path')
    args = parser.parse_args()
    result = measure()
    if args.output is None:
        check(json.loads((HERE / 'results.json').read_text()), result)
    else:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    for exact, finite in zip(result['exact'], result['finite'], strict=True):
        print(f"p={exact['p']}: modular {exact['regional_modular_slope']:.8f}, "
              f"sample {finite['measured']:.8f} ± {finite['standard_error']:.8f}")
    print('Cut trace distances (splitter, bypass):', result['boundary_trace_distance'],
          result['splitter_bypass_boundary_distance'])


if __name__ == '__main__':
    main()
