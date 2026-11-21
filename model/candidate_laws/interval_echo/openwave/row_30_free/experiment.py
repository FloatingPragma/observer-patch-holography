"""Free dispersion of a doubled load measured with local Interval Echo triangle readers."""

import argparse
from collections import deque
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import minimize_scalar
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.executable_model.carrier import seam_digest

HERE = Path(__file__).resolve().parent
TARGETS = (.12, .32, 2.9, 5., 7.)
TIMES = np.arange(513) * .375
SHOTS = 2**25
SEED = 430029


def native_hardware(car):
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    indices = [law.SRC.index((2, 0)), law.SRC.index((0, 2))]
    source = np.eye(law.DIM)[:, indices]
    block = source.T @ h @ source
    assert np.linalg.norm(h @ source - source @ block) < 1e-11
    assert np.linalg.norm(h[:, law.SRC.index((0, 0))]) < 1e-11
    degree_phase = candidate.KAPPA * contrast[indices[0], indices[0]]
    assert np.linalg.norm(candidate.KAPPA * contrast @ source - degree_phase * source) < 1e-11
    a, b = np.array(car.seams).T
    operator = sparse.coo_matrix(
        (np.repeat(block.ravel(), len(a)), (np.r_[a, a, b, b], np.r_[a, b, a, b])),
        shape=(car.n, car.n)).tocsr()
    operator += sparse.diags((6 - car.degree) * degree_phase)
    accept = np.arange(law.DIM) % 9 == law.SRC.index((1, 1))
    inputs = np.eye(law.DIM)[:, [*indices, law.SRC.index((0, 0))]]
    output = (law.echo() @ inputs)[accept]
    assert np.linalg.norm(output.T @ output - np.diag([1, 1, 0])) < 1e-12
    return operator, block, degree_phase


class Source:
    """Spread a sharp double on a fixed spanning tree and address each final phase."""

    def __init__(self, car, block, degree_phase):
        self.car = car
        self.rate = degree_phase
        self.e, self.v = np.linalg.eigh(block)
        self.order, self.children = [0], [[] for _ in range(car.n)]
        queue, visited = deque([0]), {0}
        while queue:
            p = queue.popleft()
            for q in car.neighbours[p]:
                if q not in visited:
                    self.children[p].append(q)
                    self.order.append(q)
                    queue.append(q)
                    visited.add(q)
        assert len(self.order) == car.n

    def build(self, target):
        state = np.zeros(self.car.n, complex)
        state[0] = 1.
        subtree = abs(target)**2
        for p in reversed(self.order):
            subtree[p] += sum(subtree[q] for q in self.children[p])

        def contact(a, b, duration):
            ket = state[[a, b]]
            state[[a, b]] = self.v @ (np.exp(-1j * self.e * duration) * (self.v.T @ ket))

        pulses = 0
        for p in self.order:
            children = self.children[p]
            for n, q in enumerate(children):
                retained = abs(target[p])**2 + sum(subtree[r] for r in children[n + 1:])
                duration = 2 * np.arctan2(np.sqrt(subtree[q]), np.sqrt(retained))
                contact(p, q, duration)
                pulses += 1
        for p in self.order:
            if abs(target[p]) < 1e-13:
                continue
            phase_time = (np.angle(state[p]) - np.angle(target[p])) % (2*np.pi) / self.rate
            q = p
            if self.car.degree[p] == 6:
                q = next(j for j in self.car.neighbours[p] if self.car.degree[j] == 5)
                contact(p, q, np.pi)
            state[q] *= np.exp(-1j * self.rate * phase_time)
            if q != p:
                # The doubled seam has eigenvalues 32 and 33; a second pi contact inverts the swap.
                contact(p, q, np.pi)
                pulses += 2
        error = np.linalg.norm(state - target)
        assert error < 1e-10
        return state, dict(seam_contacts=pulses, private_reference_phases=int(sum(abs(target) >= 1e-13)),
                           maximum_amplitude_error=float(np.max(abs(state - target))))


class Readers:
    def __init__(self, car):
        faces = car.faces()
        index = {edge: n for n, edge in enumerate(car.seams)}
        self.edges = np.array(car.seams)
        self.triangles = []
        self.indices = []
        for a in range(car.n):
            first = next(face for face in faces if a in face)
            b, c = [p for p in first if p != a]
            self.triangles.append([a, b, c])
            self.indices.append([index[tuple(sorted(edge))] for edge in [(a, b), (a, c), (b, c)]])
        self.indices = np.array(self.indices)

    def clicks(self, state, rng):
        probabilities = abs(state)**2
        ab = self.edges
        rates = probabilities[ab[:, 0]] + probabilities[ab[:, 1]]
        return rng.binomial(SHOTS, np.clip(rates, 0, 1))

    def ports(self, clicks):
        ab, ac, bc = clicks[self.indices].T
        return (ab + ac - bc) / (2 * SHOTS)


def spatial_measurement(readers, plus, minus, rng):
    a, b = readers.edges.T
    def dirichlet(v):
        return ((v[a] - v[b])**2).sum() / (v @ v)
    counts = [readers.clicks(state, rng) for state in [plus, minus]]
    sqrt_n = np.sqrt(len(plus))
    v = (readers.ports(counts[0]) - readers.ports(counts[1])) * sqrt_n / 2
    v_exact = (abs(plus)**2 - abs(minus)**2) * sqrt_n / 2
    return dict(incidence_lambda_exact=dirichlet(v_exact), incidence_lambda_sampled=dirichlet(v),
                plus_seam_clicks=counts[0].tolist(), minus_seam_clicks=counts[1].tolist())


def geometry_modes(lap, targets):
    energies, vectors = np.linalg.eigh(lap.toarray())
    selected = []
    for target in targets:
        e = energies[np.argmin(abs(energies - target))]
        space = vectors[:, abs(energies - e) < 1e-9]
        # The projector and the tie break fix a real mode within each eigenspace.
        visibility = (space[:12]**2).sum(axis=1)
        p = int(np.flatnonzero(visibility >= max(visibility) - 1e-12)[0])
        v = space @ space[p] / np.sqrt(visibility[p])
        assert p < 12 and np.linalg.norm(lap @ v - e*v) < 1e-10
        selected.append((e, p, v))
    assert len({round(e, 8) for e, _, _ in selected}) == len(targets)
    return selected


def fit_temporal(data):
    times, observed = TIMES[::2], data[::2]
    def fit(w):
        cols = np.column_stack((np.ones(len(times)), np.cos(w*times), np.sin(w*times)))
        coeff = np.linalg.lstsq(cols, observed, rcond=None)[0]
        delta = cols @ coeff - observed
        return delta @ delta, coeff
    grid = np.linspace(.003, 4.15, 1660)
    i = int(np.argmin([fit(w)[0] for w in grid]))
    w = minimize_scalar(lambda freq: fit(freq)[0], bounds=(grid[max(0, i-1)], grid[min(len(grid)-1, i+1)]),
                        method="bounded", options={"xatol": 1e-13}).x
    offset, cosine, sine = fit(w)[1]
    hold = TIMES[1::2]
    difference = offset + cosine*np.cos(w*hold) + sine*np.sin(w*hold) - data[1::2]
    return dict(beta=w, held_rms=float(np.sqrt(np.mean(difference**2))),
                baseline=offset, cosine=cosine, sine=sine)


def temporal_measurement(readers, port, states, rng):
    a, b, c = readers.triangles[port]
    p = abs(states[:, [a, b, c]])**2
    rates = np.column_stack((p[:, 0]+p[:, 1], p[:, 0]+p[:, 2], p[:, 1]+p[:, 2]))
    clicks = rng.binomial(SHOTS, np.clip(rates, 0, 1))
    def reconstruct(s):
        return (s[:, 0] + s[:, 1] - s[:, 2]) / 2
    return dict(triangle=[a, b, c], exact=fit_temporal(reconstruct(rates)),
                sampled=fit_temporal(reconstruct(clicks/SHOTS)), triangle_clicks=clicks.tolist())


def run():
    car = federation(0)
    h, block, degree_phase = native_hardware(car)
    lap = car.laplacian()
    rest = 6*degree_phase
    shifted = h - rest*sparse.eye(car.n)
    identity_error = abs(shifted - lap/2).max()
    assert identity_error < 1e-10
    reference = np.ones(car.n)/np.sqrt(car.n)
    source, readers = Source(car, block, degree_phase), Readers(car)
    modes = geometry_modes(lap, TARGETS)
    preparations, rows = [], []
    rng = np.random.default_rng(SEED)
    for target, (e, p, spatial) in zip(TARGETS, modes):
        plus, plus_circuit = source.build((reference + spatial)/np.sqrt(2))
        minus, minus_circuit = source.build((reference - spatial)/np.sqrt(2))
        rows.append(dict(target=target, selected_lambda=e, port=p,
                         preparation=dict(plus=plus_circuit, minus=minus_circuit),
                         spatial=spatial_measurement(readers, plus, minus, rng)))
        preparations.append(plus)
    flight = expm_multiply(-1j*shifted, np.column_stack(preparations), start=0.,
                           stop=TIMES[-1], num=len(TIMES))
    for n, row in enumerate(rows):
        temporal = temporal_measurement(readers, row["port"], flight[:, :, n], rng)
        row["temporal"] = temporal
        for kind in ["exact", "sampled"]:
            observed_lambda = row["spatial"]["incidence_lambda_" + kind]
            temporal[kind]["beta_over_lambda"] = temporal[kind]["beta"] / observed_lambda
        assert abs(temporal["exact"]["beta_over_lambda"] - .5) < 1e-6
        assert abs(temporal["sampled"]["beta_over_lambda"] - .5) < .007
    frequencies = np.array([r["temporal"]["sampled"]["beta"] for r in rows])
    measured = np.array([r["spatial"]["incidence_lambda_sampled"] for r in rows])
    comparisons = []
    for mass in (1., 10., 192., 1000.):
        c2 = frequencies[0]*(2*mass + frequencies[0])/measured[0]
        predicted = c2*measured/(np.sqrt(mass**2 + c2*measured) + mass)
        comparisons.append(dict(assumed_rest_gap=mass, speed_squared_from_first_mode=c2,
                                held_relative_errors=(abs(predicted[1:] - frequencies[1:])/frequencies[1:]).tolist()))
    return dict(candidate=candidate.NAME, cells=car.cells, ports=car.n, seams=len(car.seams),
                seam_list_sha256=seam_digest(car), kappa=candidate.KAPPA, starting_port=0,
                shots_per_seam_setting=SHOTS, random_seed=SEED, exposure_step=.375,
                samples=len(TIMES), even_fit=len(TIMES[::2]), odd_hold=len(TIMES[1::2]),
                unshifted_uniform_energy=rest, matrix_identity_max_error=float(identity_error),
                modes=rows, finite_gap_comparisons=comparisons)


def portable(x):
    if isinstance(x, dict):
        return {k: portable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [portable(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return round(float(x), 9) + 0.
    return x


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a separate receipt")
    output = parser.parse_args().output
    record = portable(run())
    if output:
        with output.open("x") as dst:
            dst.write(json.dumps(record, indent=1, allow_nan=False) + "\n")
    elif record != json.loads((HERE / "results.json").read_text()):
        raise SystemExit("free doubled-load dispersion differs from the archived run")
    for mode in record["modes"]:
        print(mode["spatial"]["incidence_lambda_sampled"], mode["temporal"]["sampled"]["beta"])
