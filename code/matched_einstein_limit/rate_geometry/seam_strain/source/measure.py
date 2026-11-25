"""Native two-head sources pumping the reciprocal whole-bridge strain qubit."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station.collar_pairing import (
    measure as pairing,
)
from ....source_implication.directional_stress import bridge

formation = pairing.formation
HERE = Path(__file__).resolve().parent
TIMES = (0.5, 1.0, 2.0, 4.0)
EPSILON = np.log(5 / 4)
RAILS = {
    "reference": (2 / 3, 1 / 3),
    "ordinary": (1 / 3, 2 / 3),
    "decorrelated": (2 / 3, 1 / 3),
}


def system(board, route, charge, coupled=True):
    """Return integer 40H and twice the native whole connecting bridge."""
    h2 = bridge.native.generator(board, route, charge, False).astype(np.int64)
    b2 = bridge.bridges(board, [route[1:3]])
    identity = sparse.eye(h2.shape[0], dtype=np.int64, format="csr")
    blocks = [
        [20 * h2 - (4 * b2 if coupled else 0), -20 * identity],
        [-20 * identity, 20 * h2 + (5 * b2 if coupled else 0)],
    ]
    a = sparse.bmat(blocks, format="csr")
    assert (a.T != a).nnz == 0
    return a, h2, b2


def launch(board, route):
    """Twice the isometry for two rail modes, four native head words and |+X>."""
    n = 4 * board.n
    half = np.zeros((n, 8), dtype=np.int64)
    eye = np.eye(4, dtype=np.int64)
    for port, sign in zip(route[1:3], (1, -1)):
        half[4 * port : 4 * port + 4, :] = np.column_stack([eye, sign * eye])
    return np.vstack([half, half])


def form(board, route, charge):
    rho = pairing.prepared_states(board, route, charge)[0]
    rho = rho / np.trace(rho).real
    return rho, np.kron(*formation.marginals(rho))


def bare_records(h2, b2, initial, states):
    n = h2.shape[0]
    # The integer launch is twice V, so its upper block is sqrt(2) times the source isometry.
    rails = initial[:n] / np.sqrt(2)
    positive = (2 / 3, 1 / 3)
    means = {}
    log_rho = pairing.log_density(states[0])
    for name, weights in RAILS.items():
        rho = states[name == "decorrelated"]
        source = np.kron(np.diag(weights), rho)
        values = [np.trace(rails.T @ (op @ rails) @ source).real / 2 for op in (h2, b2)]
        modular = -np.trace(log_rho @ rho).real - np.dot(weights, np.log(positive))
        means[name] = dict(
            bare_energy=values[0],
            whole_bridge_energy=values[1],
            coupled_energy=values[0] - 0.5 + (41 / 40 - 1) * values[1],
            modular_nats=modular,
        )
    return means


def compressed(wave):
    """Pairwise qubit blocks evaluated on the entire unobserved environment."""
    qubits = wave.reshape(2, -1, 8)
    return np.array([[right.conj().T @ left for right in qubits] for left in qubits])


def read(blocks, rho, weights):
    source = np.kron(np.diag(weights), rho)
    density = np.einsum("abij,ji->ab", blocks, source)
    bloch = np.array(
        [
            2 * density[0, 1].real,
            -2 * density[0, 1].imag,
            (density[0, 0] - density[1, 1]).real,
        ]
    )
    return dict(
        bloch=bloch.tolist(),
        density_real=density.real.tolist(),
        density_imag=density.imag.tolist(),
    )


def fast_flow(a, initial):
    previous = 0
    for time in TIMES:
        dt = time - previous
        initial = expm_multiply(
            (-1j * dt / 40) * a, initial, traceA=(-1j * dt / 40) * a.diagonal().sum()
        )
        yield time, initial
        previous = time


def run():
    runs, controls, formations = [], [], []
    for level in (0, 1):
        board = formation.federation(level)
        route = formation.region(board)
        start = launch(board, route)
        states = {s: form(board, route, s) for s in (-1, 1)}
        for s, (rho, product) in states.items():
            formations.append(
                dict(
                    level=level,
                    cells=board.cells,
                    forming_charge=s,
                    source_support_minimum=float(np.linalg.eigvalsh(rho)[0]),
                    head_pauli_probabilities=dict(
                        formed=formation.probabilities(rho).real.tolist(),
                        decorrelated=formation.probabilities(product).real.tolist(),
                    ),
                )
            )
        for q in (-1, 1):
            control = system(board, route, q, coupled=False)[0]
            for time, wave in fast_flow(control, start / 2):
                plus, minus = np.split(wave, 2)
                controls.append(
                    dict(
                        level=level,
                        cells=board.cells,
                        pump_charge=q,
                        time=time,
                        largest_minus_x_projection_norm=float(
                            np.linalg.norm((plus - minus) / np.sqrt(2), axis=0).max()
                        ),
                    )
                )
            a, h2, b2 = system(board, route, q)
            bare = {s: bare_records(h2, b2, start, pair) for s, pair in states.items()}
            for time, wave in fast_flow(a, start / 2):
                effects = compressed(wave)
                for s, pair in states.items():
                    arms = {
                        name: dict(
                            **read(effects, pair[name == "decorrelated"], weights),
                            **bare[s][name],
                        )
                        for name, weights in RAILS.items()
                    }
                    baseline = arms["reference"]
                    o, c = [
                        np.array(arms[name]["bloch"]) - baseline["bloch"]
                        for name in ("ordinary", "decorrelated")
                    ]
                    ratio = (
                        arms["decorrelated"]["modular_nats"] - baseline["modular_nats"]
                    ) / (arms["ordinary"]["modular_nats"] - baseline["modular_nats"])
                    runs.append(
                        dict(
                            cells=board.cells,
                            level=level,
                            route=list(route),
                            forming_charge=s,
                            pump_charge=q,
                            time=time,
                            epsilon=EPSILON,
                            omega=0.5,
                            arms=arms,
                            transfer=dict(
                                modular_ratio=ratio,
                                ordinary=o.tolist(),
                                decorrelated=c.tolist(),
                                held_residual=(c - ratio * o).tolist(),
                                orthogonal_distance=float(
                                    np.linalg.norm(np.cross(c, o)) / np.linalg.norm(o)
                                ),
                            ),
                        )
                    )
    return dict(runs=runs, zero_epsilon_controls=controls, formations=formations)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(out)
    else:
        print(out, end="")
