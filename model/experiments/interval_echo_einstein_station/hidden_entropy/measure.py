"""Regional entropy and invisible boundary histories of the registered free one-load repair."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import compile_wave
from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.measure import angular_frame
from model.executable_model import federation

SHOTS = 262144
STEP = 1/16
FRACTIONS = tuple(p + step for p in (0.5, 0.75) for step in (-STEP, 0, STEP))
HERE = Path(__file__).resolve().parent


def entropy(p):
    return -p*np.log(p) - (1-p)*np.log1p(-p)


def two_dark_modes(frame, ports):
    constraints = frame[ports]
    # Project onto the gateway kernel, then choose the first two independent columns.
    kernel = np.eye(5) - np.linalg.pinv(constraints) @ constraints
    choices = []
    for col in kernel.T:
        vector = col.copy()
        for previous in choices:
            vector -= (previous @ vector)*previous
        if np.linalg.norm(vector) > 1e-9:
            choices.append(vector / np.linalg.norm(vector))
        if len(choices) == 2:
            break
    modes = frame @ np.column_stack(choices)
    np.testing.assert_allclose(modes[ports], 0., atol=1e-13)
    np.testing.assert_allclose(modes.T @ modes, np.eye(2), atol=1e-13)
    return modes


def make_states(board, descendants, cell, frame):
    states = np.zeros((board.n, 3*descendants), dtype=complex)
    readers = []
    for site in range(descendants):
        start = 12*site
        ports = np.flatnonzero(board.degree[start:start+12] == 6)
        modes = two_dark_modes(frame, ports)
        first, second = modes.T
        readers.append([compile_wave(cell, target) for target in
                        (first, second, (first+second)/np.sqrt(2), (first+1j*second)/np.sqrt(2))])
        edge_ports = [port for port in ports if any(q//12 >= descendants
                      for q in board.neighbours[start+port])]
        pole = edge_ports[0] if edge_ports else ports[0]
        bright = frame @ frame[pole] / np.linalg.norm(frame[pole])
        for j, wave in enumerate((first, second, bright)):
            states[start:start+12, 3*site+j] = compile_wave(cell, wave)
    return states, np.asarray(readers)


def changed_generator(board, kinetic):
    internal = np.asarray([edge for edge in board.seams if edge[0]//12 == edge[1]//12])
    weights = (board.degree[internal] == 6).sum(axis=1)
    incidence = coo_matrix((np.tile([1., -1.], len(internal)),
                            (np.repeat(np.arange(len(internal)), 2), internal.ravel())),
                           shape=(len(internal), board.n)).tocsr()
    return kinetic + incidence.T @ incidence.multiply(weights[:, None])/2


def terminal_records(states, readers, descendants, rng):
    # Input columns repeat (mode 0, mode 1, bright) at each privately selected site.
    outside = abs(states[12*descendants:])**2
    escape = outside.sum(axis=0).reshape(descendants, 3).mean(axis=0)
    outer_counts = [int(rng.binomial(SHOTS, np.clip(p, 0, 1))) for p in escape]
    report = {"exterior_probability": escape.tolist(), "exterior_clicks": outer_counts}
    regional = states[:12*descendants].reshape(descendants, 12, -1)
    probabilities = np.abs(np.einsum("cjp,cps->cjs", readers.conj(), regional))**2
    # Average over the private source addresses, retaining each detector's address.
    probabilities = probabilities.reshape(descendants, 4, descendants, 3).mean(axis=2)
    rows = []
    for p in FRACTIONS:
        readings = p*probabilities[:,:,0] + (1-p)*probabilities[:,:,1]
        counts = []
        for per_site in readings.T:
            distribution = np.append(per_site, max(0., 1-per_site.sum()))
            counts.append(rng.multinomial(SHOTS, distribution / distribution.sum()).tolist())
        yes = np.array([1-row[-1]/SHOTS for row in counts])
        # On the validated native code space the three components give its density matrix.
        bloch = np.array([2*yes[2]-1, 2*yes[3]-1, yes[0]-yes[1]])
        eig = (1+min(1., np.linalg.norm(bloch)))/2
        rows.append({"fraction": p, "counts_modes_by_cell_then_complement": counts,
                     "yes_fractions": yes.tolist(),
                     "two_mode_entropy_assuming_closed_support": float(entropy(eig)) if eig < 1 else 0.})
    report["regional_readings"] = rows
    return report


def run():
    cell, frame, _, _ = angular_frame()
    rng = np.random.default_rng(821467)
    boards = []
    for level in range(3):
        board = federation(level)
        descendants = 4**level
        states, readers = make_states(board, descendants, cell, frame)
        kinetic = candidate.generator(board) - 6*candidate.KAPPA*eye(board.n, format="csr")
        dark = states[:,np.arange(states.shape[1])%3 != 2]
        np.testing.assert_allclose(kinetic @ dark, 3*dark, atol=3e-13)
        runs = {}
        for name, generator in (("registered", kinetic), ("gateway_weighted", changed_generator(board, kinetic))):
            evolution = expm_multiply(-1j*generator, states, start=0, stop=4, num=3)
            observations = []
            for exposure, current in zip((0,2,4), evolution, strict=True):
                outcome = terminal_records(current, readers, descendants, rng)
                if name == "registered":
                    np.testing.assert_allclose(current[:,np.arange(states.shape[1])%3 != 2],
                                               np.exp(-3j*exposure)*dark, atol=3e-13)
                    assert max(outcome["exterior_probability"][:2]) < 1e-23
                    for row in outcome["regional_readings"]:
                        p = row["fraction"]
                        np.testing.assert_allclose(row["yes_fractions"], [p, 1-p, .5, .5], atol=.01)
                    rows = outcome["regional_readings"]
                    derivatives = [(rows[j+2]["two_mode_entropy_assuming_closed_support"]
                                    - rows[j]["two_mode_entropy_assuming_closed_support"])/(2*STEP)
                                   for j in (0, 3)]
                    exact = [(entropy(p+STEP)-entropy(p-STEP))/(2*STEP) for p in (.5,.75)]
                    np.testing.assert_allclose(derivatives, exact, atol=.04)
                    outcome["sampled_entropy_slopes"] = derivatives
                observations.append(dict(exposure=exposure, **outcome))
            runs[name] = observations
        assert runs["registered"][-1]["exterior_probability"][2] > .001
        assert max(runs["gateway_weighted"][-1]["exterior_probability"][:2]) > 1e-7
        boards.append(dict(cells=board.cells, descendants=descendants, runs=runs))
    return dict(shots_per_instrument=SHOTS, seed=821467, half_step=STEP, boards=boards)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the complete observer records to a fresh path")
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, separators=(",", ":"), allow_nan=False)
            stream.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import check_replay
        check_replay(json.loads((HERE/"results.json").read_text()), result)
    for board in result["boards"]:
        print(board["cells"], board["runs"]["registered"][-1]["sampled_entropy_slopes"],
              {name: rows[-1]["exterior_clicks"] for name, rows in board["runs"].items()})
