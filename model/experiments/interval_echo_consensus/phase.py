"""Read relative phase after a strict priority repair, with both-probe arbitration as a control."""
import argparse
import hashlib
import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import federation
from .__main__ import (DIMENSIONS, HERE, INPUTS, TOL, act, branches, check_outputs,
                       distance, layout, local_ops, source_map)


def priority_window(source, arrival, priority, hit, probes):
    """Buffer both requests. Only interrogate the other address if the priority address declines."""
    if sorted(arrival) != [0, 1]:
        raise ValueError('both distinct requests are required in this service window')
    pending = set(arrival)
    pending.remove(priority)
    second, = pending
    outputs = {None: [], 0: [], 1: []}
    undecided = [source]
    for seam in (priority, second):
        quiet = []
        for state in undecided:
            for ready, op in probes:
                branch = act(state, op, seam)
                if np.linalg.norm(branch) < 1e-13:
                    continue
                if ready:
                    committed = act(branch, hit, seam)
                    outputs[seam].append(committed.reshape(-1, len(INPUTS)))
                else:
                    quiet.append(branch)
        undecided = quiet
    outputs[None] = [v.reshape(-1, len(INPUTS)) for v in undecided]
    return outputs


def rail_read(outputs, seam, phase, undo_and_mix):
    basis = np.eye(len(INPUTS))
    ket = (basis[seam] + phase * basis[seam+1]) / np.sqrt(2)
    result = 0.
    for branch in outputs[seam]:
        state = (branch @ ket).reshape(DIMENSIONS)
        out = act(state, undo_and_mix, seam)
        # The first endpoint is a sharp double, read against a supplied blank neighbouring port.
        result += np.sum(np.abs(np.take(out, 2, axis=2+seam))**2)
    return float(result)


def run():
    echo, _, hit, probes = local_ops()
    h0, h1 = native_law.local_term_parts()
    # On a blank collar the agreement read against a vacant rail is exactly the double projector.
    pointer = echo.reshape(11, 3, 3, 99)[:, 1, 1, :]
    columns = [native_law.SRC.index((n, 0)) for n in range(3)]
    observed_effect = pointer[:, columns].T @ pointer[:, columns]
    if np.linalg.norm(observed_effect - np.diag([0., 0., 1.])) > TOL:
        raise ArithmeticError('the native vacancy reader differs from a sharp double read')
    mix = expm(-0.5j * np.pi * (h0 + KAPPA*h1)) @ echo.T
    source = source_map()
    # Check that the rail/vacancy read is available in the actual observer patch.
    boards = []
    for level in range(3):
        info = layout(level)
        car = federation(level)
        unused = set(info['path']) | set(info['token_rails'])
        vacancies = [next(p for p in car.neighbours[a] if p not in unused and p//12 == a//12)
                     for a in info['path'][:2]]
        boards.append(dict(cells=car.cells, rail_vacancies=vacancies))
    rows = []
    for priority in (0, 1):
        first, reverse = [priority_window(source, order, priority, hit, probes)
                          for order in ((0, 1), (1, 0))]
        check_outputs(first, probes)
        check_outputs(reverse, probes)
        channel_error = distance(first, reverse)
        full = branches(source, (0, 1), priority, hit, probes)
        single = [rail_read(first, priority, phase, mix) for phase in (1j, -1j)]
        both = [rail_read(full, priority, phase, mix) for phase in (1j, -1j)]
        if channel_error > TOL or abs(single[1]-single[0]) < 1-TOL or abs(both[1]-both[0]) > TOL:
            raise ArithmeticError('phase control or order comparison failed')
        rows.append(dict(priority=priority, choi_frobenius_distance=round(channel_error, 12),
                         priority_probe_rail_probabilities=np.round(single, 12).tolist(),
                         both_probes_rail_probabilities=np.round(both, 12).tolist()))
    return dict(protocol_sha256=hashlib.sha256((HERE/'FOLLOWUP.md').read_bytes()).hexdigest(),
                rail_vacancies=boards, phase_order=['+i','-i'], priority_rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', action='store_true', help='create the local phase receipt')
    result = run()
    path = HERE/'phase.json'
    if parser.parse_args().record:
        with path.open('x', encoding='utf-8') as f:
            f.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    elif result != json.loads(path.read_text()):
        raise SystemExit('the phase test differs from the reference')
    print(json.dumps(result, indent=2))
