"""Compare two probe orders on the complete two-unit, blank-collar intercell wedge."""
from argparse import ArgumentParser
from itertools import product
import hashlib
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import consensus
from model.candidate_laws.interval_echo.candidate import native_law
from model.candidate_laws.interval_echo.clock import Registers
from model.executable_model import federation


HERE = Path(__file__).resolve().parent
INPUTS = ((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1))
DIMENSIONS = (11, 11, 3, 3, 3)
TOL = 1e-10


def source_map():
    # The final index labels an arbitrary input amplitude on the six-dimensional domain.
    source = np.zeros((*DIMENSIONS, len(INPUTS)))
    for col, loads in enumerate(INPUTS):
        source[(0, 0, *loads, col)] = 1
    return source


def local_ops():
    echo = native_law.echo()
    blank = np.diag(np.repeat(np.arange(11) == 0, 9).astype(float))
    pair = np.tile(np.array([a == b == 1 for a, b in native_law.SRC]), 11)
    strict_collar = [False, *(abs(contrast) == 2 for _, contrast in native_law.LABELS)]
    strict = pair & np.repeat(strict_collar, 9)
    write = echo @ blank
    hit = strict[:, None] * write
    # The prepared blank flag is read before the local echo. An occupied collar is left in place.
    probes = ((False, np.eye(99) - blank), (False, echo.T @ (write - hit)), (True, echo.T @ hit))
    return echo, pair, hit, probes


def act(amplitudes, op, seam):
    return consensus._apply(amplitudes, op, consensus.AXES[seam])


def branches(source, order, token, hit, probes):
    paths = [(source, {})]
    for seam in order:
        expanded = []
        for state, flags in paths:
            for ready, op in probes:
                branch = act(state, op, seam)
                if np.linalg.norm(branch) > 1e-13:
                    expanded.append((branch, {**flags, seam: ready}))
        paths = expanded
    outputs = {None: [], 0: [], 1: []}
    for state, flags in paths:
        ready = [seam for seam in order if flags[seam]]
        if ready:
            winner = token if len(ready) == 2 and token is not None else ready[0]
            result = act(state, hit, winner)
        else:
            winner, result = None, state
        outputs[winner].append(result.reshape(-1, len(INPUTS)))
    return outputs


def density_difference(a, b):
    """Restrict to occupied coordinates before forming the difference of two sums of outer products."""
    vectors = np.stack([*a, *b])
    supported = np.any(vectors != 0, axis=0)
    reduced = vectors[:, supported]
    x, y = reduced[:len(a)], reduced[len(a):]
    return x.T @ x.conj() - y.T @ y.conj()


def distance(a, b, ket_a=None, ket_b=None):
    # Accepted address is a classical direct-sum label in the output channel.
    size = 0.
    for which in a:
        if ket_a is None:
            delta = density_difference([k.ravel() for k in a[which]], [k.ravel() for k in b[which]])
            size += np.linalg.norm(delta)**2
        else:
            delta = density_difference([k @ ket_a for k in a[which]], [k @ ket_b for k in b[which]])
            size += np.sum(np.abs(np.linalg.eigvalsh(delta))) / 2
    return float(np.sqrt(size) if ket_a is None else size)


def with_transcript_distance(outputs, ket_a, ket_b):
    # Distinct outcomes of both sharp probes are classical orthogonal records.
    distance_sum = 0.
    for paths in outputs.values():
        for k in paths:
            delta = density_difference([k @ ket_a], [k @ ket_b])
            distance_sum += np.sum(np.abs(np.linalg.eigvalsh(delta))) / 2
    return float(distance_sum)


def check_outputs(outputs, probes):
    completeness = np.zeros((len(INPUTS), len(INPUTS)))
    material_energy = np.sum(np.indices((3, 3, 3))**2, axis=0)
    energy = np.broadcast_to(material_energy, DIMENSIONS).ravel()
    initial = np.sum(np.array(INPUTS)**2, axis=1)
    for winner, paths in outputs.items():
        for k in paths:
            completeness += k.T @ k
            if np.linalg.norm((energy[:, None] - initial + 2 * (winner is not None)) * k) > TOL:
                raise ArithmeticError('a transaction violates strict descent on the sharp input columns')
            state = k.reshape(*DIMENSIONS, len(INPUTS))
            for seam in (0, 1):
                if np.linalg.norm(act(state, probes[-1][1], seam)) > TOL:
                    raise ArithmeticError('a terminal branch still permits a strict write')
                if sum(np.linalg.norm(act(state, op, seam) - state) < TOL for _, op in probes[:-1]) != 1:
                    raise ArithmeticError('a terminal retry changes the material state')
    if np.linalg.norm(completeness - np.eye(len(INPUTS))) > TOL:
        raise ArithmeticError('the instrument omits a physical branch')


def layout(level):
    car = federation(level)
    middle, right = min(car.external)
    local = middle // 12
    left = min(p for p in car.neighbours[middle] if p // 12 == local)
    if right in car.neighbours[left]:
        raise ValueError('the declared wedge must have no outer link')
    rails = sorted(p for p in range(local*12, (local+1)*12) if p not in (left, middle))
    rail = next((a, b, c) for a, b, c in product(rails, repeat=3)
                if a != b and b != c and a != c and b in car.neighbours[a] and c in car.neighbours[a])
    observer = {local, right // 12}
    regs = Registers(car, seed=0)
    for a, b in ((left, middle), (middle, right)):
        regs.native(a, b, observer)
    regs.initialize_digits([rail])
    readings = []
    for value in (0, 1):
        regs.assign(rail, value, {local})
        observed = regs.read(rail, {local})
        if observed != value:
            raise ArithmeticError('priority differs from the locally retained rail')
        readings.append(observed)
    return dict(cells=car.cells, owners=sorted(observer), path=[left, middle, right],
                token_rails=list(rail), priority_readings=readings)


def run():
    source = source_map()
    echo, pair, hit, probes = local_ops()
    pair_read = pair[:, None] * echo
    false_accepts = []
    for seam in (0, 1):
        raw = act(source, pair_read, seam).reshape(-1, len(INPUTS))
        for col in range(3, 6):
            prob = float(raw[:, col] @ raw[:, col])
            if prob > TOL:
                false_accepts.append(dict(request=seam, initial=list(INPUTS[col]), probability=round(prob, 12)))
    boards = [layout(level) for level in (0, 1, 2)]
    d0, d1 = np.eye(len(INPUTS))[:2]
    plus, minus = (d0+d1)/np.sqrt(2), (d0-d1)/np.sqrt(2)
    rows = []
    for token in (None, 0, 1):
        ab, ba = [branches(source, order, token, hit, probes) for order in ((0, 1), (1, 0))]
        for outputs in (ab, ba):
            check_outputs(outputs, probes)
        error = distance(ab, ba)
        addressed = distance(ab, ba, d1, d1)
        phase = distance(ab, ab, plus, minus)
        if (error < TOL) != (token is not None) or abs(addressed - int(token is None)) > TOL:
            raise ArithmeticError('the protected endpoint differs from the declared admission test')
        transcript = with_transcript_distance(ab, plus, minus)
        rows.append(dict(priority=token, opposite_phases_with_transcript_distance=round(transcript, 12),
                         choi_frobenius_distance=round(error, 12),
                         middle_double_order_trace_distance=round(addressed, 12),
                         opposite_phases_endpoint_trace_distance=round(phase, 12)))
    return dict(protocol_sha256=hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest(),
                layouts=boards, input_basis=INPUTS, full_material_dimension=int(np.prod(DIMENSIONS)),
                legacy_ready_on_balanced_blank=false_accepts, priority_rows=rows,
                initial_opposite_phases_trace_distance=1,
                verified_tolerance=TOL)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, help='write a receipt to a fresh path')
    args = parser.parse_args()
    payload = json.dumps(run(), indent=2, allow_nan=False) + '\n'
    if args.record:
        with args.record.open('x', encoding='utf-8') as f:
            f.write(payload)
    elif json.loads(payload) != json.loads((HERE/'results.json').read_text()):
        raise SystemExit('consensus instrument disagrees with its deposited result')
    print(payload, end='')


if __name__ == '__main__':
    main()
