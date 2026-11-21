"""Try the stationary hand at the pilot detector. No exposure enters its clock."""
from collections import deque

import numpy as np
from scipy import sparse

from model.candidate_laws.interval_echo import candidate as echo, clock
from model.executable_model import federation

from .channel import SPEC, encode, readout, require


def banks(board, excluded):
    path = []
    queue = deque([0])
    seen = {0, excluded}
    while queue and len(path) < 3:
        cell = queue.popleft()
        for p in range(12 * cell, 12 * cell + 12):
            for q in board.neighbours[p]:
                neighbor = q // 12
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
                    path.append(neighbor)
    require(len(path) >= 3, "station needs three local counter banks")
    return path[:3]


def hardware(board, station, bank_cells, deadline, seed):
    owner = {0, *bank_cells}
    count = [d for c in bank_cells for d in clock._digits(c)]
    word = clock._digits(bank_cells[0], clock.EXTRA_DIGITS[:1])[0]
    reg = clock.Registers(board, seed)
    reg.prepare_front(0)
    reg.initialize_digits(count + [word])
    start = reg.reading(0, count, owner)
    near, far = station["contact"]
    require(start == 0 and reg.loads[near] == 2 and reg.loads[far] == 0, "initial clock contact")
    events = []
    status = "timeout"
    try:
        for _ in range(8):
            primary = reg.fair(0, owner)
            reg.publish_tick(word, primary, owner)
            tick = reg.post_tick(word, 0, count, owner)
            elapsed = reg.reading(0, count, owner)
            events.append(dict(tick=tick, elapsed=elapsed))
            if elapsed == deadline:
                status = "detector_requested_joint_read_unavailable"
                break
    except OverflowError:
        status = "overflow"
    except (ValueError, ArithmeticError, AssertionError) as error:
        status = f"unavailable: {error}"
    return dict(deadline=deadline, start=start, public_digit=list(word), events=events,
                status=status, committed_answer=None, gates=dict(reg.gates))


def probe():
    station = SPEC["stations"][0]
    board = federation(station["level"])
    near, far = station["contact"]
    bank_cells = banks(board, far // 12)
    trials = [hardware(board, station, bank_cells, ticks, seed)
              for ticks, seed in zip((1, 2), SPEC["seeds"]["clock"], strict=True)]
    # At release the five-position hand is a double on port 1. The direct
    # link joins it to a single for bit zero and a vacancy for bit one.
    law = echo.native_law
    h0, h1 = law.local_term_parts()
    local = h0 + echo.KAPPA * h1
    leak = [float(local[law.SRC.index((n, 2)), law.SRC.index((2, n))]) for n in (1, 0)]
    require(all(abs(v + .5) < 1e-12 for v in leak), "reassess common-sector leakage")
    # Cancel just the exchange on this link as a protection diagnostic.
    # Rank-one difference vv*/2; both endpoints and every other edge survive.
    v = sparse.csr_matrix(([1., -1.], ([near, far], [0, 0])), shape=(board.n, 1))
    h = echo.generator(board)
    source = encode(board, station)
    s = SPEC["control_exposure"]
    free, _ = readout(h, source, s)
    protected, _ = readout(h - (v @ v.T) * .5, source, s)
    delta = float(np.max(abs(free - protected)))
    require(delta >= SPEC["minimum_law_shift"], "reassess protection effect on detector")
    return dict(verdict="JOIN_INCOMPLETE", cells=board.cells, bank_cells=bank_cells,
                contact=[near, far], off_block_by_bit=leak, isolated_clock_trials=trials,
                protection_test_exposure=s, protected_one_load_law=protected.tolist(),
                protection_max_probability_change=delta, joint_message_trials=0)
