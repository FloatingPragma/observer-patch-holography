"""Release owner inventories on all exterior bridges; or open the whole dilute board."""
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier
from .measure import (PHASE_ERROR, RELEASE_SHOTS, SEED, SHOTS, TIMES, decode, detour, edges, fringe,
                      read_ports, source_only, survey)

HERE = Path(__file__).parent


def full_energy(board, doubles, blank):
    return 4*candidate.KAPPA*sum(6-board.degree[p] for p in doubles) + sum(
        blank[3*(2*(a in doubles)) + 2*(b in doubles),
                     3*(2*(a in doubles)) + 2*(b in doubles)] for a, b in board.seams)


def run():
    rng = np.random.default_rng(SEED+1)
    preparation_rng = np.random.default_rng(SEED+2)
    ex, co = candidate.native_law.local_term_parts()
    h = ex + candidate.KAPPA * co
    blank = h[:9, :9]
    idx = [2, 6]
    energies, frame = np.linalg.eigh(blank[np.ix_(idx, idx)])
    bright = frame[:, 0]
    branch_p = float(bright[0]**2)
    accepted = int(rng.binomial(RELEASE_SHOTS, branch_p))
    # A supplied exchange-parity filter prepares the even eigenvector from
    # 02 with success one half, on independently held copies of each bridge.
    calibration = fringe(float(energies[0]), rng)
    rows = []
    for level in range(3):
        board = carrier.federation(level)
        _, dual, gates, parked, _ = survey(board)
        distances = nx.single_source_shortest_path_length(dual, 0)
        reference = set(sorted(parked)[:2])
        adjacent = set(next((a, b) for a in sorted(parked) for b in board.neighbours[a]
                            if b > a and b in parked))
        if read_ports(reference, gates)[0] or read_ports(adjacent, gates)[0]:
            raise ValueError('Dilute reference has nonzero parity')
        for d in range(1, max(distances.values())+1):
            g = min(f for f, distance in distances.items() if distance == d)
            shortest = nx.shortest_path(dual, 0, g)
            count_reads = {}
            for name, route in (('least', edges(shortest)), ('detour', edges(detour(shortest, dual))),
                                ('sources_only', source_only(0, g, dual))):
                # Matched contacts are sharp blanks; every occupied symmetric
                # bridge closure is a separately measured 2x2 native block.
                energy = sum(np.vdot(bright, blank[np.ix_(idx, idx)] @ bright) for _ in route)
                counts = fringe(float(energy), rng)
                ready = 2 * len(TIMES) * SHOTS
                accepted_rate = branch_p**len(route)
                refused = int(preparation_rng.negative_binomial(ready, accepted_rate))
                count_reads[name] = dict(x_y=counts, measured_work=round(decode(counts), 8),
                                         launches_accepted=ready, launches_refused=refused,
                                         all_bridges_ready_probability=accepted_rate)
            necessary = source_only(0, g, dual)
            doubles = {gates[e][0] for e in necessary}
            if len(doubles) < 2:
                doubles.add(min(parked))
            odd, _ = read_ports(doubles, gates)
            if not {0, g} <= set(odd) or len(doubles) != 2:
                raise ValueError('Dilute source has the wrong source/inventory')
            upper = full_energy(board, doubles, blank)
            even_upper = full_energy(board, reference, blank)
            rows.append(dict(level=level, distance=d, coherent_work=count_reads,
                             dilute_upper=round(upper, 8),
                             dilute_even_original_trial=round(even_upper, 8),
                             dilute_even_original_ports=sorted(reference),
                             dilute_even_upper=round(full_energy(board, adjacent, blank), 8),
                             dilute_even_adjacent_ports=sorted(adjacent),
                             dilute_double_ports=sorted(doubles), dilute_odd=odd))
    return dict(local_eigenvalues=np.round(energies, 10).tolist(),
                projected_even_amplitudes=np.round(bright, 10).tolist(),
                filter_shots=RELEASE_SHOTS, filter_accepted=accepted,
                filter_rejected=RELEASE_SHOTS - accepted,
                calibration_x_y=calibration, calibration_work=decode(calibration),
                phase_half_width=PHASE_ERROR, rows=rows)


if __name__ == '__main__':
    result = run()
    data = json.dumps(result, indent=2, allow_nan=False) + '\n'
    filename = HERE/'released.json'
    if sys.argv[1:] == ['--write']:
        filename.write_text(data)
    elif sys.argv[1:] or filename.read_text() != data:
        raise SystemExit('Released report replay changed')
    print('Unguarded exterior bridge eigenvalues', result['local_eigenvalues'])
    print('Full-board dilute sharp-trial upper bound', max(row['dilute_upper'] for row in result['rows']))
