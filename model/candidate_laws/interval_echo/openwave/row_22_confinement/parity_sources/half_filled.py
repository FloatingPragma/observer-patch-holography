"""State-independent probability and work bound for the two local parity pins."""
from fractions import Fraction
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier
from .measure import survey


def character(sites, support):
    """Mean product of the selected signs under uniform exact half filling."""
    if support % 2:
        return Fraction(0)
    half = support // 2
    return Fraction((-1)**half * comb(sites//2, half), comb(sites, support))


def both_odd(sites, a, b, difference):
    return (1-character(sites, a)-character(sites, b)+character(sites, difference))/4


def describe(level, spectral_width):
    board = carrier.federation(level)
    rings, _, gates, _, _ = survey(board)
    supports = [set() for _ in rings]
    for (a, b), ports in gates.items():
        supports[a].update(ports)
        supports[b].update(ports)
    classes = {}
    for x, y in combinations(supports, 2):
        a, b = sorted((len(x), len(y)))
        sig = (a, b, len(x ^ y))
        classes[sig] = classes.get(sig, 0) + 1
    n = board.n
    width = max((a+b+diff)//2 for a, b, diff in classes)
    probabilities = {sig: both_odd(n, *sig) for sig in classes}
    lower = Fraction(1, 4) * (1-Fraction(3, n-1))
    if any(diff == 0 or diff > n//2 for _, _, diff in classes) or min(probabilities.values()) < lower:
        raise ValueError('Parity family is outside the certificate')
    degree = int(max(board.degree))
    service_count = max(int(d) + (6-int(d)) for d in board.degree)
    commutator_bound = 2*width*service_count*spectral_width
    return dict(ports=n, total_doubles=n//2, faces=len(rings), pairs=sum(classes.values()),
                classes=[dict(support_sizes=list(key), pairs=classes[key],
                              uniform_both_odd_exact=str(probabilities[key]))
                         for key in sorted(classes)],
                largest_union=width, degree=degree, native_services_per_port=service_count,
                private_service_width=4*candidate.KAPPA, smallest_exact_probability=str(min(probabilities.values())),
                common_probability_lower=str(lower), commutator_upper=commutator_bound,
                source_work_upper=commutator_bound/np.sqrt(float(lower)))


def run():
    ex, co = candidate.native_law.local_term_parts()
    binary = [0, 2, 6, 8]
    h = (ex + candidate.KAPPA*co)[np.ix_(binary, binary)]
    width = float(np.ptp(np.linalg.eigvalsh(h)))
    if not np.isclose(width, 33, atol=1e-12) or 4 * candidate.KAPPA > width:
        raise ValueError('Unexpected local spectral diameter')
    return dict(native_binary_seam_width=round(width, 10),
                levels=[describe(level, round(width, 10)) for level in range(3)])


if __name__ == '__main__':
    result = run()
    data = json.dumps(result, indent=2, allow_nan=False) + '\n'
    filename = Path(__file__).with_name('half_filled.json')
    if sys.argv[1:] == ['--write']:
        filename.write_text(data)
    elif sys.argv[1:] or filename.read_text() != data:
        raise SystemExit('Half-filled certificate replay changed')
    for row in result['levels']:
        print(row['ports'], row['pairs'], row['largest_union'], row['source_work_upper'])
