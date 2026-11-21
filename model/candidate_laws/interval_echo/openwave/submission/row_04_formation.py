"""Can ordinary free motion of initially separated units prepare the released pair?"""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from row_04_check import evolve, screen


def encounter(board):
    _, index, exchange, adjacent = screen(board)
    generator = exchange - 16 * sparse.diags(adjacent)
    if (generator != generator.T).nnz:
        raise ArithmeticError("Launch and herald amplitudes require a symmetric generator")
    rows = []
    for cell, target in ((0, (0, 1)), (13, (156, 157))):
        states = list(itertools.combinations(range(12 * cell, 12 * cell + 12), 2))
        separated = [s for s in states if adjacent[index[s]] == 0]
        if len(separated) != 36 or adjacent[index[target]] != 1:
            raise ValueError("Unexpected preparation domain")
        # exp(-itH) is complex symmetric. Each entry equals the amplitude
        # from its separated launch to the target seam at the same time.
        terminal = next(evolve(generator, index[target], [1]))
        probabilities = [float(abs(terminal[index[s]]) ** 2) for s in separated]
        example = int(np.argmax(probabilities))
        forward = next(evolve(generator, index[separated[example]], [1]))
        if abs(forward[index[target]] - terminal[index[separated[example]]]) > 1e-12:
            raise ArithmeticError("Directly launched separated pair disagrees with reciprocity")
        p = float(np.mean(probabilities))
        rows.append(dict(cell=cell, target=list(target), exposure=1,
                         conditional_on_herald="exact target pair ket, up to phase",
                         uniformly_mixed_separated_starts=36, herald=round(p, 12),
                         quiet=round(1 - p, 12),
                         preparations=[dict(ports=list(s), herald=round(prob, 12),
                                            quiet=round(1 - prob, 12))
                                       for s, prob in zip(separated, probabilities, strict=True)]))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--board", type=Path, required=True)
    ap.add_argument("--record", action="store_true")
    args = ap.parse_args()
    report = encounter(args.board)
    receipt = Path(__file__).resolve().with_name("row_04_formation.json")
    if args.record:
        with receipt.open("x") as f:
            f.write(json.dumps(report, indent=2) + "\n")
    if report != json.loads(receipt.read_text()):
        raise SystemExit("The separated-start preparation probabilities changed")
    for row in report:
        print(row["cell"], row["target"], row["herald"], row["quiet"])
