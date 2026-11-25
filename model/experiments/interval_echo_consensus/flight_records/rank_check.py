"""Independent finite-field rank lower bounds for the 20-cell derivative tests."""

from collections import defaultdict

from model.executable_model import CELL_SEAMS, federation


def mod_rank(rows, pivots):
    prime = 101
    for source in rows:
        remaining = {k: v % prime for k, v in source.items() if v % prime}
        while remaining:
            i = min(remaining)
            coefficient = remaining[i]
            if i not in pivots:
                inverse = pow(coefficient, -1, prime)
                pivots[i] = {k: inverse * v % prime for k, v in remaining.items()}
                break
            for k, v in pivots[i].items():
                new = (remaining.get(k, 0) - coefficient * v) % prime
                if new:
                    remaining[k] = new
                else:
                    remaining.pop(k, None)
    return len(pivots)


def check():
    board = federation(0)
    neighbors = [set() for _ in range(board.n)]
    for x, y in board.seams:
        neighbors[x].add(y)
        neighbors[y].add(x)
    winners = []
    for x, ys in enumerate(neighbors):

        def rank(y):
            if x // 12 != y // 12:
                return 0
            return 1 + CELL_SEAMS.index(tuple(sorted((x % 12, y % 12))))

        y = min(ys, key=rank)
        winners.append(tuple(sorted((x, y))))
    basis = [
        (x, y) for x, edge in enumerate(winners) for y in edge if winners[y] == edge
    ]
    lap = [{x: len(ys)} | {y: -1 for y in ys} for x, ys in enumerate(neighbors)]
    lap2 = []
    for row in lap:
        square = defaultdict(int)
        for k, v in row.items():
            for j, u in lap[k].items():
                square[j] += u * v
        lap2.append(square)

    first, second = defaultdict(dict), defaultdict(dict)
    for col, (i, j) in enumerate(basis):

        def add(rows, x, y, coefficient):
            if winners[x] != winners[y]:
                rows[x, y][col] = rows[x, y].get(col, 0) + coefficient

        for a, v in lap[i].items():
            add(first, a, j, v)
            for b, u in lap[j].items():
                add(second, a, b, -2 * v * u)
        for b, u in lap[j].items():
            add(first, i, b, -u)
        for a, v in lap2[i].items():
            add(second, a, j, v)
        for b, u in lap2[j].items():
            add(second, i, b, u)
    pivots = {}
    dims = [
        len(basis),
        len(basis) - mod_rank(first.values(), pivots),
        len(basis) - mod_rank(second.values(), pivots),
    ]
    if dims != [316, 20, 1]:
        raise ValueError(f"unexpected dimensions modulo 101: {dims}")
    return dims


if __name__ == "__main__":
    print("20-cell nullities modulo 101:", check())
