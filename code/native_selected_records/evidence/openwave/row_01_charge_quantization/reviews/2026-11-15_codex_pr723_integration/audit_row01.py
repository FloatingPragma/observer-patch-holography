"""Root same-group audit: custody, replay, and independent proof challenges.

The generalized path checks are independent finite challenges to the written
proof, not a substitute for that proof and not literal producer-native calls.
"""
from collections import deque
from hashlib import sha256
import json
from pathlib import Path
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1] / 'work/2026-11-15_codex_pr723_unit_selection'


def repair(x, a, b, upper_on_a):
    result = x[:]
    total = x[a] + x[b]
    lower, upper = total // 2, total - total // 2
    result[a], result[b] = (upper, lower) if upper_on_a else (lower, upper)
    return result


def decrease_path(x, adjacency):
    maximum = max(x)
    start = x.index(maximum)
    paths = deque([[start]])
    seen = {start}
    while paths:
        path = paths.popleft()
        if x[path[-1]] <= maximum - 2:
            break
        for v in sorted(adjacency[path[-1]]):
            if v not in seen:
                seen.add(v)
                paths.append(path + [v])
    else:
        raise AssertionError('No balancing path in an unbalanced connected graph')
    y, steps = x[:], 0
    for a, b in zip(path[:-2], path[1:-1]):
        assert y[a] == maximum
        if y[b] == maximum - 1:
            y = repair(y, a, b, False)
            steps += 1
        else:
            assert y[b] == maximum
    a, b = path[-2:]
    assert y[a] - y[b] >= 2
    y = repair(y, a, b, False)
    steps += 1
    assert steps <= len(x) - 1
    assert sum(y) == sum(x)
    assert sum(v*v for v in x) - sum(v*v for v in y) >= 2
    return y, steps


def challenge_proof():
    rng = random.Random(723991)
    cases = events = 0
    for n in range(2, 10):
        for trial in range(32):
            # A random labelled tree and extra edges exercise plateaus and
            # connected graphs unlike the particular twelve-port experiment.
            adjacency = [set() for _ in range(n)]
            for b in range(1, n):
                a = rng.randrange(b)
                adjacency[a].add(b)
                adjacency[b].add(a)
            for a in range(n):
                for b in range(a + 1, n):
                    if rng.randrange(5) == 0:
                        adjacency[a].add(b)
                        adjacency[b].add(a)
            x = [rng.randrange(-9, 10) for _ in range(n)]
            original = x[:]
            base, rem = divmod(sum(x), n)
            minimum = (n-rem)*base*base + rem*(base+1)*(base+1)
            decreases = (sum(v*v for v in x) - minimum) // 2
            for _ in range(decreases + 1):
                if max(x) - min(x) <= 1:
                    break
                x, used = decrease_path(x, adjacency)
                events += used
            assert sorted(x) == [base]*(n-rem) + [base+1]*rem
            assert min(original) <= min(x) <= max(x) <= max(original)
            for a in range(n):
                for b in adjacency[a]:
                    for coin in (False, True):
                        y = repair(x, a, b, coin)
                        assert max(y) - min(y) <= 1
            cases += 1
    assert repair([2, 2], 0, 1, False) == [2, 2]
    return {'signed_graph_cases': cases, 'independent_floor_ceiling_events': events,
            'uniform_two_counterexample_retained': True}


def main():
    manifest = PACKET / 'MANIFEST.sha256'
    files = []
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        assert sha256((PACKET/name).read_bytes()).hexdigest() == digest, name
        files.append(name)
    start = time.monotonic()
    replay = subprocess.run([sys.executable, '-B', str(PACKET/'headline_strict.py')],
                            capture_output=True, text=True, timeout=60, check=True)
    (HERE/'ROW01_REPLAY.log').write_text(replay.stdout + replay.stderr)
    result = {'same_group': True, 'external_G6': False,
              'bound_manifest_sha256': sha256(manifest.read_bytes()).hexdigest(),
              'custody_files': len(files), 'headline_seconds': time.monotonic()-start,
              'proof_challenges': challenge_proof(), 'passed': True,
              'boundary': 'Unit-selection theorem on declared integer domain; electric/topological identification remains open.'}
    (HERE/'ROW01_AUDIT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
