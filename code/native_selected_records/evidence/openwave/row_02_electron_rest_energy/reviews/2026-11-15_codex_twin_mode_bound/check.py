"""Independent finite checks of the written twin-mode variational theorem."""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np

CELL = ((0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),
        (2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),
        (4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),
        (7,8),(7,10),(8,9),(10,11))


def run():
    cases = []
    for blocks in (1, 2, 3):
        n = 12*blocks
        mu = min(F(2), F(4, blocks*blocks))
        w, v = 2*n, 1
        lam = w+2*v
        subtraction = F(2*w*(2*n+3)*(n-1), n)
        c = 1
        while c*mu <= lam+subtraction:
            c *= 2
        edges = [(12*k+a, 12*k+b, c) for k in range(blocks) for a,b in CELL]
        edges += [(12*k+a, 12*(k+1)+a, c)
                  for k in range(blocks-1) for a in range(12)]
        edges += [(0,n,w), (0,n+1,w), (n,n+1,v)]
        size = n+2
        lap = [[0]*size for _ in range(size)]
        counts = [[0]*size for _ in range(size)]
        evaluations = 0
        for a,b,weight in edges:
            lap[a][a] += weight
            lap[b][b] += weight
            lap[a][b] -= weight
            lap[b][a] -= weight
            for p in range(size):
                q = [int(i==p) for i in range(size)]
                for coin in (0,1):
                    total = q[a]+q[b]
                    y = q.copy()
                    y[a],y[b] = ((total//2,total-total//2) if coin==0
                                 else (total-total//2,total//2))
                    destination = y.index(1)
                    counts[destination][p] += weight
                    evaluations += 1
        weight_sum = sum(weight for _,_,weight in edges)
        f = [0]*n+[1,-1]
        for i in range(size):
            assert sum(lap[i][j]*f[j] for j in range(size)) == lam*f[i]
            for j in range(size):
                assert counts[i][j] == 2*weight_sum*int(i==j)-lap[i][j]
        bound = min(F(w+3*v), c*mu-subtraction)
        assert bound > lam
        spectrum = np.linalg.eigvalsh(np.array(lap, dtype=float))
        assert abs(spectrum[0]) < 1e-7
        assert abs(spectrum[1]-lam) < 1e-7
        assert spectrum[2] > float(bound)-1e-7
        cases.append(dict(blocks=blocks, core_ports=n, core_multiplier=c,
                          leaf_weight=w, leaf_bridge=v, compact_eigenvalue=lam,
                          other_mode_lower_bound=str(bound),
                          numerical_next_eigenvalue=float(spectrum[2]),
                          transfer_eigenvalue=str(1-F(lam,2*weight_sum)),
                          independent_endpoint_evaluations=evaluations))
    return dict(passed=True, same_group=True, auxiliary_theorem_only=True,
                physical_energy_demonstrated=False, cases=cases)


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
