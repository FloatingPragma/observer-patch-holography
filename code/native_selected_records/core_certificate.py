"""Exact cell polynomial establishing the core-gap input without an eigensolver."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE/'evidence/openwave/row_01_charge_quantization/work/2026-11-15_codex_pr723_unit_selection/seams.json'
EXPECTED = '7ce94f057fa081c9902caed98151e7a2cfeb7735e2dc5c8dd79f9052bc30fdf2'


def product(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(12))
             for j in range(12)] for i in range(12)]


def check():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED
    edges = json.loads(SOURCE.read_text())
    lap = [[0]*12 for _ in range(12)]
    reached = {0}
    for a,b in edges:
        lap[a][a] += 1
        lap[b][b] += 1
        lap[a][b] -= 1
        lap[b][a] -= 1
    for _ in range(12):
        for a,b in edges:
            if a in reached or b in reached:
                reached.update((a,b))
    assert len(reached) == 12
    square = product(lap,lap)
    shift = [[lap[i][j]-6*int(i==j) for j in range(12)] for i in range(12)]
    quadratic = [[square[i][j]-10*lap[i][j]+20*int(i==j)
                  for j in range(12)] for i in range(12)]
    residual = product(product(lap,shift),quadratic)
    assert all(value==0 for row in residual for value in row)
    # Symmetry, connectedness and the annihilating polynomial put every
    # nonconstant eigenvalue in {6,5-sqrt(5),5+sqrt(5)}, all strictly above2.
    return {'passed':True, 'exact_integer_polynomial_zero':True,
            'connected_ports':12, 'seams':len(edges),
            'polynomial':'L(L-6I)(L^2-10L+20I)=0',
            'nonzero_gap_strictly_greater_than':2,
            'scope':'Finite source graph certificate; no physical energy claim'}


if __name__ == '__main__':
    print(json.dumps(check(),indent=2))
