"""Numerical cross-check in a different Clifford representation.

The exact positivity proof belongs to gravity/exact_witness.py. This script
independently evaluates its projector invariant without a polynomial sign(H).
It is not an entropy/Einstein coefficient calculation.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1])
    I = np.eye(4)
    # Tensor factors are swapped relative to the producer's representation.
    alpha = [np.kron(X, p) for p in (X, Y, Z)]
    beta = np.kron(Z, np.eye(2))
    mass = math.sqrt(399)/20*I - 1j*beta/20
    f = [(I-1j*a)/math.sqrt(2) for a in alpha]
    u = f[2]@f[1]@f[0]@mass
    du = f[2]@f[1]@(-I-1j*alpha[0])/math.sqrt(2)@mass
    h = .5j*(u-u.conj().T)
    dh = .5j*(du-du.conj().T)
    k = h@dh+dh@h
    energy, vectors = np.linalg.eigh(h)
    positive = vectors[:, energy > 0]@vectors[:, energy > 0].conj().T
    negative = I-positive
    invariant = 4*float(np.linalg.norm(positive@k@negative, 'fro')**2)
    expected = 80199/32240000 + 79401*math.sqrt(319201)/10291040240000
    assert abs(invariant-expected) < 1e-14
    vertex = vectors.conj().T@dh@vectors
    c2 = sum(.5*abs((eb+ea)*vertex[ib, ia])**2/(eb-ea)**3
             for ib, eb in enumerate(energy) if eb > 0
             for ia, ea in enumerate(energy) if ea < 0)
    assert c2 >= invariant/64 > 0
    paths = [HERE/'gravity/EXACT_WITNESS.json', HERE/'gravity/exact_witness.py', Path(__file__)]
    result = {'scope': 'Independent single-momentum positive-interband diagnostic',
              'representation': 'tensor factors swapped; direct spectral projectors',
              'interband_invariant': invariant, 'exact_expression_decimal': expected,
              'c2_local': c2, 'expression_lower_bound_decimal': expected/64,
              'einstein_area_coefficient': False, 'finite_mass_mu': math.asin(.05),
              'status': 'positive_susceptibility_scope_check_pass',
              'pins': {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (HERE/'GRAVITY_WITNESS_REVIEW.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('interband_invariant', 'c2_local', 'status')}))


if __name__ == '__main__':
    main()
