"""Exact spin-three-halves composite probe and its complete reducing block."""
from functools import lru_cache
from math import gcd, prod

import numpy as np
from scipy import linalg

from ..centered_particle_spectrum import source as parent


def block(columns=None):
    if columns is not None and (type(columns) is not tuple or
                               any(type(x) is not int for x in columns)):
        raise ValueError('exact tuple of projector columns required')
    return _block(columns)


@lru_cache(maxsize=None)
def _block(columns=None):
    up, down = 5, 2
    q, representatives, sizes = parent.invariant_orbits(up, down, 'local')
    n = q.shape[1]
    eye = np.eye(n, dtype=np.int64)
    projectors = []
    denominator = 1
    for name, roots in (('spin_plus', (5, 12)), ('eta_minus', (3, 8))):
        ladder = parent.ladder(up, down, name)
        ell = parent.restriction(ladder.T @ ladder, q, representatives).toarray()
        p = eye.copy()
        for root in roots:
            p = parent.int_product(p, root * eye - ell)
        d = prod(roots)
        if np.any(parent.int_product(ell, p)):
            raise ArithmeticError('ladder kernel not exhausted')
        if not np.array_equal(parent.int_product(p, p), d * p):
            raise ArithmeticError('ladder polynomial not a projector')
        projectors.append(p)
        denominator *= d
    p = parent.int_product(*projectors)
    if not np.array_equal(p, parent.int_product(projectors[1], projectors[0])):
        raise ArithmeticError('spin and eta projectors do not commute')
    common = denominator
    for value in p.flat:
        common = gcd(common, int(value))
    p //= common
    denominator //= common
    if not np.array_equal(parent.int_product(p, p), denominator * p):
        raise ArithmeticError('joint projector not idempotent')
    if not np.array_equal(p.T * sizes[None, :], sizes[:, None] * p):
        raise ArithmeticError('joint projector not metric self-adjoint')
    trace = int(np.trace(p))
    if trace % denominator:
        raise ArithmeticError('nonintegral projector rank')
    rank = trace // denominator
    if columns is None:
        _, _, pivots = linalg.qr(np.sqrt(sizes)[:, None] * p.astype(float),
                                 mode='economic', pivoting=True)
        columns = tuple(int(i) for i in pivots[:rank])
    if (len(columns) != rank or len(set(columns)) != rank or
            any(i < 0 or i >= n for i in columns)):
        raise ValueError('complete distinct projector columns required')
    x = p[:, columns].copy()
    for j in range(rank):
        divisor = 0
        for value in x[:, j]:
            divisor = gcd(divisor, int(value))
        if not divisor:
            raise ArithmeticError('zero proposed column')
        x[:, j] //= divisor
    if parent.modular_rank(x) != rank:
        raise ArithmeticError('proposed columns do not span exact projector')
    full = np.asarray(q @ x)
    for name in ('spin_plus', 'eta_minus'):
        if np.any(parent.ladder(up, down, name) @ full):
            raise ArithmeticError('joint kernel membership failed')
    t = parent.restriction(parent.centered(up, down), q, representatives).toarray()
    if not np.array_equal(parent.int_product(t, p), parent.int_product(p, t)):
        raise ArithmeticError('source leaves composite block')
    gram = parent.int_product(x.T, sizes[:, None] * x)
    numerator = parent.int_product(x.T, sizes[:, None] * parent.int_product(t, x))
    if not np.array_equal(numerator, numerator.T):
        raise ArithmeticError('source compression not symmetric')
    return dict(up=up, down=down, q=q, representatives=representatives, sizes=sizes,
                basis=x, full_basis=full, gram=gram, numerator=numerator,
                columns=np.asarray(columns, dtype=np.int64), dimension=rank,
                orbit_dimension=n, projector=p, denominator=denominator,
                coefficient_source=t)


@lru_cache(maxsize=1)
def probe():
    """A_0=(sum over neighbors y of S_y^+) c_0,down, norm at most three."""
    a = parent.ground.band_data()[0]
    neighbors = tuple(int(y) for y in np.flatnonzero(a[0]))
    if neighbors != (1, 2, 4):
        raise ArithmeticError('source neighbor support changed')
    remove = parent.car(4, 4, 1)
    return sum((parent.car(4, 2, 2*y, True) @
                parent.car(4, 3, 2*y+1) @ remove for y in neighbors))


def transition(ground, composite):
    image = probe() @ ground['full_basis']
    coefficients = image[composite['representatives']]
    if not np.array_equal(composite['q'] @ coefficients, image):
        raise ArithmeticError('probe leaves origin stabilizer fixed space')
    if not np.array_equal(parent.int_product(composite['projector'], coefficients),
                          composite['denominator'] * coefficients):
        raise ArithmeticError('probe leaves spin-three-halves eta-half block')
    return parent.int_product(composite['basis'].T,
                              composite['sizes'][:, None] * coefficients)


def selection_certificate(ground, composite):
    transition(ground, composite)
    image = probe() @ ground['full_basis']
    if not np.any(image):
        raise ArithmeticError('composite probe identically zero on ground space')
    for name in ('spin_plus', 'eta_minus'):
        if np.any(parent.ladder(5, 2, name) @ image):
            raise ArithmeticError('highest spin / lowest eta selection failed')
    k = parent.centered(5, 2)
    for name, target in (('spin_plus', (6, 1)), ('eta_minus', (4, 1))):
        op = parent.ladder(5, 2, name)
        if (parent.centered(*target) @ op - op @ k).nnz:
            raise ArithmeticError('source spin/eta conservation failed')
    return {'source_number': 7, 'source_hidden_spin': '3/2', 'source_eta': '1/2',
            'reference_number': 8, 'reference_hidden_spin': '0', 'reference_eta': '0',
            'probe_neighbors': [1, 2, 4], 'probe_operator_norm_upper': '3',
            'complete_reducing_block_dimension': composite['dimension'],
            'origin_stabilizer_dimension': composite['orbit_dimension'],
            'all_linear_CAR_responses_orthogonal': True,
            'physical_Lorentz_spin_identified': False}
