"""Exact rational timelike tomography and finite-radius error propagation."""
from fractions import Fraction as F

PAIRS = ((0, 0), (0, 1), (0, 2), (0, 3), (1, 1),
         (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))
FRAME = ((F(1), F(0), F(0), F(0)),
         (F(5, 3), F(4, 3), F(0), F(0)),
         (F(5, 3), F(-4, 3), F(0), F(0)),
         (F(5, 3), F(0), F(4, 3), F(0)),
         (F(5, 3), F(0), F(-4, 3), F(0)),
         (F(5, 3), F(0), F(0), F(4, 3)),
         (F(5, 3), F(0), F(0), F(-4, 3)),
         (F(3), F(2), F(2), F(0)),
         (F(3), F(2), F(0), F(2)),
         (F(3), F(0), F(2), F(2)))


def design(frame=FRAME):
    return tuple(tuple(u[i]*u[j]*(1 if i == j else 2)
                       for i, j in PAIRS) for u in frame)


def decode(q):
    if len(q) != 10:
        raise ValueError('ten matched timelike charges required')
    return (q[0], F(9, 80)*(q[1]-q[2]), F(9, 80)*(q[3]-q[4]),
            F(9, 80)*(q[5]-q[6]),
            F(9, 32)*(q[1]+q[2])-F(25, 16)*q[0],
            q[7]/8+F(7, 16)*q[0]-F(99, 320)*(q[1]+q[3])+F(9, 320)*(q[2]+q[4]),
            q[8]/8+F(7, 16)*q[0]-F(99, 320)*(q[1]+q[5])+F(9, 320)*(q[2]+q[6]),
            F(9, 32)*(q[3]+q[4])-F(25, 16)*q[0],
            q[9]/8+F(7, 16)*q[0]-F(99, 320)*(q[3]+q[5])+F(9, 320)*(q[4]+q[6]),
            F(9, 32)*(q[5]+q[6])-F(25, 16)*q[0])


def apply(matrix, vector):
    return tuple(sum(a*b for a, b in zip(row, vector, strict=True)) for row in matrix)


def basis(n, j):
    return tuple(F(i == j) for i in range(n))


def decoder_matrix():
    columns = tuple(decode(basis(10, j)) for j in range(10))
    return tuple(tuple(columns[j][i] for j in range(10)) for i in range(10))


def inverse(matrix):
    """Independent exact Gauss-Jordan inversion; no numerical rank threshold."""
    n = len(matrix)
    rows = [list(map(F, row))+list(basis(n, i)) for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            raise ValueError('tomography frame is singular')
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        rows[j] = [entry/divisor for entry in rows[j]]
        for i in range(n):
            if i != j:
                multiplier = rows[i][j]
                rows[i] = [x-multiplier*y for x, y in zip(rows[i], rows[j], strict=True)]
    return tuple(tuple(row[n:]) for row in rows)


def normalized_error(eta, bulk_over_l4, area_over_l4, stationarity_over_l4):
    """Return pi times the component bound, without approximating pi.

    Inputs are nonnegative errors divided by the SAME proper ell^4.
    """
    eta = F(eta)
    errors = tuple(map(F, (bulk_over_l4, area_over_l4, stationarity_over_l4)))
    if eta <= 0 or any(value < 0 for value in errors):
        raise ValueError('positive area coefficient and nonnegative errors required')
    eb, ea, es = errors
    return F(17, 8)*F(15, 4)*(eb+eta*ea+es)/eta


def finite_profile(n):
    """Supplied ell=1/(n+1), errors=ell^5; no inhabited source tower."""
    if not isinstance(n, int) or n < 0:
        raise ValueError('nonnegative integer index required')
    ell = F(1, n+1)
    return {'radius': str(ell), 'area_coefficient': '2',
            'normalized_error_each': str(ell),
            'pi_times_component_bound': str(normalized_error(2, ell, ell, ell)),
            'source_premise_instance': False}


def certificate():
    forward, backward = design(), decoder_matrix()
    if inverse(forward) != backward:
        raise ValueError('explicit decoder disagrees with independent inversion')
    for u in FRAME:
        if -u[0]**2+sum(v*v for v in u[1:]) != -1 or u[0] <= 0:
            raise ValueError('probe is not future unit timelike')
    row_norms = tuple(sum(abs(entry) for entry in row) for row in backward)
    if max(row_norms) != F(17, 8):
        raise ValueError('sharp inverse norm changed')
    missing = []
    for j in range(10):
        witness = decode(basis(10, j))
        if apply(forward, witness) != basis(10, j):
            raise ValueError('missing-probe null witness failed')
        missing.append({'deleted_probe': j, 'nonzero_tensor': list(map(str, witness))})
    sharp = tuple(F(-1) if i == 0 else F(1) if i in (1, 2) else F(0) for i in range(10))
    if decode(sharp)[4] != F(17, 8):
        raise ValueError('sharpness witness failed')
    return {
        'signature': [-1, 1, 1, 1],
        'coordinate_pairs': [list(pair) for pair in PAIRS],
        'timelike_frame': [list(map(str, row)) for row in FRAME],
        'design': [list(map(str, row)) for row in forward],
        'inverse': [list(map(str, row)) for row in backward],
        'inverse_row_l1_norms': list(map(str, row_norms)),
        'sharp_coordinate_infinity_norm': '17/8',
        'sharp_charges': list(map(str, sharp)),
        'deleted_probe_witnesses': missing,
        'normalizations': {
            'radial_kernel_integral': '4*pi*ell^4/15',
            'bulk_entropy_coefficient': '8*pi^2*ell^4/15',
            'fixed_volume_area_coefficient': '-4*pi*ell^4/15',
            'einstein_coupling': '2*pi/eta_star',
            'Newton_constant_in_hbar_c_equal_one_units': '1/(4*eta_star)',
            'component_error_bound': '255*(e_bulk+eta*e_area+e_stationarity)/(32*pi*eta*ell^4)',
        },
        'finite_algebraic_profiles': [finite_profile(n) for n in (0, 1, 3, 7, 15, 31)],
        'profile_checks_prove_a_source_refinement': False,
    }
