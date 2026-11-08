"""A cutoff-uniform resolvent for the implemented, distinct U(1) rotor source."""
from fractions import Fraction as F

from model.interacting_extension.chiral_lattice.rotor_response import build_flux_ring
from .action import interval_packet


def integer(value, name, minimum=None):
    if type(value) is not int or (minimum is not None and value < minimum):
        raise ValueError('exact integer '+name+' outside its domain')
    return value


def diagonal(q, position, flux, sites=4):
    return F(position*(flux+q)**2+(sites-position)*flux**2, 2)


def add(out, label, value):
    value += out.get(label, F(0))
    if value:
        out[label] = value
    else:
        out.pop(label, None)


def apply_v(q, vector, sites=4, cutoff=None):
    out = {}
    for (x, ell), value in vector.items():
        targets = (((x+1)%sites, ell+q if x == sites-1 else ell, F(-1)),
                   ((x-1)%sites, ell-q if x == 0 else ell, F(-1)),
                   (x, ell+1, F(-1, 2)), (x, ell-1, F(-1, 2)))
        for y, flux, coefficient in targets:
            if cutoff is None or abs(flux) <= cutoff:
                add(out, (y, flux), coefficient*value)
    return out


def cutoff_error(q, cutoff, mu=F(8), initial_flux=0, sites=4):
    integer(q, 'charge')
    integer(cutoff, 'cutoff', 0)
    integer(initial_flux, 'initial flux radius', 0)
    integer(sites, 'sites', 3)
    if type(mu) not in (int, F):
        raise TypeError('exact rational spectral resolution required')
    mu = F(mu)
    if mu <= 3 or initial_flux > cutoff:
        raise ValueError('spectral margin and supported input required')
    step = max(1, abs(q))
    count = (cutoff-initial_flux+step)//step
    vector_error = 2/(mu-3)*(3/mu)**count
    return_error = 2/(mu-3)*(3/mu)**(2*count)
    return {'maximum_flux_step':step, 'minimum_exit_insertions':count,
            'minimum_exit_and_return_insertions':2*count,
            'full_to_cutoff_resolvent_on_supported_inputs':str(vector_error),
            'local_compressed_resolvent_error':str(return_error),
            'uniform_in_all_larger_cutoffs':True}


def local_resolvent(q, mu=F(8), degree=24, sites=4, cutoff=None):
    integer(q, 'charge')
    integer(degree, 'Neumann degree', 0)
    integer(sites, 'sites', 3)
    if cutoff is not None:
        integer(cutoff, 'cutoff', 0)
    if type(mu) not in (int, F) or mu <= 3:
        raise ValueError('exact spectral resolution above bounded interaction required')
    mu = F(mu)
    matrix = [[F(0) for _ in range(sites)] for _ in range(sites)]
    for column in range(sites):
        term = {(column, 0):1/(mu+diagonal(q, column, 0, sites))}
        for n in range(degree+1):
            for row in range(sites):
                matrix[row][column] += term.get((row, 0), 0)
            if n < degree:
                term = {label:-value/(mu+diagonal(q, *label, sites))
                        for label, value in apply_v(q, term, sites, cutoff).items()}
    if matrix != list(map(list, zip(*matrix))):
        raise ArithmeticError('complete local resolvent polynomial must be Hermitian')
    error = (3/mu)**(degree+1)/(mu-3)
    return matrix, error


def source_replay(q, cutoff, sites=4):
    """Every finite source matrix column, including all cutoff boundaries."""
    source = build_flux_ring(q, cutoff, L=sites)
    labels = source['labels']
    index = {label:i for i, label in enumerate(labels)}
    n = len(labels)
    for j, label in enumerate(labels):
        column = apply_v(q, {label:F(1)}, sites, cutoff)
        add(column, label, diagonal(q, *label, sites))
        actual = {labels[i]:F(float(source['H'][i, j])) for i in range(n) if source['H'][i, j]}
        if actual != column:
            raise ArithmeticError('literal rotor matrix differs')
        if any(source['Gauss'][j]):
            raise ArithmeticError('complete source Gauss constraint differs')
        x, ell = label
        for edge in range(sites):
            if source['flux'][j, edge] != ell+q*int(edge < x):
                raise ArithmeticError('source electric flux differs')
        if source['initial_state'][j] != int(label == (0, 0)):
            raise ArithmeticError('source input differs')
    if index.get((0, 0)) is None:
        raise ArithmeticError('complete rotor input missing')
    return {'charge6':q, 'cutoff':cutoff, 'sites':sites,
            'complete_columns_checked':n, 'Gauss_rows_checked':sites*n}


def action_witness(q=1, mu=F(8), degree=24):
    import sympy as s
    matrix, error = local_resolvent(q, mu, degree)
    # PHP has diagonal at most3q²/2 and bounded interaction norm at most3.
    floor = 1/(F(mu)+F(3*q*q, 2)+3)
    if error >= floor:
        raise ValueError('local action inversion lacks a positive margin')
    effective = s.Matrix(matrix).inv()
    inverse = [[F(effective[i, j]) for j in range(4)] for i in range(4)]
    action_error = error/(floor*(floor-error))
    return {'charge6':q, 'spectral_resolution_mu':str(mu), 'series_degree':degree,
            'coarse_space':'all four positions with loop flux zero',
            'full_Hilbert_space':'l2(Z) tensor C4 of the distinct U1 rotor completion',
            'local_resolvent_intervals':[[interval_packet(x, error) for x in row] for row in matrix],
            'local_action_intervals':[[interval_packet(x, action_error) for x in row] for row in inverse],
            'resolvent_operator_error_upper':interval_packet(error, F(0))[1],
            'action_operator_error_upper':interval_packet(action_error, F(0))[1],
            'unbounded_electric_operator_retained':True,
            'infinite_flux_limit_constructed_on_this_fixed_graph':True,
            'native60_Higgs_source_replaced_or_identified':False,
            'nonabelian_or_spatial_continuum_identification':False}
