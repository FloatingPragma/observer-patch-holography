"""Exact full-CAR coefficient ranks, including anomalous pairing charges.

A modular rank is a lower bound over the rationals. An explicitly checked
independent commuting family gives the matching upper bound; neither an
eigensolver tolerance nor a selected particle-number sector enters the proof.
"""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt

from ..coarse_energy_ward import ward as w
from ..full_fock_field.field import FieldGeometry


PRIME = 65521
PAULI = (((0, 1), (1, 0)), ((0, -1j), (1j, 0)), ((1, 0), (0, -1)))


def checked_geometry(geometry):
    if type(geometry) is not FieldGeometry or geometry.volume > 8:
        raise ValueError('an unchanged open box with at most eight sites is required')
    if any(F(value).denominator > 2**16
           for value in (geometry.speed, geometry.interaction)):
        raise ValueError('source coefficients must be binary rationals of denominator <= 2**16')
    return geometry


def coefficient(value):
    return F(value.real), F(value.imag)


def source(geometry, *, centered):
    """Reconstruct the source links independently of its emitted hopping groups."""
    checked_geometry(geometry)
    if type(centered) is not bool:
        raise ValueError('a Boolean energy convention is required')
    ladders = [w.annihilator(p) for p in range(geometry.modes)]
    number = w.total(a.adjoint()@a for a in ladders)
    h = w.Clifford({})
    sites = {site: index for index, site in enumerate(geometry.sites)}
    c, u = F(geometry.speed), F(geometry.interaction)
    for x, site in enumerate(geometry.sites):
        for axis in range(3):
            neighbor = list(site)
            neighbor[axis] += 1
            y = sites.get(tuple(neighbor))
            if y is None:
                continue
            for s in range(2):
                for t in range(2):
                    z = coefficient(1j*PAULI[axis][s][t])
                    term = (ladders[2*x+s].adjoint()@ladders[2*y+t]).scale(z).scale(c/2)
                    h = h+term+term.adjoint()
        n0 = ladders[2*x].adjoint()@ladders[2*x]
        n1 = ladders[2*x+1].adjoint()@ladders[2*x+1]
        h = h+(n0@n1).scale(u)
    if h != w.total(w.energies(geometry)):
        raise ArithmeticError('the reconstructed action differs from the actual source')
    return h-number.scale(u/2) if centered else h


def quadratic_basis(modes):
    if type(modes) is not int or not 2 <= modes <= 16:
        raise ValueError('two through sixteen fermionic modes are required')
    return tuple(w.Clifford({(1 << p)|(1 << q): (0, 1)})
                 for p, q in combinations(range(2*modes), 2))


def modular_rank(columns):
    """Exact elimination over F_65521; all rational denominators must invert."""
    if type(columns) not in (list, tuple) or not columns:
        raise ValueError('a nonempty finite list of Clifford columns is required')
    if PRIME < 2 or any(PRIME % divisor == 0 for divisor in range(2, isqrt(PRIME)+1)):
        raise ArithmeticError('the rank modulus must be prime')
    pivots = {}
    for column in columns:
        if type(column) is not w.Clifford:
            raise ValueError('exact Clifford columns are required')
        vector = {}
        for mask, values in column.terms.items():
            for part, value in enumerate(values):
                if value.denominator % PRIME == 0:
                    raise ValueError('coefficient denominator is singular modulo the prime')
                entry = value.numerator*pow(value.denominator, -1, PRIME) % PRIME
                if entry:
                    vector[2*mask+part] = entry
        while vector:
            pivot = min(vector)
            if pivot not in pivots:
                inverse = pow(vector[pivot], -1, PRIME)
                pivots[pivot] = {key: value*inverse % PRIME for key, value in vector.items()}
                break
            factor = vector[pivot]
            for key, value in pivots[pivot].items():
                entry = (vector.get(key, 0)-factor*value) % PRIME
                if entry:
                    vector[key] = entry
                else:
                    vector.pop(key, None)
    return len(pivots)


def rational_rank(columns):
    """Independent rational row reduction for the small boundary controls."""
    import sympy as s
    keys = sorted({(mask, part) for column in columns for mask, values in column.terms.items()
                   for part, value in enumerate(values) if value})
    rows = {key: index for index, key in enumerate(keys)}
    entries = {(rows[mask, part], j): s.Rational(value)
               for j, column in enumerate(columns) for mask, values in column.terms.items()
               for part, value in enumerate(values) if value}
    return int(s.MutableSparseMatrix(len(keys), len(columns), entries).rank())


def generators(geometry, *, centered):
    """Traceless N-V, hidden spin, and (for K) the two eta-pair generators."""
    checked_geometry(geometry)
    if type(centered) is not bool:
        raise ValueError('a Boolean energy convention is required')
    a = [w.annihilator(p) for p in range(geometry.modes)]
    charges = [w.total(b.adjoint()@b for b in a)-w.Clifford.scalar(geometry.volume)]
    for axis in range(3):
        terms = []
        for x, site in enumerate(geometry.sites):
            sign = (-1)**sum(site[j] for j in range(3) if j != axis)
            for s in range(2):
                for t in range(2):
                    terms.append((a[2*x+s].adjoint()@a[2*x+t])
                                 .scale(coefficient(PAULI[axis][s][t])).scale(sign))
        charges.append(w.total(terms))
    if centered:
        pair = w.total((a[2*x].adjoint()@a[2*x+1].adjoint()).scale((-1)**sum(site))
                       for x, site in enumerate(geometry.sites))
        charges.extend((pair+pair.adjoint(), (pair-pair.adjoint()).scale((0, -1))))
    return tuple(charges)


def classification(geometry=FieldGeometry((2, 2, 2)), *, centered=True):
    """A matching exact lower/upper rank certificate on the full quadratic space."""
    checked_geometry(geometry)
    if geometry.speed <= 0 or geometry.interaction <= 0:
        raise ValueError('the classification requires strictly positive hopping and contact')
    h = source(geometry, centered=centered)
    basis = quadratic_basis(geometry.modes)
    charges = generators(geometry, centered=centered)
    for q in charges:
        if q != q.adjoint() or any(mask.bit_count() != 2 for mask in q.terms):
            raise ArithmeticError('a charge is not a traceless self-adjoint quadratic')
        if not w.commutator(h, q).is_zero():
            raise ArithmeticError('an exhibited charge fails complete-source conservation')
    lower = modular_rank(tuple(w.commutator(h, q) for q in basis))
    independent = modular_rank(charges)
    if independent != len(charges) or lower+independent != len(basis):
        raise ArithmeticError('commutator rank and explicit kernel do not exhaust the algebra')
    return {
        'sides': list(geometry.sides), 'fock_dimension': 1 << geometry.modes,
        'clock': 'K=H-U*N/2' if centered else 'H',
        'speed': str(F(geometry.speed)), 'interaction': str(F(geometry.interaction)),
        'quadratic_dimension_excluding_identity': len(basis),
        'prime': PRIME, 'modular_rank_lower_bound': lower,
        'independent_exact_kernel_dimension': independent,
        'rational_and_real_commutator_rank': lower,
        'identity_dimension_additional': 1,
        'full_fock_operator_identities': True,
        'anomalous_pairings_in_search': True,
    }


def spectator_witnesses(sites=8):
    """Every offsite pair is isolated by vacuum and one-spectator transitions."""
    if type(sites) is not int or not 2 <= sites <= 64:
        raise ValueError('two through 64 sites are required')
    def create(mask, p):
        if mask >> p & 1:
            return None, 0
        return mask | (1 << p), (-1)**(mask & ((1 << p)-1)).bit_count()
    def pair(mask, p, q):
        target, sign = create(mask, q)
        if target is None:
            return None, 0
        target, second = create(target, p)
        return target, sign*second
    def doubles(mask):
        return sum(((mask >> (2*x)) & 3) == 3 for x in range(sites))
    count = 0
    for p, q in combinations(range(2*sites), 2):
        if p//2 == q//2:
            continue
        spectator = 1 << (p ^ 1)
        final, sign = pair(spectator, p, q)
        vacuum_final, vacuum_sign = pair(0, p, q)
        matches = [(i, j) for i, j in combinations(range(2*sites), 2)
                   if pair(spectator, i, j)[0] == final]
        if matches != [(p, q)] or abs(sign) != 1 or abs(vacuum_sign) != 1:
            raise ArithmeticError('an offsite pair coefficient is not isolated')
        if doubles(final)-doubles(spectator) != 1 or doubles(vacuum_final) != 0:
            raise ArithmeticError('the pair and spectator have incorrect contact-energy differences')
        count += 1
    return {'sites': sites, 'isolated_complex_offsite_pair_coefficients': count,
            'vacuum_contact_difference': 0, 'spectator_contact_difference': 1,
            'both_pair_creation_signs_retained': True}


def full_fock_check(geometry=FieldGeometry((2, 2, 1)), *, centered=True):
    """Independent occupation-matrix commutator rank on every Fock input.

    Source Gaussian-integer matrices use their own occupation routing. Their
    integer columns are read as rational Clifford containers only for exact
    rank elimination; these container indices label matrix entries here.
    """
    from ..field_stress.energy import GaussianInteger as GI, commutator
    from ..field_symmetry.symmetry import operators, _annihilator
    checked_geometry(geometry)
    if geometry.volume > 4:
        raise ValueError('independent full-matrix ranks admit at most four sites')
    data = operators(geometry)
    h = data.K if centered else data.H
    annihilators = [_annihilator(geometry.modes, p) for p in range(geometry.modes)]
    gammas = []
    for a in annihilators:
        gammas.extend((a+a.adjoint(), (a.adjoint()-a).times_i()))
    columns = []
    dimension = 1 << geometry.modes
    for p, q in combinations(range(2*geometry.modes), 2):
        charge = (gammas[p]@gammas[q]).times_i()
        residual = commutator(h, charge)
        real, imag = residual.real.tocoo(), residual.imag.tocoo()
        entries = {}
        for part, matrix in enumerate((real, imag)):
            for row, col, value in zip(matrix.row, matrix.col, matrix.data, strict=True):
                key = int(row)*dimension+int(col)
                pair = list(entries.get(key, (0, 0)))
                pair[part] = int(value)
                entries[key] = tuple(pair)
        columns.append(w.Clifford(entries))
    rank = modular_rank(columns)
    expected = classification(geometry, centered=centered)
    if rank != expected['rational_and_real_commutator_rank']:
        raise ArithmeticError('independent full-Fock and symbolic ranks differ')
    return {'sides': list(geometry.sides), 'clock': expected['clock'],
            'all_input_output_matrix_entries': dimension**2,
            'quadratic_columns': len(columns), 'exact_rank': rank}


def boundaries():
    geometry = FieldGeometry((2, 1, 1), interaction=0.)
    h = source(geometry, centered=True)
    free_rank = rational_rank(tuple(w.commutator(h, q) for q in quadratic_basis(4)))
    if 28-free_rank <= 6:
        raise ArithmeticError('the free control failed to enlarge the conserved algebra')
    stopped = FieldGeometry((2, 1, 1), speed=0.)
    stopped_rank = rational_rank(tuple(w.commutator(source(stopped, centered=True), q)
                                       for q in quadratic_basis(4)))
    if 28-stopped_rank != 12:
        raise ArithmeticError('disconnected onsite centered charges were lost')
    local_pair = (w.annihilator(0).adjoint()@w.annihilator(1).adjoint())
    local_pair = local_pair+local_pair.adjoint()
    actual = source(FieldGeometry((2, 1, 1)), centered=True)
    if w.commutator(actual, local_pair).is_zero():
        raise ArithmeticError('a strictly localized eta pair was incorrectly conserved')
    number = generators(FieldGeometry((2, 1, 1)), centered=True)[0]
    detuned = actual-number.scale(F(1, 16))
    detuned_rank = rational_rank(tuple(w.commutator(detuned, q) for q in quadratic_basis(4)))
    if 28-detuned_rank != 4:
        raise ArithmeticError('a detuned chemical shift incorrectly retains eta pairing')
    # A separate real-hopping triangle has no alternating eta profile.
    ladders = [w.annihilator(p) for p in range(6)]
    triangle = w.Clifford({})
    for x in range(3):
        n0 = ladders[2*x].adjoint()@ladders[2*x]
        n1 = ladders[2*x+1].adjoint()@ladders[2*x+1]
        triangle = triangle+n0@n1-(n0+n1).scale(F(1, 2))
        y = (x+1) % 3
        for spin in range(2):
            hop = ladders[2*x+spin].adjoint()@ladders[2*y+spin]
            triangle = triangle+(hop+hop.adjoint()).scale(F(1, 16))
    triangle_rank = rational_rank(tuple(w.commutator(triangle, q) for q in quadratic_basis(6)))
    if 66-triangle_rank != 4:
        raise ArithmeticError('the odd-cycle control incorrectly retains eta pairing')
    return {'U_zero_centered_dimer_quadratic_kernel_dimension': 28-free_rank,
            'c_zero_centered_two_site_kernel_dimension': 28-stopped_rank,
            'detuned_centered_dimer_quadratic_kernel_dimension': 28-detuned_rank,
            'separate_real_hopping_triangle_kernel_dimension': 66-triangle_rank,
            'proper_site_eta_pair_conserved': False}


def certificate():
    return {
        'classifications': [classification(FieldGeometry(sides), centered=centered)
                            for sides in ((2, 1, 1), (1, 2, 1), (1, 1, 2),
                                          (3, 1, 1), (2, 2, 1), (2, 2, 2))
                            for centered in (False, True)],
        'spectator_witnesses': spectator_witnesses(),
        'independent_full_fock': [full_fock_check(centered=centered)
                                  for centered in (False, True)],
        'boundary_controls': boundaries(),
        'physical_species_identified': False,
        'continuum_stress_or_interacting_limit_proved': False,
        'higher_degree_or_nonlocal_charges_excluded': False,
    }
