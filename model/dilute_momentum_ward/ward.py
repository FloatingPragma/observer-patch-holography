"""Exact periodic local Ward columns and rigorous dilute-band error bounds.

The band projector describes initial states, not a modified Hamiltonian.
All subsequent evolution and contact forces use the complete source action.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import isqrt

import sympy as s

from ..coarse_energy_ward import ward as w
from ..lattice_momentum_balance.momentum import lift, site_restriction, contact_edge_force
from ..periodic_momentum import momentum as p


def columns(geometry=p.PeriodicGeometry()):
    """All local sine-momentum/stress/force identities on the periodic source."""
    h = p.one_particle(geometry)
    H, H0, V = p.action(geometry)
    index = {site: x for x, site in enumerate(geometry.sites)}
    oriented = []
    for x, site in enumerate(geometry.sites):
        for axis, length in enumerate(geometry.sides):
            if length == 1:
                continue
            target = list(site); target[axis] = (target[axis]+1) % length
            oriented.append((x, index[tuple(target)], axis))
    edges = sorted({tuple(sorted((x, y))) for x, y, _ in oriented})
    result = []
    for axis in range(3):
        permutation = geometry.shift(axis)
        t = s.ImmutableSparseMatrix(geometry.modes, geometry.modes,
                                    {(target, source): 1 for source, target in enumerate(permutation)})
        sine = (t-t.T)/(2*s.I)
        if h*sine != sine*h:
            raise ArithmeticError('periodic free sine momentum is not conserved')
        local = tuple(lift(site_restriction(sine, x)) for x in range(len(geometry.sites)))
        stress = {}
        for x, y in edges:
            a = s.MutableSparseMatrix(geometry.modes, geometry.modes, {})
            for i, j in product(range(2), repeat=2):
                a[2*x+i, 2*y+j] = h[2*x+i, 2*y+j]
                a[2*y+j, 2*x+i] = -h[2*y+j, 2*x+i]
            value = lift(s.I*(a*sine+sine*a)/2)
            if not value.is_zero():
                stress[x, y] = value
        force = tuple(w.commutator(V, value).scale((0, 1)) for value in local)
        explicit = [w.Clifford({}) for _ in local]
        for x, y, edge_axis in oriented:
            if edge_axis == axis:
                value = contact_edge_force(x, y, geometry.interaction).scale(F(-1, 2))
                explicit[x] = explicit[x]+value
                explicit[y] = explicit[y]+value
        if tuple(explicit) != force:
            raise ArithmeticError('periodic sine-force orientation changed')
        divergence = w.divergence(stress, len(local))
        for x, density in enumerate(local):
            if not (w.commutator(H, density).scale((0, 1))+divergence[x]-force[x]).is_zero():
                raise ArithmeticError('exact periodic local Ward identity failed')
        for widths in product(*(range(1, n+1) for n in geometry.sides)):
            labels = w.block_labels(geometry.sites, widths)
            if w.divergence(w.aggregate_flows(stress, labels), max(labels)+1) != w.aggregate(divergence, labels):
                raise ArithmeticError('periodic Ward divergence failed block aggregation')
        result.append({'local': local, 'stress': stress, 'force': force, 'sine': sine})
    return {'H': H, 'H0': H0, 'V': V, 'columns': tuple(result), 'oriented_edges': tuple(oriented)}


def band_counts(side=64, radius=1, dimension=3):
    """Exact contact principal angles for a cubic spin-complete momentum set."""
    if (type(side) is not int or not 3 <= side <= 4096 or
            type(radius) is not int or not 0 <= radius < side//2 or
            type(dimension) is not int or not 1 <= dimension <= 3 or
            (2*radius+1)**dimension > 4096):
        raise ValueError('bounded cubic momentum grid and proper centered band required')
    band = tuple(product(range(-radius, radius+1), repeat=dimension))
    convolution = Counter(tuple((x+y) % side for x, y in zip(left, right, strict=True))
                          for left in band for right in band)
    volume = side**dimension
    fractions = {key: F(count, volume) for key, count in sorted(convolution.items())}
    fraction = F(len(band), volume)
    if max(fractions.values()) != fraction or sum(convolution.values()) != len(band)**2:
        raise ArithmeticError('complete contact overlap convolution failed')
    return {'side': side, 'dimension': dimension, 'band': band, 'volume': volume,
            'fraction': fraction, 'contact_eigenvalues': fractions,
            'contact_norm_squared': max(fractions.values()),
            'initial_escape_norm_squared': max(value*(1-value) for value in fractions.values())}


def sqrt_upper(value, bits=48):
    if type(value) not in (int, F) or value < 0 or type(bits) is not int or not 1 <= bits <= 256:
        raise ValueError('nonnegative exact rational and bounded precision required')
    value = F(value)
    scale = 1 << bits
    numerator = value.numerator*scale*scale
    root = isqrt(numerator//value.denominator)
    if root*root*value.denominator < numerator:
        root += 1
    return F(root, scale)


def budget(fraction, maximum_number, duration, interaction=F(1), *, cutoff=None, block_sites=1, start_time=F()):
    """Exact outward local-cell and optional two-particle momentum bounds."""
    if (type(fraction) not in (int, F) or not 0 < fraction <= 1 or
            type(maximum_number) is not int or not 1 <= maximum_number <= 4096 or
            any(type(v) not in (int, F) or v < 0 for v in (duration, interaction, start_time)) or
            type(block_sites) is not int or block_sites < 1):
        raise ValueError('exact nonnegative data, band fraction in (0,1], and positive counts required')
    f, t, u, start = F(fraction), F(duration), F(interaction), F(start_time)
    d = maximum_number//2
    root = sqrt_upper(f)
    collision = u*d*root
    leakage = min(F(1), (start+t)*collision)
    local = min(u*t, 4*u*f*t+u*collision*(start*t+t*t/2)) if d else F()
    result = {'sqrt_fraction_upper': root, 'interaction_on_initial_band_norm_upper': collision,
              'evolved_band_leakage_upper': leakage, 'local_Ward_cell_defect_upper': local,
              'coarse_block_defect_upper': block_sites*local,
              'previous_leakage_upper': min(F(1), (start+t)*u*d),
              'interval_start_after_preparation': start, 'interval_end_after_preparation': start+t}
    if cutoff is not None:
        if type(cutoff) not in (int, F) or cutoff < 0 or 2*cutoff >= 3 or maximum_number != 2:
            raise ValueError('two-particle cutoff requires exact kappa>=0 and 2*kappa<3')
        kappa = F(cutoff)
        result['two_particle_local_momentum_error_upper'] = min(4*kappa, kappa**3/3+4*kappa*leakage)
        result['rescaled_by_cutoff_momentum_error_upper'] = (
            result['two_particle_local_momentum_error_upper']/kappa if kappa else F())
    return result


def occupation_matrix(operator, modes, number=2):
    """Exact fixed-number matrix for diagnostics; the Ward proof is full CAR."""
    if type(modes) is not int or not 2 <= modes <= 20 or number != 2:
        raise ValueError('two-particle diagnostic on two through20 modes required')
    masks = tuple(sum(1 << j for j in pair) for pair in combinations(range(modes), 2))
    indices = {mask: i for i, mask in enumerate(masks)}
    result = {}
    for column, initial in enumerate(masks):
        outputs = {}
        for mask, (a, b) in operator.terms.items():
            if mask.bit_length() > 2*modes:
                raise ValueError('operator exceeds the diagnostic carrier')
            current, coefficient = initial, s.Rational(a)+s.I*s.Rational(b)
            for gamma in reversed(range(mask.bit_length())):
                if not mask >> gamma & 1:
                    continue
                mode, kind = divmod(gamma, 2)
                coefficient *= (-1)**((current & ((1 << mode)-1)).bit_count())
                if kind:
                    coefficient *= -s.I if current >> mode & 1 else s.I
                current ^= 1 << mode
            outputs[current] = outputs.get(current, 0)+coefficient
        for mask, coefficient in outputs.items():
            coefficient = s.expand(coefficient)
            if coefficient:
                if mask not in indices:
                    raise ArithmeticError('diagnostic operator does not preserve number two')
                result[indices[mask], column] = coefficient
    return s.ImmutableSparseMatrix(len(masks), len(masks), result), masks


def zero_pair(masks, sites):
    values = {}
    for j, mask in enumerate(masks):
        occupied = [i for i in range(2*sites) if mask >> i & 1]
        if occupied[0] % 2 != occupied[1] % 2:
            values[j, 0] = s.Rational(1 if occupied[0] % 2 == 0 else -1, sites)
    return s.ImmutableSparseMatrix(len(masks), 1, values)
