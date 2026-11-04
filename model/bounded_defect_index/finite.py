"""Independent finite original-source controls; no infinite proof via a fit."""
from itertools import product

import numpy as np
from scipy.sparse import coo_matrix

from ..native_mass_channels.algebra import SIGMA, rational
from .algebra import orientation
from .bounds import parameters


def sites(radius):
    if type(radius) is not int or radius < 0:
        raise ValueError('nonnegative integer source cube radius required')
    return tuple(product(range(-radius, radius+1), repeat=3))


def source_matrix(radius, spacing=1, velocity=1, mass=1, sign=1, hopping_fraction=1):
    a, v, m = map(float, parameters(spacing, velocity, mass))
    sign = orientation(sign)
    fraction = rational(hopping_fraction)
    if not 0 <= fraction <= 1:
        raise ValueError('the source hopping homotopy parameter must lie in [0,1]')
    points = sites(radius)
    addresses = {point: j for j, point in enumerate(points)}
    sigma = [np.array(x, complex) for x in SIGMA]
    rows, columns, values = [], [], []
    def add_block(row, column, block):
        for spin in range(2):
            for other in range(2):
                if block[spin, other] != 0:
                    rows.append(2*row+spin)
                    columns.append(2*column+other)
                    values.append(block[spin, other])
    for point, address in addresses.items():
        norm = np.linalg.norm(point)
        if norm:
            block = sign*m*sum(point[j]*(-1)**(point[j] % 2)*sigma[j] for j in range(3))/norm
            add_block(address, address, block)
        for axis in range(3):
            other = tuple(n+int(j == axis) for j, n in enumerate(point))
            if other in addresses:
                block = float(fraction)*1j*v/(2*a)*sigma[axis]
                add_block(address, addresses[other], block)
                add_block(addresses[other], address, block.conj().T)
    return coo_matrix((values, (rows, columns)), shape=(2*len(points), 2*len(points))).tocsr()


def inversion_bases(radius):
    points = sites(radius)
    addresses = {point: j for j, point in enumerate(points)}
    pairs = [point for point in points if point < tuple(-n for n in point)]
    plus_rows, plus_columns, plus_values = [], [], []
    minus_rows, minus_columns, minus_values = [], [], []
    origin = addresses[(0, 0, 0)]
    for spin in range(2):
        plus_rows.append(2*origin+spin)
        plus_columns.append(spin)
        plus_values.append(1.)
    for number, point in enumerate(pairs):
        opposite = tuple(-n for n in point)
        for spin in range(2):
            for target, sign in ((point, 1), (opposite, -1)):
                plus_rows.append(2*addresses[target]+spin)
                plus_columns.append(2+2*number+spin)
                plus_values.append(1/np.sqrt(2))
                minus_rows.append(2*addresses[target]+spin)
                minus_columns.append(2*number+spin)
                minus_values.append(sign/np.sqrt(2))
    dimension = 2*len(points)
    plus = coo_matrix((plus_values, (plus_rows, plus_columns)), shape=(dimension, len(points)+1)).tocsr()
    minus = coo_matrix((minus_values, (minus_rows, minus_columns)), shape=(dimension, len(points)-1)).tocsr()
    return plus, minus


def compact_odd_density(radius, coupling=1):
    strength = float(rational(coupling))
    points = sites(radius)
    diagonal = [strength*point[0] if max(map(abs, point)) <= 1 else 0.
                for point in points for _ in range(2)]
    return coo_matrix((diagonal, (range(len(diagonal)), range(len(diagonal))))).tocsr()
