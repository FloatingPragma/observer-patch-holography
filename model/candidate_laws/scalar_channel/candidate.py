"""The interface of the scalar channel to the executable model and to the evaluation.

The scalar channel is written on one cell in its own labelling of the ports (``carrier.EDGES``); the evaluation
maps the rotations of the fixed cell into that labelling. Its one-object law is the matrix ``D`` on the 30 seam
modes and the 120 wedge modes of one cell.
"""
from __future__ import annotations

import numpy as np

from . import carrier as scalar

NAME = 'scalar_channel'
TITLE = 'A-scalar-CAR-pair-transit-v2 (the scalar channel)'
ROLE = 'second'
SUMMARY = ('One object or none on every seam and on every wedge, a wedge being two seams at one port. An object '
           'hops between a seam and a wedge, which sit at the two opposite levels +m and -m.')
DOCUMENTS = {'readme': 'README.md', 'proof': 'PROOF.md'}
LABELLING = scalar.EDGES
FEDERATION_LEVELS = ('cell',)

PARAMETERS = {
    'm': {'value': 1.0, 'provenance': 'designed', 'quantity': 'rest level of the seam and wedge modes'},
    'half_link_coefficient': {'value': '1/sqrt(2)', 'provenance': 'designed',
                              'quantity': 'amplitude of the seam-wedge hop'},
    'services': {'value': 'protocol.json', 'provenance': 'designed',
                 'quantity': 'loop instruments, reference ramps, pulses and pointers'},
}


def run(argv=None):
    """The replay of this candidate: preparation, pulses and pointers against results.json."""
    from .experiment import main
    return main(argv)


def _source():
    return scalar.load_source('cell')


def basis(car=None):
    data = _source()
    return [('seam', tuple(map(int, e))) for e in data['edges']] + [('wedge', tuple(map(int, w))) for w in data['wedges']]


def generator(car=None):
    if car is not None and car.cells > 1:
        raise NotImplementedError('the scalar channel is written on one cell')
    return np.array(_source()['D'], dtype=float)


def response(car=None):
    """The levels of the seams: omega^2 = m^2 + lambda with lambda the eigenvalues of B B^T."""
    B = _source()['B']
    return B @ B.T


def represent(rotation, car=None):
    """A port permutation acting on seams and wedges; a wedge read in the opposite direction changes sign."""
    data = _source()
    labels = basis()
    index = {lab: k for k, lab in enumerate(labels)}
    n = len(labels)
    P = np.zeros((n, n))
    for k, (kind, where) in enumerate(labels):
        if kind == 'seam':
            image = ('seam', tuple(sorted(rotation[v] for v in where)))
            P[index[image], k] = 1.0
        else:
            a, b, c = (rotation[v] for v in where)
            if ('wedge', (a, b, c)) in index:
                P[index[('wedge', (a, b, c))], k] = 1.0
            else:
                P[index[('wedge', (c, b, a))], k] = -1.0
    del data
    return P


PREMISES = {
    'PR-43': {'status': 'fulfilled',
              'evidence': ['PROOF.md', 'carrier.py'],
              'note': 'The law is the conjugation flow of the quadratic generator D, a pointwise-continuous group of '
                      'star-automorphisms of the public algebra.'},
    'PR-02': {'status': 'partial',
              'evidence': ['PROOF.md', 'instruments.py'],
              'note': 'Loop instruments and the complete law of two pointers are finite instruments with their '
                      'back-action; the representation of the records is supplied with them.'},
}

JOINT = {}
