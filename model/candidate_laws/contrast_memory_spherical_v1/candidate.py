"""The interface of A-contrast-memory-spherical-v1 to the executable model and to the evaluation.

The law runs on the fixed carrier of ``model.executable_model``. Its one-load generator is assembled here from the
local seam term of ``law.py`` and the degree term of SPEC.md section 4, so the evaluation checks the implemented
law itself.
"""
from __future__ import annotations

import numpy as np

from model.executable_model import CELL_SEAMS

from . import law

NAME = 'contrast_memory_spherical_v1'
TITLE = 'A-contrast-memory-spherical-v1'
ROLE = 'reference'
SUMMARY = ('Loads 0, 1 or 2 on the ports and an eleven-state memory on every seam. A seam event writes the '
           'contrast of its two loads into the memory, gives it a phase and erases it again; on loads this is '
           'an exchange plus an agreement energy kappa.')
DOCUMENTS = {'readme': 'README.md', 'specification': 'SPEC.md', 'proof': 'PROOF.md'}
LABELLING = CELL_SEAMS
FEDERATION_LEVELS = ('cell', 0, 1, 2)
KAPPA = 8.0

PARAMETERS = {
    'kappa': {'value': KAPPA, 'provenance': 'designed', 'quantity': 'agreement energy in units of the exchange',
              'note': 'tested at 0 to 128; the headline values are 8 and 128'},
    'emitter_strength': {'value': 64.0, 'provenance': 'designed', 'quantity': 'coupling of the local emitter',
                         'note': 'enters the error budget of the emitted pair'},
    'reader_apparatus': {'value': 'SPEC.md section 7', 'provenance': 'designed',
                         'quantity': 'reader registers and their complete outcome set'},
}


def run(argv=None):
    """The replay of this candidate: its headline results against RESULTS.json (``--full`` recomputes all)."""
    from .verify import main
    main(argv)
    return 0


def basis(car):
    return [('port', p) for p in range(car.n)]


def generator(car, kappa: float = KAPPA):
    """The law on one load, ``sum_e h_e + kappa sum_v (6 - d_v) n_v^2``, from the implemented seam term."""
    from scipy import sparse
    h = law.seam_term(kappa)
    one = {(1, 0): 3, (0, 1): 1}
    rows, cols, vals = [], [], []
    diag = np.zeros(car.n)
    for a, b in car.seams:
        diag[a] += h[one[(1, 0)], one[(1, 0)]]
        diag[b] += h[one[(0, 1)], one[(0, 1)]]
        rows += [b, a]
        cols += [a, b]
        vals += [h[one[(0, 1)], one[(1, 0)]], h[one[(1, 0)], one[(0, 1)]]]
    diag += kappa * (6 - car.degree)
    H = sparse.coo_matrix((vals, (rows, cols)), shape=(car.n, car.n)) + sparse.diags(diag)
    return H.tocsr()


def response(car):
    """The operator whose slowest nonconstant level is the carrier response: the one-load law."""
    return generator(car).toarray()


def represent(rotation, car):
    """A port permutation acting on the one-load sector."""
    P = np.zeros((car.n, car.n))
    for p, q in enumerate(rotation):
        P[q, p] = 1.0
    return P


# Evidence for the premises of the premise register (tracking/premise_register.json). Unlisted premises are open.
PREMISES = {
    'PR-43': {'status': 'fulfilled',
              'evidence': ['SPEC.md', 'law.py'],
              'note': 'The law is the conjugation flow of one Hamiltonian on the federation, a pointwise-continuous '
                      'group of star-automorphisms of the public algebra.'},
    'PR-44': {'status': 'partial',
              'evidence': ['SPEC.md', 'events.py'],
              'note': 'Ports, seams and seam memories give a finite subsystem split, and seam events are normalized '
                      'local operations. The schedule of the events is supplied.'},
    'PR-02': {'status': 'partial',
              'evidence': ['reader.py', 'presence.py', 'PROOF.md'],
              'note': 'The reader and the presence readings are complete finite instruments on the state of the law; '
                      'the algebra-state representation of the records is supplied with them.'},
}

# The joint requirements of the scorecard: premises fulfilled together, on the same objects.
JOINT = {
    'relativistic-matter': {
        'status': 'fails',
        'evidence': ['../../experiments/material_inertia/README.md'],
        'note': 'At kappa = 8 the kinetic ratios of clusters of two, three and four loads (about 34, 1,155 and '
                '39,236) differ from their rest-energy ratios (about 1.686, 2.354 and 3.020); a common relativistic '
                'mass would make them equal.',
        'counterexample': 'experiments/material_inertia, binding challenge and flight records'},
}
