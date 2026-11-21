"""The executable model: everything that is fixed for every candidate law.

The carrier is the spherical echosahedral federation. One cell, the echosahedron, has 12 ports, 30 seams and
20 faces with the incidence and the sixty rotations of the icosahedron. The federations of 20, 80 and 320 cells
glue complete cells by external seams. A candidate law lives in ``model/candidate_laws/<name>`` and declares its
local generator on this carrier, its parameters with their provenance, and its evidence for the premises of the
finite theorems. ``python -m model`` runs a candidate and ``python -m model evaluate`` scores every candidate.
"""
from . import carrier
from .carrier import CELL_SEAMS, EMITTER_SEAM, NAMES, Carrier, cell, federation
from .symmetry import incident, rotations

__all__ = ['CELL_SEAMS', 'EMITTER_SEAM', 'NAMES', 'Carrier', 'carrier', 'cell', 'federation', 'incident', 'rotations']
