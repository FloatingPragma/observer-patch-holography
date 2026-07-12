"""Real lattice QCD engine, diagnostic scale, for the OPH hadron lane.

Quenched SU(3) gauge generation, clover-Wilson valence propagators, pion and
nucleon correlators, and a two-flavor pseudofermion HMC branch. Every number
this package emits is computed from the lattice dynamics; nothing is anchored
to target masses. Its output class is diagnostic and non-promoting: the #425
production closure additionally requires HPC-scale ensembles and the full
systematics program, which no local run can supply.
"""

from .core import average_plaquette, cold_start, sweep, wilson_gauge_action
from .dirac import WilsonClover, cg_normal, point_propagator
from .hmc import TwoFlavorHMC
from .spectroscopy import effective_mass, nucleon_correlators, pion_correlator

__all__ = [
    "average_plaquette",
    "cold_start",
    "sweep",
    "wilson_gauge_action",
    "WilsonClover",
    "cg_normal",
    "point_propagator",
    "TwoFlavorHMC",
    "effective_mass",
    "nucleon_correlators",
    "pion_correlator",
]
