"""Measure conserved, additive regional energy on the complete sphere in the distinct-port two-load sector."""
import argparse
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model import CELL_SEAMS, federation
from .experiment import (ANGLES, HERE, KAPPA, STEPS, WEIGHTS, Bond, CellEnergyMeter,
                         check_native_gate, portable, states)


class InternalMeter:
    """Local word probabilities read the reference energy at every degree-five port."""

    def __init__(self, sector, owner):
        offset = 12 * owner
        bonds = [Bond(sector, (offset+a, offset+b)) for a, b in CELL_SEAMS]
        self.seams = CellEnergyMeter(bonds)
        reference = KAPPA * (6-sector.car.degree) * (np.arange(sector.n)//12 == owner)
        self.reference = reference[sector.aa] + reference[sector.bb]
        self.operator = sum(b.h for b in bonds) + sparse.diags(self.reference)

    def read(self, wave):
        return self.seams.read(wave) + float(np.dot(self.reference, abs(wave)**2))


class Region:
    def __init__(self, owners, meters, board, bridges):
        self.owners = owners
        self.meters = [meters[c] for c in owners]
        self.boundary = {}
        self.storage = []
        self.operator = sum(m.operator for m in self.meters)
        for edge in board.external:
            left, right = (port//12 in owners for port in edge)
            if left or right:
                share = 1. if left and right else .5
                self.storage.append((share, bridges[edge]))
                self.operator = self.operator + share * bridges[edge].h
            if left != right:
                self.boundary[edge] = 1 if left else -1

    def energies(self, state):
        internal = sum(m.read(state) for m in self.meters)
        total = internal
        for share, bond in self.storage:
            left, right = state[bond.left], state[bond.right]
            total += share * (KAPPA * (np.vdot(left, left).real + np.vdot(right, right).real)
                              + .5 * np.vdot(left-right, left-right).real)
        return np.array([total, internal])


def commutator(a, b):
    return a @ b - b @ a


def run():
    board = federation(0)
    sector = PairSector(board)
    meters = [InternalMeter(sector, c) for c in range(board.cells)]
    bridges = {edge: Bond(sector, edge) for edge in board.external}
    if len({p for edge in bridges for p in edge}) != 2 * len(bridges):
        raise ArithmeticError("the external contacts must form a matching")
    full = sum(m.operator for m in meters) + sum(b.h for b in bridges.values())
    flight = sector.hamiltonian(KAPPA)
    missing = full - flight - (12*KAPPA)*sparse.eye(sector.P)
    if missing.count_nonzero():
        raise ArithmeticError("the cell meters and bridge energies do not sum to the native law")
    adjacent = min(board.external)[1]//12
    regions = [Region(s, meters, board, bridges) for s in [(0,), (adjacent,), (0, adjacent)]]
    addresses = sorted(set().union(*(r.boundary for r in regions)))
    currents = {}
    for edge in addresses:
        a, b = (port//12 for port in edge)
        currents[edge] = .5 * commutator(bridges[edge].h, meters[a].operator-meters[b].operator)
    for region in regions:
        total = sum(sign * currents[e] for e, sign in region.boundary.items())
        if (commutator(flight, region.operator)-total).count_nonzero():
            raise ArithmeticError("the region and its oriented boundary differ on a binary column")
    if (regions[0].operator+regions[1].operator-regions[2].operator).count_nonzero():
        raise ArithmeticError("the fused region has a different energy operator")

    rows = []
    for port, sign, launch in states(sector):
        if (port, sign) not in ((0, 1), (9, -1)):
            continue
        for time in (.375, 1.125):
            wave = expm_multiply((-1j*time) * flight, launch)
            responses = {}
            for edge in addresses:
                a, b = (p//12 for p in edge)
                samples = []
                for angle in ANGLES:
                    exposed = bridges[edge].expose(wave, angle)
                    samples.append([meters[a].read(exposed), meters[b].read(exposed)])
                left, right = WEIGHTS @ samples
                responses[edge] = float((left-right)/2)
            slopes = []
            for step in STEPS:
                before = expm_multiply((1j*step)*flight, wave)
                after = expm_multiply((-1j*step)*flight, wave)
                slopes.append(np.array([(r.energies(after)-r.energies(before))/(2*step)
                                        for r in regions]))
            regional = []
            for i, region in enumerate(regions):
                incoming = [orientation*responses[e] for e, orientation in region.boundary.items()]
                prediction = sum(incoming)
                free, internal = slopes[-1][i]
                regional.append(dict(cells=list(region.owners), inflow=prediction, free_slope=free,
                                     residual=prediction-free, step_change=abs(free-slopes[0][i,0]),
                                     omit_first_residual=prediction-incoming[0]-free,
                                     omit_storage_residual=prediction-internal))
            rows.append(dict(port=port, sign=sign, exposure=time,
                             currents=[dict(ports=list(edge), inflow_to_first=responses[edge])
                                       for edge in addresses], regions=regional))
    all_regions = [r for row in rows for r in row["regions"]]
    summary = {metric: max(abs(r[metric]) for r in all_regions) for metric in
               ("residual", "step_change", "omit_first_residual", "omit_storage_residual")}
    if summary["residual"] >= 1e-6 or min(summary[k] for k in
                                           ("omit_first_residual", "omit_storage_residual")) <= 1e-3:
        raise ArithmeticError("the comparison or its sensitivity control failed")
    return portable(dict(cells=board.cells, active_seams=len(board.seams),
                         degree_reference_ports=int(np.sum(board.degree==5)), states=sector.P,
                         checked_regions=[dict(cells=list(r.owners), boundary=len(r.boundary),
                                               nonzero_operator_residual_entries=0) for r in regions],
                         native_column_error=check_native_gate(), readings=rows, maximum=summary))


def verify(reference, actual):
    if isinstance(reference, dict):
        if reference.keys() != actual.keys():
            raise ValueError("receipt fields changed")
        for key, value in reference.items():
            verify(value, actual[key])
    elif isinstance(reference, list):
        for left, right in zip(reference, actual, strict=True):
            verify(left, right)
    elif isinstance(reference, float):
        np.testing.assert_allclose(actual, reference, rtol=0, atol=8e-9)
    elif actual != reference:
        raise ValueError("receipt changed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="replace the regional measurement receipt")
    args = parser.parse_args()
    result = run()
    path = HERE/"regional.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        verify(json.loads(path.read_text()), result)
    print(json.dumps(result["maximum"], indent=2))


if __name__ == "__main__":
    main()
