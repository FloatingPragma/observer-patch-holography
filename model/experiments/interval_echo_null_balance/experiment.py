"""Compare 19 physical gateway exposures with independent free-flight energy readings."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import currents
from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import cell, federation
from model.executable_model.carrier import Carrier

HERE = Path(__file__).resolve().parent
ANGLES = 2 * np.pi * np.arange(19) / 19
FREQUENCIES = np.arange(1, 10)
WEIGHTS = 2 * np.sin(ANGLES[:, None] * FREQUENCIES) @ FREQUENCIES / len(ANGLES)
TIMES = (1., 4., 16.)
STEPS = (1e-4, 5e-5)


class Bond:
    """One actual seam and the two one-end-occupied sectors of the pair basis."""

    def __init__(self, sector, edge):
        a, b = edge
        spectators = [j for j in range(sector.n) if j not in edge]
        self.left = np.array([sector.index(a, j) for j in spectators])
        self.right = np.array([sector.index(b, j) for j in spectators])
        i, j = self.left, self.right
        entries = np.repeat([KAPPA + .5, KAPPA + .5, -.5, -.5], len(i))
        self.h = sparse.csr_matrix(
            (entries, (np.r_[i, j, i, j], np.r_[i, j, j, i])),
            shape=(sector.P, sector.P),
        )

    def expose(self, state, angle):
        x, y = state[self.left], state[self.right]
        even = .5 * np.exp(-1j * KAPPA * angle) * (x + y)
        odd = .5 * np.exp(-1j * (KAPPA + 1) * angle) * (x - y)
        out = state.copy()
        out[self.left], out[self.right] = even + odd, even - odd
        return out


class CellEnergyMeter:
    """Sum of the 30 diagonal and antisymmetric-port readout expectations."""

    def __init__(self, bonds):
        self.left = np.stack([bond.left for bond in bonds])
        self.right = np.stack([bond.right for bond in bonds])

    def read(self, state):
        x, y = state[self.left], state[self.right]
        # Quiet and equally occupied endpoints give zero in both settings.
        mixed = (abs(x)**2 + abs(y)**2).sum()
        antisymmetric = .5 * (abs(x - y)**2).sum()
        return float(KAPPA * mixed + antisymmetric)


def hardware():
    board = federation(0)
    links = [edge for edge in board.external if edge[0] < 12][:2]
    nodes = [0, *(b // 12 for _, b in links)]
    port_map = {12 * node + p: 12 * i + p for i, node in enumerate(nodes) for p in range(12)}
    cells = [[(12 * i + a, 12 * i + b) for a, b in cell().seams] for i in range(3)]
    bridges = [(port_map[a], port_map[b]) for a, b in links]
    active = Carrier("three serviced cells", 36, sum(cells, []) + bridges)
    sector = PairSector(active)
    meter_bonds = [[Bond(sector, edge) for edge in edges] for edges in cells]
    link_bonds = [Bond(sector, edge) for edge in bridges]
    cell_h = [sum(bond.h for bond in stack) for stack in meter_bonds]
    source_operator = sum(cell_h) + sum(bond.h for bond in link_bonds)
    # The deliberately closed degree-reference services are absent from both readers.
    expected = sector.exchange() + sparse.diags(
        KAPPA * (active.degree[sector.aa] + active.degree[sector.bb] - 2 * sector.adjacent)
    )
    delta = source_operator - expected
    if delta.nnz and max(abs(delta.data)) > 1e-12:
        raise ArithmeticError("gate decomposition differs from the native pair flight")
    return sector, nodes, links, CellEnergyMeter(meter_bonds[0]), cell_h, link_bonds


def check_native_gate():
    exchange, contrast = native_law.local_term_parts()
    term = exchange + KAPPA * contrast
    ids = [native_law.SRC.index(p) for p in ((0, 0), (1, 0), (0, 1), (1, 1))]
    columns = np.eye(len(term))[:, ids]
    target = np.zeros((4, 4))
    target[1:3, 1:3] = [[KAPPA + .5, -.5], [-.5, KAPPA + .5]]
    err = np.linalg.norm(term @ columns - columns @ target)
    if err > 1e-11 or np.max(abs(np.linalg.eigvalsh(target) - [0, 0, 8, 9])) > 1e-11:
        raise ArithmeticError("19-angle reconstruction requires the native binary spectrum")
    return float(err)


def expectation(state, operator):
    return float(np.vdot(state, operator @ state).real)


def states(sector):
    local = PairSector(cell())
    local_indices = np.array([sector.index(p, q) for p, q in zip(local.aa, local.bb, strict=True)])
    # Reuse the registered inversion-odd source, including its out-of-band weight.
    _, source, _, _, _, _, _, _, _, _ = currents._service(cell(), True, 1.)
    for port in (0, 9):
        for sign in (1, -1):
            start = np.zeros(sector.P, dtype=complex)
            start[local_indices] = currents.pulse(port, sign * currents.STRENGTH) @ source[:, currents.SOURCE_PORT]
            yield port, sign, start


def measure():
    sector, cells, links, meter, cell_h, bonds = hardware()
    closed = sum(cell_h)
    full = closed + sum(b.h for b in bonds)
    energy, vectors = eigh(full.toarray())
    closed_energy, closed_vectors = eigh(closed.toarray())
    residuals, one_link, fluxes = [], [], []
    rows = []
    max_exact = max_total = max_bond_energy = max_meter = max_closed = 0.
    for port, sign, start in states(sector):
        coefficients = vectors.T @ start
        closed_coefficients = closed_vectors.T @ start

        def flight(t):
            return vectors @ (np.exp(-1j * t * energy) * coefficients)

        initial = expectation(start, full)
        initial_cells = expectation(start, closed)
        initial_central = meter.read(start)
        for t in TIMES:
            state = flight(t)
            local = meter.read(state)
            sample = [[meter.read(b.expose(state, a)) for a in ANGLES] for b in bonds]
            independent = np.asarray(sample) @ WEIGHTS
            slopes = [(meter.read(flight(t + dt)) - meter.read(flight(t - dt))) / (2 * dt) for dt in STEPS]
            # Control omitting the second gateway keeps the same measured trajectory.
            error = independent.sum() - slopes[-1]
            residuals.append(error)
            one_link.append(independent[0] - slopes[-1])
            fluxes.append(slopes[-1])
            # d<E>/dt = -2 Im <H psi | E psi>.
            direct = -2 * np.vdot(full @ state, cell_h[0] @ state).imag
            max_exact = max(max_exact, abs(independent.sum() - direct))
            max_meter = max(max_meter, abs(local - expectation(state, cell_h[0])))
            link_energy = [expectation(state, b.h) for b in bonds]
            region_energy = [expectation(state, h) for h in cell_h]
            max_total = max(max_total, abs(sum(region_energy) + sum(link_energy) - initial))
            max_bond_energy = max(max_bond_energy, abs(sum(region_energy) - initial_cells))
            isolated = closed_vectors @ (np.exp(-1j * t * closed_energy) * closed_coefficients)
            max_closed = max(max_closed, abs(meter.read(isolated) - initial_central))
            rows.append(dict(current_port=port, sign=sign, time=t, cell_energies=region_energy,
                             bridge_energies=link_energy, exposed_cell_energy=sample,
                             measured_inflow=independent.tolist(), free_slopes=slopes,
                             analytic_slope=float(direct), residual=float(error)))
    def rms(values):
        return float(np.sqrt(np.mean(np.square(values))))

    summary = dict(rms_free_slope=rms(fluxes), rms_two_gateway_error=rms(residuals),
                   rms_omit_second_gateway_error=rms(one_link), max_two_gateway_error=max(abs(np.array(residuals))),
                   max_two_gateway_commutator_error=float(max_exact), max_meter_error=max_meter,
                   max_closed_cell_drift=max_closed, max_total_energy_drift=max_total,
                   max_sum_of_cell_energy_change=max_bond_energy,
                   max_halved_step_change=max(abs(r["free_slopes"][0]-r["free_slopes"][1]) for r in rows))
    return dict(cells=cells, gateway_ports=links, kappa=KAPPA, sector_dimension=sector.P,
                native_binary_column_error=check_native_gate(), times=TIMES, steps=STEPS,
                angles=ANGLES.tolist(), derivative_weights=WEIGHTS.tolist(), rows=rows, summary=summary)


def portable(tree):
    if isinstance(tree, dict):
        return {k: portable(v) for k, v in tree.items()}
    if isinstance(tree, (tuple, list)):
        return [portable(v) for v in tree]
    if isinstance(tree, (np.floating, float)):
        return round(float(tree), 9) + 0.
    return tree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write the reference measurement")
    args = parser.parse_args()
    report = portable(measure())
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    else:
        # Compare archived data numerically because centered differences amplify eigensolver roundoff.
        reference = json.loads(path.read_text())
        if reference.keys() != report.keys() or reference["cells"] != report["cells"]:
            raise ValueError("registered hardware differs")
        for actual, expected in zip(report["rows"], reference["rows"], strict=True):
            if actual.keys() != expected.keys():
                raise ValueError("registered fields differ")
            for key in actual:
                np.testing.assert_allclose(actual[key], expected[key], rtol=0, atol=2e-8)
        for key in set(report) - {"rows", "summary"}:
            np.testing.assert_allclose(report[key], reference[key], rtol=0, atol=2e-8)
        if report["summary"].keys() != reference["summary"].keys():
            raise ValueError("registered summary fields differ")
        for key, expected in reference["summary"].items():
            np.testing.assert_allclose(report["summary"][key], expected, rtol=0, atol=2e-8)
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
