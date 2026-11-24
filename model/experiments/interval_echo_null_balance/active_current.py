"""Recover hidden circulation from exterior occupation after one local energy pulse."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as contrast
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from .hidden_current import compile_wave, exact_modes
from .regional import verify


EDGE = (4, 5)
ANGLE = np.pi / 3
TIMES = (1/8, 1/4, 1/2, 1.)
RECEIPT = Path(__file__).with_name("active_current.json")


def leading_coefficients(board, bridges, u, v):
    """Signed t^4 coefficients from the exact two-hop amplitude at each receiver."""
    wave = (u + sp.I * v) / sp.sqrt(2)
    local = sp.Matrix([wave[p] for p in EDGE])
    even_phase = (-1 - sp.I * sp.sqrt(3)) / 2  # exp(-8 i pi/3)
    odd_phase = -1  # exp(-9 i pi/3)
    exchange = sp.Matrix([[1, -1], [-1, 1]]) / 2
    kick = (even_phase - 1) * (sp.eye(2) - exchange) + (odd_phase - 1) * exchange
    plus, minus = kick * local, kick * local.conjugate()
    differences = []
    for gateway, receiver in bridges:
        if gateway in EDGE or set(board.neighbours[receiver]).intersection(range(12)) != {gateway}:
            raise ValueError("receiver requires a unique home gateway outside the marked seam")
        incident = [i for i, port in enumerate(EDGE) if port in board.neighbours[gateway]]
        # Each two-hop product has weight 1/4 and the Taylor factorial supplies 1/2.
        first = sum((plus[i] for i in incident), sp.S.Zero) / 8
        second = sum((minus[i] for i in incident), sp.S.Zero) / 8
        differences.append(sp.simplify(abs(first)**2 - abs(second)**2))
    return differences


def pulsed_records(board, initial, bridges, strength, angle):
    """Local pulse, unrestricted flight and a sharp three-port occupation read."""
    exchange, squared_contrast = candidate.native_law.local_term_parts()
    term = exchange + strength * squared_contrast
    ids = [candidate.native_law.SRC.index(pair) for pair in ((1, 0), (0, 1))]
    restricted = term[np.ix_(ids, ids)]
    columns = np.eye(len(term))[:, ids]
    np.testing.assert_allclose(term @ columns, columns @ restricted, rtol=0, atol=1e-13)
    launched = initial.copy()
    launched[list(EDGE)] = expm(-1j * angle * restricted) @ initial[list(EDGE)]
    # Remove one global phase from the native generator, including its degree references.
    free = contrast.generator(board, strength) - 6 * strength * eye(board.n, format="csr")
    np.testing.assert_allclose(free @ initial, 3 * initial, rtol=0, atol=2e-13)
    ports = [receiver for _, receiver in bridges]
    readings = []
    for t in TIMES:
        flight = expm_multiply(-1j * t * free, launched)
        probability = abs(flight[ports])**2
        other = np.ones(board.n, dtype=bool)
        other[ports] = False
        distribution = np.vstack([probability, (abs(flight[other])**2).sum(axis=0)])
        np.testing.assert_allclose(distribution.sum(axis=0), 1., rtol=0, atol=3e-13)
        readings.append(distribution)
    return np.array(readings)


def run():
    local = cell()
    boards = []
    for level in range(3):
        board = federation(level)
        bridges = sorted((a, b) for a, b in board.external if a < 12)
        if len(bridges) != 3 or candidate.KAPPA != 8:
            raise ValueError("the exact coefficient uses the registered three gateways and contact")
        u, v = exact_modes(local, [a for a, _ in bridges])
        coefficients = leading_coefficients(board, bridges, u, v)
        if sorted(coefficients) != [-sp.Rational(1, 128), 0, sp.Rational(1, 128)]:
            raise ArithmeticError("whole-energy pulse lost the exact exterior current contrast")
        psi = np.asarray((u + sp.I*v) / sp.sqrt(2), dtype=complex).ravel()
        initial = np.zeros((board.n, 2), dtype=complex)
        initial[:12] = np.column_stack([compile_wave(local, w) for w in (psi, psi.conj())])
        arms = {}
        for arm, strength, angle in (("whole_energy", 8, ANGLE), ("passive", 8, 0),
                                     ("exchange_only", 0, ANGLE)):
            records = pulsed_records(board, initial, bridges, strength, angle)
            signed = records[:, :, 0] - records[:, :, 1]
            if arm == "whole_energy":
                if not (np.all(signed[:, 0] < -1e-6) and np.all(signed[:, 1] > 1e-6)):
                    raise ArithmeticError("the gateway ports do not resolve the predicted signs")
            elif np.max(abs(signed)) > 1e-13:
                raise ArithmeticError("a control distinguishes the conjugate preparations")
            arms[arm] = {"plus": records[:, :, 0].tolist(), "minus": records[:, :, 1].tolist()}
        boards.append(dict(cells=board.cells, bridges=[list(pair) for pair in bridges],
                           signed_t4=[str(c) for c in coefficients],
                           arms=arms))
    return dict(edge=list(EDGE), flight_exposures=list(TIMES), boards=boards)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        RECEIPT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        verify(json.loads(RECEIPT.read_text()), result)
    for row in result["boards"]:
        print(row["cells"], row["signed_t4"])
