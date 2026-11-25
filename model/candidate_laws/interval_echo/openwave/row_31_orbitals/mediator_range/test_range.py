"""Independent finite checks of the ground certificate, range and clock records."""

from fractions import Fraction
import json
from math import cos, pi

import numpy as np
import pytest
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from . import measure as m


def assembled_native(board, written, amount):
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    h = h0 + candidate.KAPPA * h1
    blank = np.eye(law.DIM)[:, :9]
    vacuum = law.SRC.index((0, 0))
    memory = (law.echo() @ blank[:, vacuum]).reshape(11, 9)[:, vacuum]
    written_loads = np.kron(memory[:, None], np.eye(9))
    # Build each column, including the registered compensation for missing neighbors.
    matrix = np.diag(candidate.KAPPA * (6 - board.degree) * amount**2).astype(float)
    for edge in board.seams:
        a, b = edge
        basis = written_loads if edge in written else blank
        generator = basis.T @ h @ basis
        np.testing.assert_allclose(h @ basis, basis @ generator, atol=2e-13)
        sub = [law.SRC.index((amount, 0)), law.SRC.index((0, amount))]
        matrix[np.ix_([a, b], [a, b])] += generator[np.ix_(sub, sub)]
    return matrix - 6 * candidate.KAPPA * amount**2 * np.eye(board.n)


@pytest.mark.parametrize("amount", (1, 2))
def test_full_native_cell_ground_and_range(amount):
    board = cell()
    lap = board.laplacian().astype(np.int64)
    distance = shortest_path(board.adjacency(), directed=False, unweighted=True)
    first = tuple(e for e in board.seams if 0 in e)[:2]
    available = [e for e in board.seams if set(e).isdisjoint(m.supports(first))]
    for second in ((available[0],), tuple(available[:3])):
        d = int(distance[np.ix_(m.supports(first), m.supports(second))].min())
        levels = []
        for written in (first, second, first + second):
            values = np.linalg.eigvalsh(assembled_native(board, written, amount))
            interval = m.binding_interval(lap, written, candidate.KAPPA * amount**2)
            low, high = interval["binding_nanounits"]
            assert low / m.SCALE < -values[0] < high / m.SCALE
            levels.append(values[0])
        gain = min(levels[:2]) - levels[2]
        assert (
            -1e-12
            <= gain
            <= m.envelope(d, m.max_touch(first), m.max_touch(second), amount)
        )


def test_wrong_ground_bracket_fails_on_both_sides():
    board = cell()
    lap = board.laplacian().astype(np.int64)
    edges = (board.seams[0],)
    moments = m.local_powers(lap, m.supports(edges))
    w2 = m.defect_twice(edges, 8).astype(object)
    for t, upper in ((Fraction(4), True), (Fraction(9), False)):
        with pytest.raises(
            ArithmeticError, match="inconclusive ground-energy certificate"
        ):
            m.certify_side(t, w2, moments, upper)


def test_deposited_clocks_and_calibration_against_universal_cutoff():
    records = json.loads((m.HERE / "results.json").read_text())
    board = federation(2)
    for shape in records["shapes"].values():
        rows = shape["distance_panels"]
        first = shape["first_edges"]
        dist = shortest_path(
            board.adjacency(),
            directed=False,
            unweighted=True,
            indices=m.supports(first),
        ).min(axis=0)
        cal = shape["calibration"]
        for d, row in rows.items():
            assert min(dist[m.supports(row["second_edges"])]) == int(d)
            for bank in [
                v for v in row.values() if isinstance(v, dict) and "plus" in v
            ]:
                assert bank["plus"] + bank["minus"] == records["shots_per_bank"]
        assert rows["1"]["dark_bright"]["plus"] > 0
        for d in ("6", "9"):
            row = rows[d]
            assert row["dark_bright"]["plus"] == 0
            assert row["zero_contrast_dark_bright"]["plus"] == 0
            assert row["receiver_pi_zero_contrast"]["minus"] == 0
            assert (
                row["coulomb_dark_floor"]
                == cos(pi * cal["allowance"] / (2 * cal["gain"])) ** 2
            )
            assert row["coulomb_dark_floor"] > 0.9988
        degree = m.max_touch(first)
        cutoff = 7 if degree == 1 else 8
        coefficient = 3 * (cal["gain"] - cal["allowance"])
        assert m.envelope(cutoff, degree, degree, 1) < coefficient / (cutoff + 1)
        assert m.envelope(cutoff - 1, degree, degree, 1) >= coefficient / cutoff
