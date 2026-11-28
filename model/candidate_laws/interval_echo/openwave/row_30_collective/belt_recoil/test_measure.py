import json
from collections import Counter

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.executable_model import carrier
from . import measure


def archive():
    return json.loads((measure.HERE / "counts.json").read_text())


def test_native_incidence_reproduces_two_held_distributions():
    raw = archive()
    ports = raw["belt_ports"]
    seams = [(ports.index(a), ports.index(b)) for a, b in carrier.cell().seams
             if a in ports and b in ports]
    board = carrier.Carrier("native belt", len(ports), seams)
    for loads, harmonic in ((3, 3), (4, -4)):
        sector = BinarySector(board, loads)
        walls = []
        for start in range(10):
            steps = [0, 4, 7] if loads == 3 else [0, 1, 4, 7]
            word = tuple(sorted((start+j) % 10 for j in steps))
            walls.append(sector.index[word])
        launch = np.zeros(len(sector.states), dtype=complex)
        launch[walls] = (1+np.exp(2j*np.pi*harmonic*np.arange(10)/10)) / np.sqrt(20)
        prob = abs(expm(-2j*sector.hamiltonian(8)) @ launch)**2
        row = next(r for r in raw["rows"] if (r["loads"], r["kappa"], r["harmonic"])
                   == (loads, 8, harmonic))
        panel = np.array(row["counts"][raw["waits"].index(2)]) / raw["shots"]
        assert max(abs(panel - np.r_[prob[walls], 1-sum(prob[walls])])) < .004

        projected = sector.K[np.ix_(walls, walls)]
        assert np.all(np.diag(projected) == projected[0, 0])
        shifts = np.roll(np.eye(10), 3, axis=0) + np.roll(np.eye(10), -3, axis=0)
        assert np.array_equal(projected, projected[0, 0]*np.eye(10) - shifts/2)


def test_planted_positive_and_negative_massive_bands():
    q = 2*np.pi/10
    for sign in (-1, 1):
        curves = []
        for harmonic in (1, 2, 3, 4):
            delta = sign*(np.sqrt(1.2**2 + (q*harmonic)**2) - 1.2)
            x = np.arange(10)
            # Twenty percent loss tests that the reader keeps complement outcomes.
            probs = .8*(1+np.cos(q*harmonic*x[None, :] - np.outer(measure.WAITS, delta))) / 10
            counts = np.rint(probs*measure.SHOTS).astype(int)
            rows = np.column_stack((counts, measure.SHOTS-counts.sum(axis=1)))
            result = measure.recoil(dict(counts=rows, weights=[1, 1], harmonic=harmonic))
            lo, hi = result["interval"]
            assert lo <= delta <= hi
            curves.append((lo, hi))
        for harmonic, (lo, hi) in enumerate(curves[1:], 2):
            minimum, maximum = measure.forecast(curves[0], harmonic, "continuum")
            assert maximum >= lo and minimum <= hi


def test_archived_comparisons_and_counts_integrity():
    raw = archive()
    rows = measure.read(raw)
    for kind in ("continuum", "lattice"):
        verdicts = Counter((r["kappa"], r[kind]["status"]) for r in rows if kind in r)
        assert verdicts == {(8, "excluded"): 12, (8, "compatible"): 2, (0, "unresolved"): 14}
    row = raw["rows"].pop()
    with pytest.raises(ValueError, match="Incomplete launch"):
        measure.read(raw)
    row["counts"][1][-1] += 1
    with pytest.raises(ValueError, match="Incomplete wall"):
        measure.recoil(row)
