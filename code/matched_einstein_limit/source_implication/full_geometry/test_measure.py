import json
from pathlib import Path

import numpy as np
import pytest

from .measure import fit, run


def test_fifth_port_and_anchored_gram():
    data = run()
    saved = json.loads(Path(__file__).with_name("results.json").read_text())
    np.testing.assert_allclose(data["range_matrix"],saved["range_matrix"],atol=2e-11,rtol=0)
    extended = data["extended_five"]
    assert extended["last_to_first"] == pytest.approx(.66437169,abs=1e-7)
    d2 = np.asarray(extended["distances"])**2
    # Difference vectors from the new station provide an independent, uncentered Euclidean test.
    anchor = (d2[:-1,-1,None] + d2[None,-1,:-1] - d2[:-1,:-1])/2
    assert np.linalg.eigvalsh(anchor).min() > 3
    assert np.prod(extended["gains"]) == pytest.approx(1.)


def test_genuine_three_dimensional_survey_passes():
    points = np.array([[0,0,0],[1,0,0],[0,2,0],[0,0,3],[1,1,1]])
    distances = np.linalg.norm(points[:,None]-points,axis=2)
    rates = np.array([.8,1.1,1.2,.95,1.3])
    fitted = fit(rates[:,None]*distances)
    assert abs(fitted["last_to_first"]) < 1e-14
    np.testing.assert_allclose(fitted["gains"],rates/np.prod(rates)**.2)
