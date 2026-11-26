import json
from fractions import Fraction

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ..measure import fixture, material
from . import certificate, experiment


def test_source_gates_native_contacts_and_swap_symmetry():
    board = federation(0)
    dimers, ports, links = fixture(board, 2)
    assert tuple(dimers) == certificate.DIMERS
    assert tuple(ports) == certificate.CONTACTS
    assert {82 // 12, 74 // 12} == {6}
    # Addressed exposure pi/(4*kappa) is CZ up to an overall phase
    # in the endpoint occupation basis. Blank dimer exchange plus this
    # gate with one neighbor flip yields both one-site rotation axes.
    law = candidate.native_law
    u = np.eye(law.DIM)[:, : law.N_SRC]
    donor = (u[:, law.SRC.index((2, 0))] + u[:, law.SRC.index((0, 2))]) / np.sqrt(2)
    collar = (law.echo() @ donor).reshape(-1, law.N_SRC)[:, law.SRC.index((1, 1))]
    exchange, contrast = law.local_term_parts()
    for a in (0, 1):
        for b in (0, 1):
            word = np.kron(collar, np.eye(law.N_SRC)[:, law.SRC.index((a, b))])
            energy = 4 * candidate.KAPPA * (1 - a * b)
            np.testing.assert_allclose(
                (exchange + candidate.KAPPA * contrast) @ word,
                energy * word,
                atol=1e-13,
            )
    neighbors = {frozenset((a, b)) for a, _, b, _ in links}
    assert all(
        frozenset(edge) in neighbors
        for path in experiment.PATHS
        for edge in zip(path, path[1:])
    )
    swap = np.arange(128)
    for b, c in ((3, 5), (4, 6)):
        swap = experiment.swap(swap, b, c)
    for coupling in experiment.CONTRAST:
        h, wall = material(links, coupling)
        np.testing.assert_array_equal(h.toarray()[np.ix_(swap, swap)], h.toarray())
        assert all(swap[w] == w for w in np.flatnonzero(np.isin(wall, [1, 2])))
        for r in (0, 1):
            src, cal = experiment.sources(wall, r)
            np.testing.assert_array_equal(src[swap, 2], src[:, 2])
            np.testing.assert_array_equal(src[swap, 3], -src[:, 3])
            np.testing.assert_allclose(cal[2:, 0], np.eye(2), atol=1e-14)
            np.testing.assert_allclose(cal[:, 1], 0.5, atol=1e-14)


def test_integer_native_six_radian_witness():
    with np.load(experiment.HERE / "laws.npz") as stored:
        cdf = stored["local"]
        t = int(np.flatnonzero(experiment.TIMES == 6)[0])
        for j, phase in ((2, 1), (3, -1)):
            lower, upper = certificate.remote(phase)
            prob = Fraction(
                int(cdf[0, 1, j, t, 2] - cdf[0, 1, j, t, 0]), experiment.SCALE
            )
            assert (
                lower - Fraction(1, experiment.SCALE)
                < prob
                < upper + Fraction(1, experiment.SCALE)
            )
            assert upper - lower < Fraction(1, 10**15)
            if j == 2:
                assert Fraction(59201, 100000) < lower < upper < Fraction(59202, 100000)
            else:
                assert Fraction(11383, 100000) < lower < upper < Fraction(11384, 100000)


def test_keeping_all_records_and_discriminating_global_from_local():
    with np.load(experiment.HERE / "records.npz") as records:
        for cohort in records.values():
            assert np.all(cohort.sum(axis=-1) == experiment.SHOTS)
        report = experiment.read(records)
    stored = json.loads((experiment.HERE / "results.json").read_text())
    assert stored["settings"]["k8_r1"]["remote"]["plus"] == 158926139
    for r in (0, 1):
        row = report["settings"][f"k8_r{r}"]
        assert row["remote"]["earliest"] is not None
        assert row["scalar_early"]["forecast"][0] > row["scalar_early"]["observed"][1]
        control = report["settings"][f"k0_r{r}"]
        assert control["remote"]["earliest"] is None
        assert control["near"]["earliest"] is not None
