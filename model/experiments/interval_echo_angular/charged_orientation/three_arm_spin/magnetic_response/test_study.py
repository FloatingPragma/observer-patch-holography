"""Independent checks of the port meter and the physically observable phase classes."""

from collections import deque
from math import gcd

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from . import study
from ..apparatus import fixture
from ..recoil_search.guides import detour


def chart():
    _, ports, legs, h, actions, _, _ = detour(0)
    f = fixture()
    paths = [[f["hub"], *leg, f["tip"]] for leg in legs]
    return ports, h, paths, study.bonds(h, actions, ports, paths, f["marked"])[0]


def test_read_all_ports_after_literal_contact_pulses():
    ports, h, _, bonds = chart()
    rng = np.random.default_rng(51105)
    psi = rng.normal(size=h.shape[0]) + 1j * rng.normal(size=h.shape[0])
    psi /= np.linalg.norm(psi)
    index = {p: 64 * i + np.arange(64) for i, p in enumerate(ports)}
    measured = study.meter(psi, bonds, index, 1.0)["edges"]
    check = [0, len(bonds) - 1]
    for k in check:
        (u, v), sources, destinations, magnitude = bonds[k]
        a = coo_matrix(
            (np.full(len(sources), -magnitude), (destinations, sources)), shape=h.shape
        ).tocsr()
        for pulse, row in zip((1, 3), measured[k]["original_start_laws"], strict=True):
            out = expm_multiply(-1j * pulse * np.pi / (4 * magnitude) * (a + a.T), psi)
            ps, pv = (np.linalg.norm(out[index[p]]) ** 2 for p in (u, v))
            np.testing.assert_allclose(row[:3], [ps / 3, pv / 3, (1 - ps - pv) / 3])
        current = 2 * np.vdot(psi, a @ psi).imag
        assert (
            abs(current / magnitude - measured[k]["conditional_current_score"]) < 1e-12
        )


def test_a_spanning_tree_recovers_zero_and_nonzero_loop_flux():
    _, h, paths, bonds = chart()
    labels = {}
    first_left_edges = {(r[0], r[1]) for r in paths}
    marked = set(fixture()["marked"])
    for edge, source, dest, _ in bonds:
        labels.update({(a, b): edge for a, b in zip(source, dest, strict=True)})
    off = h.copy()
    off.setdiag(0)
    off.eliminate_zeros()

    def check(one_arm):
        potentials = {0: 0}
        search = deque([0])
        loop_phase_unit = 0
        while search:
            src = search.popleft()
            for dest in off.getrow(src).indices:
                orientation = 1 if (src, dest) in labels else -1
                edge = labels[(src, dest) if orientation == 1 else (dest, src)]
                phase = 3 if edge in marked else int(edge in first_left_edges)
                if one_arm and edge == (paths[0][0], paths[0][1]):
                    phase += 2
                wanted = potentials[src] + orientation * phase
                if dest in potentials:
                    loop_phase_unit = gcd(loop_phase_unit, wanted - potentials[dest])
                else:
                    potentials[dest] = wanted
                    search.append(dest)
        assert len(potentials) == h.shape[0]
        return loop_phase_unit

    assert check(False) == 0
    assert check(True) == 2
