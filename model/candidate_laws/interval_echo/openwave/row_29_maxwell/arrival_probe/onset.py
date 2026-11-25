"""Exact derivative test for an independently injected double on the receiving seam."""

import json

import numpy as np

from model.candidate_laws.interval_echo import clock
from model.experiments.material_inertia.neutral_defect import generators
from .. import experiment as dipole
from .measure import HERE


def audit():
    rec = json.loads((HERE / "records.json").read_text())
    cell = dipole.carrier.cell()
    x, y, adjacent, exchange, _, _ = generators(cell)
    index = {(p, q): i for i, (p, q) in enumerate(zip(x, y, strict=True))}
    # 2H is integral in both settings. In the one-double sector 2H=L
    # up to a constant. L is the local site Laplacian.
    single = cell.laplacian().toarray().astype(np.int64)
    for kappa in (0, clock.KAPPA):
        twice = 2 * exchange.toarray() - np.diag(8 * kappa * adjacent)
        h = twice.astype(np.int64)
        np.testing.assert_array_equal(h, twice)
        for port in range(12):
            p = (x == port).astype(np.int64)
            comm = h * (p[None, :] - p[:, None])
            acc = comm @ h - h @ comm  # 4 times the Hermitian acceleration.
            q = (np.arange(12) == port).astype(np.int64)
            solo = single * (q[None, :] - q[:, None])
            for a, b in cell.seams:
                kept = [s for s in range(12) if s not in (a, b)]
                rows = [index[a, s] for s in kept]
                cols = [index[b, s] for s in kept]
                # Every immediate current is identical to the empty-receiver
                # instrument for these launches and any admitted arrival.
                for u, i in ((a, rows), (b, cols)):
                    for v, j in ((a, rows), (b, cols)):
                        np.testing.assert_array_equal(
                            comm[np.ix_(i, j)], solo[u, v] * np.eye(10, dtype=np.int64)
                        )
                cross = acc[np.ix_(rows, cols)]
                # Reversing the input phase gives i(cross-cross^T).
                np.testing.assert_array_equal(cross, cross.T)
    # Positive detector control. Insert a phase on a different side of a
    # native triangle. The third-port counter can see different accelerations
    # for the two launch phases even though its initial population vanishes.
    a, b = rec["probe_ports"]
    c = min(set(cell.neighbours[a]) & set(cell.neighbours[b]))
    twisted = single.astype(complex) / 2
    twisted[c, a], twisted[a, c] = -0.5j, 0.5j
    initial_acc = [abs(twisted[c, a] + phase * twisted[c, b])**2 for phase in (1j, -1j)]
    np.testing.assert_array_equal(initial_acc, [1, 0])
    # Each face outcome's unit-normal contribution to a component is
    # bounded by sqrt(3)/2. Simultaneous Hoeffding over 8*20*3 means.
    delta = np.sqrt(3 * np.log(2 * 8 * 20 * 3 / 0.01) / (2 * rec["shots"]))
    u = np.array(rec["direction"])
    for name in ("circle+", "circle-"):
        m = np.array(rec["matched_axial_reading"][name])
        print(name, "u cross M =", np.cross(u, m), "norm error <=", np.sqrt(3) * delta)
        assert np.linalg.norm(np.cross(u, m)) > np.sqrt(3) * delta
    print(
        "Both laws, all 30 seams: every double-address probability has zero velocity-odd second derivative at injection."
    )
    print("Supplied triangle phase, third-port accelerations:", initial_acc)


if __name__ == "__main__":
    audit()
