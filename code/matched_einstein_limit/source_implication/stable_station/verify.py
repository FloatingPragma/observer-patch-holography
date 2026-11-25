"""Check the selected clock record with a separately evaluated 50-digit Hessian."""

from itertools import combinations
import mpmath as mp
import numpy as np
from scipy.linalg import null_space


def precision_check(corners, records, x, eta):
    with mp.workdps(50):
        corner = [mp.matrix(a.tolist()) for a in corners.reshape(4, 4, 4)]
        means = [list(map(mp.mpf, a.ravel())) for a in records]

        def read(u, s, v):
            w = (1 - u - s + v, s - v, u - v, v)
            r = sum((weight * mat for weight, mat in zip(w, corner)), mp.zeros(4))
            gains = [
                mp.exp(sum(mp.log(r[i, j] / r[j, i]) for j in range(4) if j != i) / 4)
                for i in range(4)
            ]
            sq = mp.matrix(4)
            for i, j in combinations(range(4), 2):
                sq[i, j] = sq[j, i] = r[i, j] * r[j, i] / (gains[i] * gains[j])

            def content(sites):
                i, *other = sites
                g = mp.matrix(
                    [
                        [(sq[i, j] + sq[i, k] - sq[j, k]) / 2 for j in other]
                        for k in other
                    ]
                )
                return mp.sqrt(mp.det(g)) / mp.factorial(len(other))

            entropy = mp.log(len(means))
            for rates in means:
                m = mp.fdot(w, rates) - rates[0]
                entropy -= sum(p * mp.log(p) for p in ((1 + m) / 2, (1 - m) / 2)) / len(
                    means
                )
            return (
                sum(content(face) for face in combinations(range(4), 3)),
                content(range(4)),
                entropy,
            )

        point = tuple(map(mp.mpf, x))

        def volume(*args):
            return read(*args)[1]

        def total(*args):
            a, _, k = read(*args)
            return k + eta * a

        ident = np.eye(3, dtype=int)

        def gradient(f):
            return np.array(
                [float(mp.diff(f, point, tuple(map(int, e)))) for e in ident]
            )

        gv, gf = gradient(volume), gradient(total)
        lam = gf @ gv / (gv @ gv)

        def lagrangian(*args):
            a, v, k = read(*args)
            return eta * a + k - lam * v

        h = np.array(
            [
                [
                    float(mp.diff(lagrangian, point, tuple(map(int, i + j))))
                    for i in ident
                ]
                for j in ident
            ]
        )
        basis = null_space(gv[None, :])
        ew, ev = np.linalg.eigh(basis.T @ h @ basis)
        target = volume(*point)
        finite = []
        for direction in (basis @ ev).T:
            for step in (-1e-4, 1e-4):
                y = [a + step * b for a, b in zip(point, direction)]
                # Retain both new independent mixture coordinates. Solve for the
                # correlation probability alone using only the measured volume.
                y[2] = mp.findroot(lambda v: volume(y[0], y[1], v) - target, y[2])
                u, s, v = y
                assert 0 < v / s < 0.5 and 0 < (u - v) / (1 - s) < 0.5
                finite.append(
                    dict(
                        step=step,
                        mixture=list(map(float, y)),
                        volume_error=float(volume(*y) - target),
                        objective_change=float(total(*y) - total(*point)),
                    )
                )
        return dict(
            tangent_gradient_norm=float(np.linalg.norm(basis.T @ gf)),
            lambda_volume=float(lam),
            eigenvalues=ew.tolist(),
            area_volume_entropy=list(map(float, read(*point))),
            equal_volume_secants=finite,
        )
