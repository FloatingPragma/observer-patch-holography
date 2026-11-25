"""Joint rest, chirality, and charge read in the declared finite free-field realization."""

import numpy as np


def transfer(a, b):
    """The exterior action c†_a c_b on the five binary occupations."""
    result = np.zeros((32, 32))
    for source in range(32):
        if not source & (1 << b):
            continue
        removed = source ^ (1 << b)
        if removed & (1 << a):
            continue
        sign = (-1) ** (
            (source & ((1 << b) - 1)).bit_count()
            + (removed & ((1 << a) - 1)).bit_count()
        )
        result[removed | (1 << a), source] = sign
    return result


def port_chart(projector):
    """Fix the coordinate basis by the first three independent physical port columns."""
    columns = []
    for column in projector.real.T:
        direction = column.copy()
        for earlier in columns:
            direction -= np.dot(earlier, direction) * earlier
        size = np.linalg.norm(direction)
        if size > 1e-9:
            columns.append(direction / size)
    basis = np.column_stack(columns)
    assert basis.shape == (12, 3)
    np.testing.assert_allclose(basis @ basis.T, projector, atol=2e-14)
    return basis


def rest_and_motion(projectors, masses):
    # Use the positive-rank part of the theorem band, excluding its nine null ports.
    p = projectors.sum(axis=0)
    v = port_chart(p)
    mass = v.conj().T @ np.einsum("j,jkl->kl", masses, projectors) @ v
    occupation = np.stack([np.diag(transfer(a, a)) for a in range(5)])
    q = -occupation[:3].sum(axis=0) / 3 + occupation[3]  # Y + T3
    assert [q[k] for k in (23, 15, 24)] == [-1, 0, 1]
    raising = transfer(3, 4)
    assert abs(raising[15, 23]) == 1 and np.count_nonzero(raising[:, 24]) == 0
    # The two Cartan generators and the six off-diagonal colour actions.
    color = [transfer(a, b) for a in range(3) for b in range(3) if a != b]
    color += [transfer(0, 0) - transfer(1, 1), transfer(1, 1) - transfer(2, 2)]
    for g in color:
        np.testing.assert_array_equal(g[:, [15, 23, 24]], 0)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1])
    eye = np.eye(2)
    beta = np.kron(X, eye)
    alpha = [np.kron(-Z, axis) for axis in (X, Y, Z)]
    rows = []
    for sign, M in [(-1, mass), (1, mass.conj())]:
        for momentum in ([0, 0, 0], [0.25, 0, 0], [0, 0.25, 0], [0, 0.15, 0.2]):
            h = np.kron(beta, M) + np.kron(
                sum(k * a for k, a in zip(momentum, alpha, strict=True)), np.eye(3)
            )
            sq = np.kron(np.eye(4), M @ M) + np.dot(momentum, momentum) * np.eye(12)
            np.testing.assert_allclose(h @ h, sq, atol=2e-13)
            levels = np.linalg.eigvalsh(h)
            positive = levels[levels > 0]
            centers = np.sqrt(np.sort(masses) ** 2 + np.dot(momentum, momentum))
            np.testing.assert_allclose(positive, np.repeat(centers, 2), atol=2e-12)
            axis = np.array(momentum, dtype=float)
            if not np.any(axis):
                axis = np.array([0.0, 1.0, 0.0])
            axis /= np.linalg.norm(axis)
            spin = np.kron(
                np.kron(
                    eye, sum(k * a / 2 for k, a in zip(axis, (X, Y, Z), strict=True))
                ),
                np.eye(3),
            )
            np.testing.assert_allclose(h @ spin, spin @ h, atol=1e-13)
            for energy in centers:
                vals, cols = np.linalg.eigh(h)
                mode = cols[:, abs(vals - energy) < 1e-10]
                np.testing.assert_allclose(
                    np.linalg.eigvalsh(mode.conj().T @ spin @ mode),
                    [-0.5, 0.5],
                    atol=5e-13,
                )
            rows.append(
                dict(
                    charge=sign,
                    momentum=momentum,
                    positive_energies=centers.tolist(),
                    spin_values=[-0.5, 0.5],
                )
            )
    return dict(
        mask_e_left=23,
        mask_nu_left=15,
        mask_conjugate_right=24,
        colour_generators_checked=len(color),
        field_banks=rows,
        declared_inputs=[
            "lepton exterior attachment",
            "conjugate right field and two-component spin",
            "Dirac kinetic and mass terms",
            "fixed family orientation per history",
        ],
    )
