"""Independent Ramsey and terminal-position reads of a material-written bridge."""

from functools import lru_cache

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import eigsh, expm_multiply

from model.executable_model import federation
from model.experiments.interval_echo_einstein_station.mobile_collar_source.measure import (
    region,
)
from ...source_implication.material_clocks.measure import blank_clock

OMEGA = 0.5
EPSILON = np.log(1.25)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1])
IDENTITY = np.eye(2)


def coefficients(effect):
    """Trace against (I+xX+yY+zZ)/2 gives the complete population mean."""
    np.testing.assert_allclose(effect, effect.conj().T, atol=4e-13, rtol=0)
    return [float(np.trace(effect @ p).real / 2) for p in (IDENTITY, X, Y, Z)]


@lru_cache(maxsize=3)
def modes(level):
    board = federation(level)
    a, b = region(board)[1:3]
    h = board.laplacian() / 2
    # On blank (0,0) seams a single excess sees 8D + (D-Xport)/2.
    # The physical degree references complete the diagonal offset to 6*8.
    start = np.random.default_rng(7325).normal(size=board.n)
    energies, axes = eigsh(h, k=5, which="SA", v0=start, tol=1e-13)
    slow = axes[:, 1:4]
    weights = slow[b] - slow[a]
    weights /= np.linalg.norm(weights)
    u = slow @ weights
    g = np.ones(board.n) / np.sqrt(board.n)
    np.testing.assert_allclose(slow.T @ g, 0, atol=5e-12)
    np.testing.assert_allclose(h @ slow, slow * energies[1:4], atol=1e-12)
    if not (energies[0] < 1e-12 < energies[1] and energies[3] < energies[4]):
        raise ValueError("the selected three-mode band must be separated")
    energy = float(np.dot(weights**2, energies[1:4]))
    score = u / max(abs(u))
    # The outside rail sees zero contact energy, fixed across source arms.
    bond = sparse.lil_matrix(h.shape)
    bond[np.ix_([a, b], [a, b])] = blank_clock(-1, 8)
    return h, bond.tocsr(), g, u, score, energy, energies[1:4], weights, (a, b)


def reciprocal(h, bond, epsilon, omega):
    size = h.shape[0]
    return (
        sparse.kron(IDENTITY, h, format="csr")
        - sparse.kron(omega * X, sparse.eye(size), format="csr")
        + sparse.kron(np.diag(np.expm1([-epsilon, epsilon])), bond, format="csr")
    )


def local_effect(energy_left, energy_right, time, epsilon, omega):
    """Signed interference quadrature for two real energy sectors."""
    multiplier = np.diag(np.exp([-epsilon, epsilon]))
    left = expm(-1j * time * (energy_left * multiplier - omega * X))
    right = expm(-1j * time * (energy_right * multiplier - omega * X))
    coherence = left.conj().T @ right
    return (coherence.conj().T - coherence) / (2j)


def band_observables(waves, g, u, score):
    """Position-score difference and a separately prepared mode Ramsey quadrature."""
    # Two basis packets are evolved for each of the two initial strain columns.
    waves = waves.reshape(2, len(g), 2, 2)
    const, dipole = waves[:, :, 0], waves[:, :, 1]
    cross = np.einsum("spi,p,spj->ij", const.conj(), score, dipole, optimize=True)
    position = 0.5j * (cross - cross.conj().T)
    real_launch = (const + dipole) / np.sqrt(2)
    left = np.einsum("p,spi->si", g, real_launch)
    right = np.einsum("p,spi->si", u, real_launch)
    cross = left.conj().T @ right
    clock = 1j * (cross - cross.conj().T)
    return np.array([position, clock])


def initial_reads(bond, g, u, score, energy, ends):
    a, b = ends
    # A positive flux is a rise of the signed detector coordinate.
    overlap = float(np.dot(g * score, u))
    position_rate = energy * overlap
    edge_rate = float((score[b] - score[a]) * (g[a] * u[b] - u[a] * g[b]) / 2)
    clock_shift = float(u @ (bond @ u) - g @ (bond @ g))
    return dict(
        position_rate=position_rate,
        clock_rate=energy,
        fraction_of_position_rate_on_bridge=edge_rate / position_rate,
        fractional_clock_response_to_bridge=clock_shift / energy,
        position_score_overlap=overlap,
    )


def panel(level, epsilon=EPSILON, omega=OMEGA):
    h, bond, g, u, score, gap, energies, weights, ends = modes(level)
    extended = reciprocal(h, bond, epsilon, omega)
    gen = -1j * extended
    launch = np.kron(IDENTITY, np.column_stack([g, u])).reshape(2 * len(g), 2, 2)
    # Kronecker ordering above is initial strain then mode. The analyzer expects
    # mode then initial strain in the two trailing axes.
    waves = launch.swapaxes(1, 2).reshape(2 * len(g), 4)
    t_previous = 0
    rows = []
    for time in (1 / 32, 1 / 8, np.pi / (6 * gap), np.pi / (3 * gap)):
        dt = time - t_previous
        waves = expm_multiply(dt * gen, waves, traceA=dt * gen.diagonal().sum())
        np.testing.assert_allclose(
            waves.conj().T @ waves, np.eye(4), atol=8e-12, rtol=0
        )
        values = band_observables(waves, g, u, score)
        # The derivative of a Hermitian quadratic read Q on psi is its cross
        # term between psi and dpsi/dt. Polarization also covers signed effects.
        velocity = gen @ waves
        deriv = (
            band_observables(waves + velocity, g, u, score)
            - band_observables(waves - velocity, g, u, score)
        ) / 2
        row = dict(
            time=float(time),
            band_position=coefficients(values[0]),
            band_clock=coefficients(values[1]),
            band_position_derivative=coefficients(deriv[0]),
            band_clock_derivative=coefficients(deriv[1]),
        )
        if time <= 1 / 8:
            row["outside_bridge_clock"] = coefficients(
                local_effect(0, 8, time, epsilon, omega)
            )
            row["directed_bridge_transfer"] = coefficients(
                local_effect(8, 9, time, epsilon, omega)
            )
        rows.append(row)
        t_previous = time
    return dict(
        cells=len(g) // 12,
        ports=len(g),
        bridge=list(ends),
        slow_gaps=energies.tolist(),
        spectral_probabilities=(weights**2).tolist(),
        max_dipole_cell_probability=float(np.max((u.reshape(-1, 12) ** 2).sum(axis=1))),
        mean_gap=gap,
        rates=initial_reads(bond, g, u, score, gap, ends),
        appointments=rows,
    )
