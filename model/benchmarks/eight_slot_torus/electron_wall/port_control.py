"""An eight-socket control with the native free coin and Peierls-Pauli routing."""
import numpy as np

from model.benchmarks.eight_slot_torus import law
from .walk import PAULI


def step(state, connection, charge):
    if state.ndim != 6 or state.shape[-2:] != (8, 2) or charge not in (-1, 1):
        raise ValueError("expected four native axes, eight sockets and a spinor")
    if connection.shape != (3, *state.shape[:3]) or not np.isfinite(connection).all():
        raise ValueError("invalid public connection")
    if np.max(abs(abs(connection) - 1)) > 1e-12:
        raise ValueError("nonunitary public link")
    coin = law.coin(state, axis=-2)
    for axis in range(3):
        phase = connection[axis] ** charge
        for direction, sign in enumerate((1, -1)):
            port = 2 * axis + direction
            oriented = phase if sign == 1 else np.roll(phase.conj(), 1, axis=axis)
            coin[..., port, :] = ((sign * 1j * oriented[..., None, None])
                                  * (coin[..., port, :] @ PAULI[axis].T))
    return law.route(coin, slot_axis=-2)
