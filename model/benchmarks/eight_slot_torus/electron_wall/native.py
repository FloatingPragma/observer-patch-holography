"""Probe the torus walk and both fully mature compact rails at one address."""
import numpy as np

from model.benchmarks.eight_slot_torus import law


def half_turn(state, axis):
    """Rotate the public carrier and its slots by pi, fixing the private register."""
    out = state
    slots = np.arange(8)
    for other in range(3):
        if axis == other:
            continue
        out = np.take(out, (-np.arange(out.shape[other])) % out.shape[other], axis=other)
        slots[2 * other:2 * other + 2] = slots[2 * other:2 * other + 2][::-1]
    return out[..., slots]


def native_rows(period):
    if period < 3:
        raise ValueError("each slot must have a distinct neighbour")
    parent = law.LAWS["torus"]
    updates = [lambda state: law.step(state, parent)]
    for rail in range(2):
        def mature(state, index=rail):
            return law.step_maturing(state, parent, age=parent.maturation[-1], rail=index)
        updates.append(mature)
    result = []
    for update in updates:
        coin = np.empty((8, 8), complex)
        # Source translation generates the other columns of the full walk.
        error = 0.
        for i in range(8):
            state = np.zeros((period,) * 4 + (8,), complex)
            state[0, 0, 0, 0, i] = 1
            output = update(state)
            reconstructed = np.zeros_like(state)
            for slot in range(8):
                where = [0] * 4
                where[slot // 2] = (1 if slot % 2 == 0 else -1) % period
                index = (*where, slot ^ 1)
                coin[slot, i] = output[index]
                reconstructed[index] = output[index]
            error = max(error, float(np.max(np.abs(output - reconstructed))))
            for axis in (0, 1):
                error = max(error, float(np.max(np.abs(half_turn(output, axis)
                                      - update(half_turn(state, axis))))))
            commutator = half_turn(half_turn(half_turn(half_turn(output, 0), 1), 0), 1)
            error = max(error, float(np.max(np.abs(commutator - output))))
        if error > 1e-12 or np.any(np.abs(coin) < 1e-8):
            raise ArithmeticError("native one-slot address transfers or geometric rotation failed")
        result.append(dict(min_entry=float(np.abs(coin).min()), offdiag=float(abs(coin[0, 1])),
                           residual=error))
    if abs(result[1]["offdiag"] - result[2]["offdiag"]) < 1e-8:
        raise ArithmeticError("compact rail inequivalence failed")
    return result
