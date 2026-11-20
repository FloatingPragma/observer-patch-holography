"""Compare every registered Fourier fiber and passive tail to the pinned free update."""
import numpy as np

from simulator import law
from . import clock, kernel
from .local_clock import source


def verify_fourier(protocol, local):
    """Check every slot basis vector at each enrolled floating-point wave vector."""
    wave_numbers = set()

    def enroll(angles):
        wave_numbers.update(map(tuple, angles.reshape(-1, 4)))

    def enroll_packet(size, center):
        modes, _, _, _ = kernel.packet(size, center, protocol["halfwidth"])
        enroll(kernel.packet_angles(size, modes))

    for size in protocol["sizes"]:
        for center in protocol["panels"].values():
            enroll_packet(size, clock.center_for(size, center))
    for size in local["sizes"]:
        enroll_packet(size, (16, 0))
        for center in local["center"].values():
            moving = clock.center_for(size, center)
            rest = (0, moving[1])
            enroll_packet(size, moving)
            enroll_packet(size, rest)
            for preparation, tags in ((rest, local["tag_modulus"]),
                                      (moving, local["tag_modulus"]), (moving, 1)):
                _, angles, _ = source.launch(size, preparation, tags)
                enroll(angles)

    ordered = np.array(sorted(wave_numbers))
    basis = np.eye(law.SLOTS, dtype=complex)
    maximum = 0.0
    for start in range(0, len(ordered), 256):
        angles = ordered[start:start + 256]
        # The input row indexes a basis vector. walk_matrix stores basis images in columns.
        actual = kernel.step(basis, kernel.multipliers(angles)[:, None, :])
        expected = np.stack([law.walk_matrix(k, kernel.FREE.equal_phase).T for k in angles])
        error = float(np.max(np.abs(actual - expected)))
        kernel.require(np.isfinite(error), "nonfinite Fourier action")
        maximum = max(maximum, error)
    kernel.require(maximum <= 2e-14, "frozen Fourier action differs from the free walk")
    return {"fibers": len(ordered), "maximum_column_error": maximum}


def passive_outcomes(waits, speeds, tails):
    """All four terminal outcomes on the strip, with exact integer Born weights.

    At each address, 4 C has diagonal -3 and off-diagonal 1. Carry its
    unscaled integer amplitudes through the literal simulator route. The
    normalized initial amplitudes have a common factor 1/sqrt(8). The
    equal-slot detector contributes a further factor 1/sqrt(8).
    """
    slots = law.SLOTS
    matrix = 4 * law.coin(np.eye(slots), kernel.FREE.equal_phase)
    kernel.require(np.array_equal(matrix, np.ones((slots, slots)) - 4 * np.eye(slots)),
                   "integer lift requires the frozen equal-slot reflection")
    origin = max(waits)
    state = np.zeros((2 * origin + 1, 1, 1, 1, slots), dtype=object)
    state[origin] = 1
    result = {}
    for t in range(1, origin + 1):
        state = law.route(state.sum(axis=-1, keepdims=True) - 4 * state)
        if t not in waits:
            continue
        amps = state[:, 0, 0, 0, :]
        total = slots * (amps * amps).sum(axis=1)
        equal = amps.sum(axis=1) ** 2
        other = total - equal
        denominator = slots ** 2 * 16 ** t
        kernel.require(total.sum() == denominator and min(other) >= 0,
                       "passive detector partition lost normalization")
        rows = {}
        for a, b in speeds:
            key = f"{a}/{b}"
            cut = origin + (a * t + b - 1) // b
            counts = [equal[cut:].sum(), equal[:cut].sum(), other[cut:].sum(), other[:cut].sum()]
            reviewed = int(tails[str(t)][key]["numerator"])
            kernel.require(counts[0] == 16 * reviewed, f"passive tail differs at {t}, {key}")
            rows[key] = dict(zip(("equal_at_or_above", "equal_below",
                                  "orthogonal_at_or_above", "orthogonal_below"), map(str, counts)))
        result[str(t)] = {"denominator": str(denominator), "thresholds": rows}
    return result
