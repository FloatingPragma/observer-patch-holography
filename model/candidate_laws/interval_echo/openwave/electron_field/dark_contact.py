"""Check the stored excitation against all nine native two-port load states."""
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo.candidate import KAPPA


def lift(contrast, contrast_parity, coin_parity):
    memory = np.zeros(1 + law.N_LAB, dtype=int)
    for coin in (0, 1):
        for signed in sorted({contrast, -contrast}):
            amplitude = (coin_parity if coin else 1) * (contrast_parity if signed < 0 else 1)
            memory[1 + law.LABELS.index((coin, signed))] = amplitude
    return np.kron(memory[:, None], np.eye(len(law.SRC), dtype=int))


def run():
    exchange, contrast_cost = law.local_term_parts()
    local = exchange + KAPPA * contrast_cost
    duration = np.pi / (4 * KAPPA)
    pulse = expm(-1j * duration * local)
    parked = law.SRC.index((1, 1))
    for magnitude, sign in ((0, 1), (2, 1), (2, -1)):
        embedded = lift(magnitude, sign, -1)
        # W times sqrt(2) is integral. Every column is orthogonal to every native write.
        np.testing.assert_array_equal(law._write_pattern().T @ embedded[law.N_SRC:], 0)
        np.testing.assert_array_equal(exchange @ embedded, ((1-sign)/2)*embedded)
        np.testing.assert_array_equal(contrast_cost @ embedded, magnitude**2*embedded)
        print(f"Coin-odd contrast {magnitude}, reversal parity {sign:+d}: "
              f"all nine load states acquire the same energy {KAPPA*magnitude**2 + (1-sign)/2:g}.")
    # A written even-coin control has the same empty-rail gap but feels a two-load contact.
    # Each signed probe is split between the marked port and an outside reference rail.
    # The rails have the same total load. Other terms and the reference evolution are isolated.
    for name, cp, expected in (("odd", -1, 1.), ("even", 1, 0.)):
        marked = lift(2, 1, cp).astype(float)
        marked /= np.linalg.norm(marked[:, 0])
        output = pulse @ marked
        phase = np.diag(marked.T @ output)
        np.testing.assert_allclose(output, marked * phase, atol=3e-13, rtol=0)
        for change in (-1, 1):
            here = law.SRC.index((1 + change, 1))
            bright = abs((phase[parked] + phase[here]) / 2)**2
            np.testing.assert_allclose(bright, expected, atol=3e-13, rtol=0)
            print(f"Coin-{name}, signed probe {change:+d}: bright probability {bright:.8f}.")


if __name__ == "__main__":
    run()
