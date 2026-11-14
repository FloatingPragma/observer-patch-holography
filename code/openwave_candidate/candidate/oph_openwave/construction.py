"""Finite first-order construction from REPORT_FOR_CODEX R4/R9.

This implements algebra and literal committed events, not a quantum/classical
hybrid channel or an OpenWave validation. Integer records are current counters;
seam amplitudes in H_m are transit amplitudes, not those counters.
"""
from dataclasses import dataclass
import operator
import numpy as np


@dataclass(frozen=True)
class Construction:
    ports: int
    seams: tuple[tuple[int, int], ...]
    mass: float = 0.0

    def __post_init__(self):
        if isinstance(self.ports, (bool, np.bool_)) or not isinstance(self.ports, (int, np.integer)):
            raise ValueError("integer port count required")
        if isinstance(self.mass, (bool, np.bool_)) or not isinstance(self.mass, (int, float, np.integer, np.floating)):
            raise ValueError("real rest parameter required")
        if self.ports < 1 or not np.isfinite(self.mass) or self.mass < 0:
            raise ValueError("positive port count and finite nonnegative mass required")
        for a, b in self.seams:
            if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer)) for v in (a, b)):
                raise ValueError("integer seam endpoints required")
            if not 0 <= a < self.ports or not 0 <= b < self.ports or a == b:
                raise ValueError("invalid oriented seam")

    @property
    def incidence(self):
        b = np.zeros((self.ports, len(self.seams)))
        for e, (a, c) in enumerate(self.seams):
            b[a, e], b[c, e] = 1, -1
        return b

    @property
    def hamiltonian(self):
        b = self.incidence
        return np.block([[self.mass * np.eye(self.ports), b],
                         [b.T, -self.mass * np.eye(len(self.seams))]])

    def evolve(self, state, time):
        """Unitary evolution in native hop units, with no global seam clock."""
        values, vectors = np.linalg.eigh(self.hamiltonian)
        state = np.asarray(state, dtype=complex)
        if state.shape != (len(values),) or not np.all(np.isfinite(state)) or not np.isfinite(time):
            raise ValueError("state dimension or time invalid")
        return vectors @ (np.exp(-1j * time * values) * (vectors.T @ state))

    def evolve_independent(self, states, time):
        """Product preparation under sum_i H_i (one transit space per load).

        Returning factors avoids allocating an exponentially large tensor. This
        API makes no claim about interacting loads or exchange statistics.
        """
        return tuple(self.evolve(state, time) for state in states)

    def committed_event(self, loads, records, seam, coin):
        """Literal integer event; return new registers and signed current."""
        n = [operator.index(v) for v in loads]
        e = [operator.index(v) for v in records]
        if len(n) != self.ports or len(e) != len(self.seams):
            raise ValueError("register dimensions differ from construction")
        if coin not in (0, 1) or not 0 <= seam < len(self.seams):
            raise ValueError("invalid event")
        a, b = self.seams[seam]
        delta = n[a] - n[b]
        current = (delta - (2 * coin - 1) * (abs(delta) % 2)) // 2
        n[a] -= current
        n[b] += current
        e[seam] += current
        return tuple(n), tuple(e), current

    def gauss_register(self, loads, records):
        """Exact n + div E, using Python integers even beyond int64."""
        n = [operator.index(v) for v in loads]
        e = [operator.index(v) for v in records]
        if len(n) != self.ports or len(e) != len(self.seams):
            raise ValueError("register dimensions differ from construction")
        for (a, b), value in zip(self.seams, e):
            n[a] += value
            n[b] -= value
        return tuple(n)


def run():
    c = Construction(3, ((0, 1), (1, 2), (2, 0)), 0.19)
    b, h = c.incidence, c.hamiltonian
    expected = np.block([[b @ b.T + c.mass**2*np.eye(3), np.zeros((3,3))],
                         [np.zeros((3,3)), b.T @ b + c.mass**2*np.eye(3)]])
    state = np.eye(6)[0]
    first, _ = c.evolve_independent((state, np.eye(6)[2]), .71)
    n, e = (13, -4, 7), (0, 0, 0)
    initial = c.gauss_register(n, e)
    for step in range(1000):
        n, e, _ = c.committed_event(n, e, step % 3, step % 2)
    probe = np.array([3, -1, 7])
    cost = sum((probe[a]-probe[d])**2 for a,d in c.seams)
    return {"scope": "Finite first-order incidence construction and literal committed current counters; no hybrid-channel or physical continuum claim.",
            "checks": {"massive_square": bool(np.allclose(h@h, expected)),
                       "independent_load_factor": bool(np.allclose(first, c.evolve(state, .71))),
                       "gauss_integer_identity": c.gauss_register(n,e) == initial,
                       "native_cost": bool(probe@b@b.T@probe == cost)},
            "mutations": {"missing_current_detected": c.gauss_register(n, (e[0]+1,*e[1:])) != initial,
                          "wrong_mass_square_detected": not bool(np.allclose(h@h, expected + np.eye(6)))},
            "full_criterion_met": False,
            "supplied": {"mass": c.mass, "graph": "oriented triangle", "time_unit": "inverse hop amplitude"}}


if __name__ == "__main__":
    import json
    print(json.dumps(run(), indent=2))
