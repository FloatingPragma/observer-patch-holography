"""Exact additive charge constraints for the registered ternary write and echo."""

import numpy as np
from sympy import Matrix

from model.candidate_laws.interval_echo.candidate import native_law as law


def write_constraints():
    """One integer conservation equation per sharp nonzero write amplitude.

    The unknowns are three common port charges followed by the ten memory
    charges. The blank memory supplies the zero of memory charge.
    """
    rows = []
    pattern = law._write_pattern()
    for target, source in zip(*pattern.nonzero(), strict=True):
        memory, output = divmod(int(target), law.N_SRC)
        row = [0] * (3 + law.N_LAB)
        for n in law.SRC[source]:
            row[n] += 1
        for n in law.SRC[output]:
            row[n] -= 1
        row[3 + memory] -= 1
        rows.append(row)
    return np.array(rows, dtype=int)


def compensating_memory(port):
    second_difference = port[0] + port[2] - 2 * port[1]
    return np.diag([second_difference if abs(d) == 2 else 0 for _, d in law.LABELS])


def dark_memory_projectors():
    """Integer multiples of three projectors orthogonal to every written word."""
    def coin_difference(contrast):
        return np.array([int(d == contrast) * (1 - 2 * coin) for coin, d in law.LABELS])

    left, middle, right = (coin_difference(d) for d in (-2, 0, 2))
    return [np.outer(v, v) for v in (middle, left + right, left - right)]


def commutes_with_seam(port, memory):
    """Check all 99-state echo and generator blocks with integer matrices."""
    pair = np.diag([port[a] + port[b] for a, b in law.SRC])
    written = np.kron(memory, np.eye(law.N_SRC, dtype=int)) + np.kron(
        np.eye(law.N_LAB, dtype=int), pair
    )
    pattern = law._write_pattern()
    lower_echo = 2 * np.eye(len(written), dtype=int) - pattern @ pattern.T
    reversal, contrast_squared = law._label_ops()
    return all((
        np.array_equal(written @ pattern, pattern @ pair),
        np.array_equal(written @ lower_echo, lower_echo @ written),
        np.array_equal(memory @ reversal, reversal @ memory),
        np.array_equal(memory @ contrast_squared, contrast_squared @ memory),
    ))


def exact_report():
    equations = write_constraints()
    neutral_rank = Matrix(equations[:, :3].tolist()).rank()
    charged_rank = Matrix(equations.tolist()).rank()
    signed = (-1, 0, 1)
    unsigned = (1, 0, 1)
    neutral_memory = np.zeros((law.N_LAB, law.N_LAB), dtype=int)
    results = {
        "nonzero_write_coefficients": len(equations),
        "port_dimensions_with_neutral_memory": 3 - neutral_rank,
        "port_and_label_dimensions_with_memory_charge": 3 + law.N_LAB - charged_rank,
        "unsigned_failures_with_neutral_memory": np.count_nonzero(equations[:, :3] @ unsigned),
        "signed_charge_survives": commutes_with_seam(signed, neutral_memory),
        "unsigned_charge_with_memory_survives": commutes_with_seam(unsigned, compensating_memory(unsigned)),
        "three_unwritten_memory_projectors_survive": all(
            commutes_with_seam((0, 0, 0), operator) for operator in dark_memory_projectors()
        ),
    }
    return {key: int(value) if not isinstance(value, bool) else value for key, value in results.items()}


if __name__ == "__main__":
    import json

    print(json.dumps(exact_report(), indent=2))
