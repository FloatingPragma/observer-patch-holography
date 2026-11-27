import numpy as np

from model.candidate_laws.interval_echo.candidate import native_law as law
from model.experiments.material_inertia.pair_spectrum.neutral_guide.charge_identity import (
    compensating_memory,
    dark_memory_projectors,
    exact_report,
    write_constraints,
)


def test_charge_spaces_and_full_echo():
    assert exact_report() == {
        "nonzero_write_coefficients": 18,
        "port_dimensions_with_neutral_memory": 2,
        "port_and_label_dimensions_with_memory_charge": 3,
        "unsigned_failures_with_neutral_memory": 4,
        "signed_charge_survives": True,
        "unsigned_charge_with_memory_survives": True,
        "three_unwritten_memory_projectors_survive": True,
    }


def test_each_port_basis_vector_and_changed_memory():
    equations = write_constraints()
    for port in np.eye(3, dtype=int):
        label = np.diag(compensating_memory(port))
        assert np.all(equations @ np.r_[port, label] == 0)
    for label in range(law.N_LAB):
        changed = np.zeros(3 + law.N_LAB, dtype=int)
        changed[3 + label] = 1
        assert np.any(equations @ changed)
    projectors = dark_memory_projectors()
    for i, first in enumerate(projectors):
        assert np.any(first)
        assert np.array_equal(first @ first, np.trace(first) * first)
        for second in projectors[i + 1:]:
            assert not np.any(first @ second)
