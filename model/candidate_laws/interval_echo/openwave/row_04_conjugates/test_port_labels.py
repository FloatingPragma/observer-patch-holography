import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo.openwave.row_04_conjugates.experiment import (
    LABELS, contrast_diagonal, verify_reader,
)
from model.executable_model import cell, federation


def test_reference_terms_break_port_relabeling_only_where_degrees_differ():
    for car in (cell(), federation(0)):
        sector = PairSector(car)
        original = contrast_diagonal(sector, 0, 1)
        equivalent = {p for p in LABELS if np.ptp(contrast_diagonal(sector, *p) - original) == 0}
        expected = {(0, 1)} if car.cells == 20 else {(0, 1), (1, 0), (1, 2), (2, 1)}
        assert equivalent == expected


def test_zero_written_contrast_resolves_equal_cores_in_every_binary_subsector():
    for background, core in LABELS:
        verify_reader(background, core)
