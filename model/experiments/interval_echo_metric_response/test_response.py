"""Independent checks of the two-load reduction and the instrument used on it."""
import json

import numpy as np
from scipy import sparse
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_null_balance.experiment import Bond

from .measure import GAPS, HERE, SHOTS, hand_columns, instruments, propagate


def test_resident_meter_is_an_exact_subspace_of_full_two_load_board():
    board = federation(0)
    sector = PairSector(board)
    base, parked, gateway, _, _ = instruments(board, hand_columns())
    mobile = np.delete(np.arange(board.n), parked)
    retained = np.array([sector.index(parked, int(port)) for port in mobile])
    background = sector.hamiltonian(candidate.KAPPA) + 12*candidate.KAPPA*sparse.eye(sector.P)
    for neighbor in board.neighbours[parked]:
        background -= Bond(sector, (parked, neighbor)).h
    outside = np.ones(sector.P, dtype=bool)
    outside[retained] = False
    assert background[outside][:, retained].count_nonzero() == 0
    offset = candidate.KAPPA*(12-board.degree[parked])
    for gap in GAPS:
        source_energy = gap*(mobile != gateway)
        direct = background[retained][:, retained] + sparse.diags(source_energy)
        reduced = base[mobile][:, mobile] + sparse.diags(offset + source_energy)
        assert (direct-reduced).count_nonzero() == 0


def test_receiver_phase_against_independent_dense_unitaries():
    board = federation(0)
    base, _, gateway, near, _ = instruments(board, hand_columns())
    source = np.eye(board.n)[:, near]
    for gap in (0., *GAPS):
        dense = base.toarray()
        dense[gateway, gateway] -= gap
        reference = expm(-2j*dense) @ source
        observed = propagate(base, gateway, source, gap)[-1]
        np.testing.assert_allclose(observed, reference, rtol=0, atol=2e-12)


def test_local_records_and_resolved_controls():
    result = json.loads((HERE/"results.json").read_text())
    assert [board["cells"] for board in result["screens"]] == [20, 80, 320]
    for board in result["screens"]:
        for rows in board["arms"].values():
            assert all(np.all(np.asarray(r["counts"]).sum(axis=(1,2)) == SHOTS) for r in rows)
        for reading in board["differences"]:
            lo, hi = reading["restored_minus_empty"]["simultaneous_interval"]
            assert lo <= 0 <= hi
            lo, hi = reading["empty_minus_vacuum"]["simultaneous_interval"]
            assert lo <= 0 <= hi
            lo, hi = reading["resident_minus_empty"]["simultaneous_interval"]
            if reading["source"] in ("contact", "one_bridge"):
                assert hi < 0
            else:
                assert lo <= 0 <= hi
