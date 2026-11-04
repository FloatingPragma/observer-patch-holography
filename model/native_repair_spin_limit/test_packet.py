"""Reject scope promotion and stale proof custody in the retained receipt."""
from copy import deepcopy

import pytest

from . import algebra, packet


def test_retained_source_and_analytic_certificate_replay():
    value = packet.report()
    assert value['all_pass']
    assert value['bounds']['native_repair_rate_fixed_positive_on_every_edge']
    assert value['bounds']['volume_independent_local_interacting_limit']
    assert not value['bounds']['arbitrary_record_feedback_covered']
    assert not value['bounds']['interacting_relativistic_continuum']
    assert not value['bounds']['finite_c_faithful_inputs_stationary']


def test_relativistic_scope_promotion_is_rejected():
    changed = deepcopy(packet.report())
    changed['bounds']['interacting_relativistic_continuum'] = True
    with pytest.raises(ValueError, match='scope or custody'):
        packet.validate(changed)


def test_source_change_after_import_is_rejected_before_replay(monkeypatch):
    pins = packet.source_bindings()
    pins['model/native_repair_spin_limit/PROOF.md'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: pins)
    with pytest.raises(ValueError, match='source changed after import'):
        packet.report()


def test_mutated_repair_that_creates_a_charge_defect_is_rejected(monkeypatch):
    matrices = list(algebra.branch_matrices(2, 'C', (0, 1)))
    matrices[0] = matrices[0].copy()
    # A protected up/up input is deliberately sent to vacuum.
    matrices[0][:, 5] = 0
    matrices[0][0, 5] = 1
    monkeypatch.setattr(algebra, 'branch_matrices', lambda *_: tuple(matrices))
    with pytest.raises(ArithmeticError, match='defect number'):
        algebra.source_audit(2, 'C', (0, 1))
