"""Actual interpreter agreement, causal ordering and retained stopping tests."""
from fractions import Fraction as F
from itertools import product
import pytest
from . import source as s


def test_complete_endpoint_instrument_matches_original_every_matrix_unit():
    g = s.graph(0)
    seam = g.inter_seams[0]
    masks = [((x & 1) << seam[0]) | ((x >> 1) << seam[1]) for x in range(4)]
    # N-changing matrix units are excluded from a fixedN admission. Every
    # same-N pair is tested through its actual declared sector.
    for law in ('C', 'E', 'mixture'):
        for left, right in product(masks, repeat=2):
            if left.bit_count() != right.bit_count():
                continue
            p = s.admitted_program(0, left.bit_count(), [s.fixed_command(seam, law)])
            output = s.native.execute(p, {(left, 0, right, 1): F(1)})
            expected = s.branch_column(left, right, seam, law)
            assert len(output.branches) == len(expected) == 25
            for actual, (_, weight, image) in zip(output.branches, expected, strict=True):
                wanted = {} if not weight else {(image[0], 0, image[1], 1): weight}
                assert actual.density == wanted


def test_partial_order_rejects_overlapping_reversal():
    g = s.graph(0); a, b = g.inter_seams[0]
    c = next(v for v in g.neighbours[a] if v != b)
    seams = (tuple(sorted((a, c))), (a, b))
    with pytest.raises(ValueError, match='causal'):
        s.execute_event_tape(0, 1, seams, {(1 << a, 0, 1 << a, 0): F(1)}, (1, 0))


def test_endpoint_identity_and_disjoint_raw_instruments():
    r = s.endpoint_certificate()
    assert r['disjoint_outcome_matrix_unit_comparisons'] == 160000
    assert any(not row.pairs for row in s.raw_local_branches())


def test_stopping_and_changed_rule_are_not_postselected():
    r = s.local_shell_and_stop_controls()
    assert r['FIFO_timeout']['status'] == 'FIFO_TIMEOUT'
    assert r['host_resource_stop']['status'] == 'RESOURCE_STOP'
    assert r['altered_rule_control']['altered_rule_arrival_probability_exact'] == '1/4'


def test_forbidden_target_instruction_rejected_by_original_admission():
    with pytest.raises(ValueError):
        s.admitted_program(commands=[{'op': 'target_Hamiltonian', 'repeat': 1}])
