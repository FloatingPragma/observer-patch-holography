"""Every actual sign survives the registered blank reference and its restrictions."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import math
import random

import pytest

from . import registered_reference as ref


class FixedTicket:
    def __init__(self, ticket, total):
        self.ticket, self.total = ticket, total

    def randrange(self, total):
        assert total == self.total
        return self.ticket


def test_all_raw_records_and_all_sampling_tickets():
    for raw in product((-1, 1, None), repeat=2):
        kernel = ref.completion_kernel(raw)
        total = math.lcm(*(p.denominator for p in kernel))
        outcomes = Counter(ref.complete(raw, FixedTicket(k, total)) for k in range(total))
        assert tuple(F(outcomes[pair], total) for pair in ref.PAIRS) == kernel
        for pair in outcomes:
            assert all(physical is None or physical == sign
                       for physical, sign in zip(raw, pair, strict=True))


def test_joint_source_stocks_exact_restriction_composition_and_positive_divergence():
    result = ref.verify()
    assert result["basis_and_inclusion_checks"] == 36
    ordinary, held = (result["stocks"][s] for s in ("ordinary", "held"))
    assert 0 < ordinary["entropy_change"] < held["entropy_change"]
    assert 0 < ordinary["divergence"] < held["divergence"]
    with pytest.raises(ValueError, match="never observed"):
        ref.restrict_record((-1, 1), 1, 2, random.Random(9))
    assert ref.restrict_record((-1, 1), 3, 3, random.Random(16)) == (-1, 1)


def test_unequal_physical_masks_share_one_registered_blank_and_all_regional_entropy():
    held = (ref.joint(F(1, 4)), ref.joint(F(3, 8)))
    held_masks, flat_masks = (3, 1), (2, 0)
    def tensor(a, b):
        return tuple(i*j for i, j in product(a, b))
    p = tensor(*(ref.projected_law(law, mask) for law, mask in zip(held, held_masks)))
    q = tensor(*(ref.projected_law(ref.BLANK, mask) for mask in flat_masks))
    assert q == tensor(ref.BLANK, ref.BLANK)
    assert math.isclose(ref.entropy(p)-ref.entropy(q),
                        ref.entropy(held[0])-ref.entropy(ref.BLANK), abs_tol=1e-15)
    assert math.isclose(ref.divergence(p, q), ref.divergence(held[0], ref.BLANK), abs_tol=1e-15)
    # Without conditioning or padding the two physical raw blank backgrounds differ.
    assert (ref.entropy(ref.BLANK)+math.log(2))-math.log(2) > 1
