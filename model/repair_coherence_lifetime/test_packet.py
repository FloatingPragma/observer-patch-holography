"""Retained evidence must reject silent scope and JSON type changes."""
from copy import deepcopy

import pytest

from .packet import compare


@pytest.mark.parametrize('field,value', [
    ('all_pass', 1), ('physical_comparison_count_added', False),
    ('integrated_gate_accepted', True), ('coherence_factor', '1'),
])
def test_changed_value_or_json_type_rejected(field, value):
    original = {'all_pass': True, 'physical_comparison_count_added': 0,
                'integrated_gate_accepted': False, 'coherence_factor': '5/6'}
    changed = deepcopy(original)
    changed[field] = value
    with pytest.raises(ArithmeticError):
        compare(changed, original)


def test_nonfinite_evidence_rejected():
    with pytest.raises(ValueError):
        compare({'gap': float('nan')}, {'gap': float('nan')})
