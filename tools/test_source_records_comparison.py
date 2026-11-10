"""Only the retained c7 field cross-scale pointer is renewable."""
from copy import deepcopy
import pytest
from tools import renew_source_records_comparison as current
from tools import strict_json


def test_current_comparison_preserves_entire_science_and_scope():
    old = strict_json.load(current.ROOT / (current.ARCHIVE + current.RECEIPT))
    fresh = strict_json.load(current.ROOT / current.RECEIPT)
    assert current.conserve(old, fresh) == fresh


@pytest.mark.parametrize('damage', ['science', 'status', 'extra_source', 'lost_source', 'unknown_pointer'])
def test_comparison_rejects_science_and_unlisted_metadata(damage):
    old = strict_json.load(current.ROOT / (current.ARCHIVE + current.RECEIPT))
    fresh = deepcopy(strict_json.load(current.ROOT / current.RECEIPT))
    if damage == 'science':
        fresh['physical_prediction'] = 1
    elif damage == 'status':
        fresh['physical_gate_accepted'] = True
    elif damage in ('extra_source', 'lost_source'):
        _, mapping, _ = next(row for row in current.retained.source_maps(fresh) if isinstance(row[1], dict))
        if damage == 'extra_source':
            mapping['unreviewed/source.py'] = 'f' * 64
        else:
            mapping.pop(next(iter(mapping)))
    else:
        current.retained.set_path(fresh, current.POINTER, 'f' * 64)
    with pytest.raises(ValueError):
        current.conserve(old, fresh)


def test_old_authority_and_intermediate_are_distinct_authenticated_inputs():
    spec = current.policy()
    assert spec['retained_path'] == current.ARCHIVE + current.RECEIPT
    original = current.retained.load_catalog()['particle_comparison']
    assert original['retained_path'] != spec['retained_path']
    assert original['producer_source_sha256'] == spec['producer_source_sha256']
