"""Protect the historical/current distinction in the dated custody record."""
import copy
import pytest

from model.completion import custody_renewal as renewal
from tools import strict_json


def test_current_custody_preserves_scientific_status():
    result = renewal.validate()
    assert result['classification_counts'] == {'domains': 156, 'realized': 151}
    assert result['new_mathematical_acceptance'] is False


@pytest.mark.parametrize('damage, message', [
    ('count', 'classification counts differ'),
    ('history', 'historical review snapshot changed'),
    ('physical_acceptance', 'cannot promote scientific acceptance'),
    ('binding', 'current custody binding changed'),
    ('coverage', 'current custody binding coverage changed'),
])
def test_new_renewal_cannot_relabel_history_or_scientific_acceptance(damage, message):
    packet = copy.deepcopy(strict_json.load(renewal.ROOT / renewal.RECORD))
    if damage == 'count':
        packet['current_classification_counts']['domains'] = 145
    elif damage == 'history':
        packet['historical_snapshot_restorations'][0]['restored_sha256'] = '0' * 64
    elif damage == 'physical_acceptance':
        packet['new_mathematical_acceptance'] = True
    elif damage == 'coverage':
        del packet['current_bindings'][renewal.INTERFACES]
    else:
        first = next(iter(packet['current_bindings']))
        packet['current_bindings'][first] = '0' * 64
    with pytest.raises(ValueError, match=message):
        renewal.validate(packet)
