"""Only the declared parent/campaign hashes may change in either diagnostic."""
from copy import deepcopy
import pytest
from tools import renew_diagnostic_parent_custody as renewal
from tools import strict_json


@pytest.mark.parametrize('owner', renewal.OWNERS)
def test_current_diagnostic_has_exact_retained_science_and_source_paths(owner):
    authority = renewal.authority()
    spec = authority['owners'][owner]
    old = strict_json.load(renewal.ROOT / spec['original_receipt'])
    current = strict_json.load(renewal.ROOT / spec['receipt'])
    assert renewal.conserve_receipt(spec, old, current, authority['parent']['current_sha256'],
                                    renewal.sha(spec['campaign'])) == current


@pytest.mark.parametrize('owner', renewal.OWNERS)
@pytest.mark.parametrize('damage', ['status', 'status_type', 'predecessor', 'source_add', 'source_delete', 'unlisted_hash'])
def test_unlisted_campaign_mutation_is_rejected(owner, damage):
    authority = renewal.authority()
    spec = authority['owners'][owner]
    old = strict_json.load(renewal.ROOT / spec['original_campaign'])
    current = renewal.projected_campaign(old, authority['parent']['current_sha256'])
    if damage == 'status':
        current['integrated_gate_accepted'] = True
    elif damage == 'status_type':
        current['integrated_gate_accepted'] = 0
    elif damage == 'predecessor':
        current['predecessor_git_commit'] = 'f' * 40
    elif damage == 'source_add':
        current['source_sha256']['unreviewed/source.py'] = 'f' * 64
    elif damage == 'source_delete':
        current['source_sha256'].pop(next(iter(current['source_sha256'])))
    else:
        target = next(p for p in current['source_sha256'] if p != renewal.PARENT)
        current['source_sha256'][target] = 'f' * 64
    with pytest.raises(ValueError, match='campaign changed outside'):
        renewal.conserve_campaign(old, current, authority['parent']['current_sha256'], allow_original=True)


@pytest.mark.parametrize('owner', renewal.OWNERS)
@pytest.mark.parametrize('damage', ['science', 'status', 'status_type', 'contract', 'source_add', 'source_delete'])
def test_changed_science_or_source_inventory_cannot_be_repaired_by_new_hashes(owner, damage):
    authority = renewal.authority()
    spec = authority['owners'][owner]
    old = strict_json.load(renewal.ROOT / spec['original_receipt'])
    current = deepcopy(strict_json.load(renewal.ROOT / spec['receipt']))
    if damage == 'science':
        current['new_mass_claim'] = '1 GeV'
    elif damage == 'status':
        current['native_family_changed'] = True
    elif damage == 'status_type':
        current['native_family_changed'] = 0
    elif damage == 'contract':
        current['contract_sha256'] = 'f' * 64
    elif damage == 'source_add':
        current['source_sha256']['unreviewed/source.py'] = 'f' * 64
    else:
        current['source_sha256'].pop(next(iter(current['source_sha256'])))
    with pytest.raises(ValueError):
        renewal.conserve_receipt(spec, old, current, authority['parent']['current_sha256'],
                                 renewal.sha(spec['campaign']))
