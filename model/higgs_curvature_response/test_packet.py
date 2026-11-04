"""Receipt scope, exact coefficient and cold-replay controls."""
from copy import deepcopy
import subprocess
import sys

import pytest

from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_fresh_report_and_scoped_boundary(fresh):
    packet.validate(fresh)
    assert fresh['canonical_physical_premises_discharged'] == []
    assert not any(name.startswith('tracking/') for name in fresh['source_bindings'])
    assert fresh['regional_first_variation']['first_variation_only']


@pytest.mark.parametrize('damage', ['volume', 'wald', 'force', 'nonlinear', 'source', 'type'])
def test_receipt_mutations_fail_closed(fresh, monkeypatch, damage):
    monkeypatch.setattr(packet, 'report', lambda: fresh)
    value = deepcopy(fresh)
    if damage == 'volume':
        value['wald_and_generalized_volume']['generalized_volume_density'] = '1'
    elif damage == 'wald':
        value['wald_and_generalized_volume']['Wald_area_coefficient'] = '4*pi*F0'
    elif damage == 'force':
        value['curved_force_balance']['nonzero_force_terms']['Yukawa_exchange'] = ['0']*4
    elif damage == 'nonlinear':
        value['scope']['nonlinear_Einstein_equation_derived'] = True
    elif damage == 'source':
        value['source_bindings'].pop(next(iter(value['source_bindings'])))
    else:
        value['scope']['source_regional_entropy_or_stationarity_constructed'] = 0
    with pytest.raises(ValueError, match='fresh calculation'):
        packet.validate(value)


def test_independent_cold_replay(tmp_path, fresh):
    path = tmp_path/'higgs_curvature.json'
    path.write_text(packet.canonical(fresh))
    result = subprocess.run([sys.executable, '-W', 'error', '-m',
                             'model.higgs_curvature_response', '--check', str(path)],
                            cwd=packet.ROOT, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout+result.stderr
    assert 'replay passed' in result.stdout
