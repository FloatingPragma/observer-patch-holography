"""Canonical renewal admits one documentary pin and no scientific changes."""
from copy import deepcopy
import pytest
from tools import renew_native_control_canonical as renewal


@pytest.mark.parametrize('damage', ['science', 'type', 'source_scope', 'executable'])
def test_pilot_rejects_every_change_outside_documentary_custody(damage):
    packet = deepcopy(renewal.retained(renewal.CONTROL))
    if damage == 'science':
        packet['scope']['physical_energy_or_mass_derived'] = True
    elif damage == 'type':
        packet['carrier']['states'] = float(packet['carrier']['states'])
    elif damage == 'source_scope':
        packet['source_sha256']['unreviewed.py'] = 'f' * 64
    else:
        packet['source_sha256']['model/common_source/native_defect_decision.py'] = 'f' * 64
    with pytest.raises(ValueError):
        renewal.compare_source_receipt(packet, renewal.CONTROL, {renewal.DOCUMENT})


def test_archived_inputs_have_their_exact_pre_renewal_bytes():
    for name in renewal.EXPECTED:
        renewal.retained(name)
