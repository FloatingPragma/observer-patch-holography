from copy import deepcopy

import pytest
from tools import strict_json

from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_complete_cold_replay(replay):
    packet.validate(strict_json.load(packet.OUTPUT),replay)


@pytest.mark.parametrize('field',['interaction_strictly_site_local','native_CE_prepares_Gibbs_reference',
    'invariant_four_mode_slice_is_full_core_vacuum','unchanged_MCPR_coefficients_are_interacting_physical_masses',
    'G1_through_G6_joint_geometric_premises_realized','Einstein_geometry_derived','integrated_gate_accepted'])
def test_false_promotions_rejected(replay,field):
    wrong = deepcopy(replay)
    wrong['scope'][field] = True
    with pytest.raises(ValueError,match='complete replay'):
        packet.validate(wrong,replay)


def test_missing_states_and_outcomes_rejected(replay):
    wrong = deepcopy(replay)
    wrong['complete_record_nonlinear_control']['retained_fine_records_per_input'] -= 1
    with pytest.raises(ValueError,match='complete replay'):
        packet.validate(wrong,replay)


def test_no_document_or_global_runner_bound_as_evidence():
    paths = packet.source_paths()
    assert not any(path.startswith('docs/') or path == 'model/verify.py' for path in paths)
