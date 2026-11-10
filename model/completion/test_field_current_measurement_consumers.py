"""Execute every unchanged consumer test against its named current custody service.

Only the packet/adapter globals of three private test namespaces are replaced.
Every test body, parametrization and fixture is loaded from authenticated
immutable source bytes; none is omitted or rewritten.
"""
import hashlib
import importlib.util
import inspect
import pytest
from pathlib import Path
from . import field_current_measurement as current

TEST_SOURCES={'local_energy_clock': ('model.native_local_energy_clock._field_current_tests', 'model/native_local_energy_clock/test_packet.py', 'b19a31769dbcc9472e4e23d37d6a615b16aef1e43af0876a051a996951585fc4'), 'fair_ground_readout': ('model.native_fair_ground_readout._field_current_tests', 'model/native_fair_ground_readout/test_packet.py', 'd038005012a782a7c2b7650891629d0ed376d799207cc287559b301ae9c56715'), 'adapter': ('model.completion._field_current_adapter_tests', 'model/completion/test_native_measurement_parent_adapter.py', 'f7c5bce2d37ea2a623de67827aa7e684d6e85826f58f7f41b029bef66a0b5051')}
EXPECTED={'local_energy_clock': ['test_changed_measurement_execution_witness_rejected', 'test_changed_measurement_owner_bound_rejected', 'test_changed_retained_receipt_rejected', 'test_current_complete_receipt_replays', 'test_duplicate_and_nonfinite_receipt_values_rejected', 'test_source_changes_during_replay_rejected', 'test_sparse_graph_counts_do_not_allocate_the_port_carrier'], 'fair_ground_readout': ['test_actual_words_and_complete_score_use_one_native_chart', 'test_bounded_local_arithmetic_does_not_allocate_large_regions', 'test_changed_prior_live_source_is_rejected', 'test_duplicate_and_nonfinite_receipts_fail', 'test_mid_replay_source_change_fails', 'test_prior_measurement_word_mutation_fails_exact_pin', 'test_retained_receipt_matches_current_complete_owner', 'test_retained_scientific_mutations_fail', 'test_source_inventory_is_local_and_has_no_research_note'], 'adapter': ['test_actual_parent_authenticates_retained_bytes_and_exact_live_science', 'test_archive_current_and_source_byte_guards', 'test_every_unreviewed_payload_or_source_transition_fails']}
for label,(name,path,digest) in TEST_SOURCES.items():
    full=current.ROOT/path
    assert hashlib.sha256(full.read_bytes()).hexdigest()==digest
    spec=importlib.util.spec_from_file_location(name,full)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    if label=='adapter':module.adapter=current.service()
    else:module.packet=current.owner('native_'+label)
    exported=[]
    for key,value in vars(module).items():
        if inspect.isfunction(value)and key.startswith('test_'):
            globals()['test_'+label+'_'+key.removeprefix('test_')]=value;exported.append(key)
        elif hasattr(value,'_pytestfixturefunction')or hasattr(value,'_fixture_function'):
            assert key not in globals();globals()[key]=value
    assert sorted(exported)==EXPECTED[label]


def test_exact_original_test_files_are_routed_to_this_complete_current_collector():
    config=(current.ROOT/'pytest.ini').read_text()
    for _,path,_ in TEST_SOURCES.values():assert '--ignore='+path in config
    assert len(TEST_SOURCES)==3


@pytest.fixture(autouse=True)
def current_replay_record_for_temporary_root_tamper_tests(tmp_path):
    # The original three mutation cases copy their original input scope.
    # Supply the one additional authenticated record required by the successor,
    # so each case reaches its intended archive/current/source guard.
    path = tmp_path / current.REPLAY
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((current.ROOT / current.REPLAY).read_bytes())
