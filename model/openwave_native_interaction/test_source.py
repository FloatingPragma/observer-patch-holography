from fractions import Fraction as F
import pytest
from .source import run,unread,trace,assay_zero_imbalance,event_maps

@pytest.fixture(scope='module')
def receipt():return run()

def test_full_operator_and_original_compiler(receipt):
 assert receipt['complete_operator_control']['matrix_units_checked']==128
 c=receipt['original_compiler_control'];assert c['status']=='COMPLETE'and c['raw_histories']==625
 assert c['all_H66_tagged_column_history_checks']==41250
 assert receipt['native_word']['full_carrier_dimension']==28680

def test_interaction_survives_identical_one_point_means(receipt):
 rows=receipt['two_recording_block_results']
 assert [r['pair_density_difference_exact']for r in rows]==['1/4','9/64']
 assert [r['operational_zero_imbalance_difference_exact']for r in rows]==['1/4','3/64']
 for r in rows:assert r['native_one_point_means_exact']==r['independent_one_load_control']['one_point_means_exact']
 assert receipt['actual_record_current_moments_exact']['connected_current_correlation']=='1/8'

def test_assay_keeps_failure_and_backaction_domain():
 # Full8-coordinate operator: no success normalization and no arbitrary
 # assertion that the terminal instrument preserves coherence.
 rho={(3,3):F(1,2),(6,6):F(1,2),(3,6):F(1,2),(6,3):F(1,2)}
 value,out=assay_zero_imbalance(rho)
 assert value==F(1,2)and trace(out)==1 and out!=rho
 assert len(event_maps((1,2)))==25

def test_unsupported_unitary_and_long_history_scope(receipt):
 assert receipt['unsupported_instruction_witness']['status']=='REJECTED'
 assert receipt['unsupported_instruction_witness']['source_operations_executed']==0
 assert not receipt['scope']['unitary_interacting_matter_field_compiled']
 assert not receipt['native_word']['expanded_five_event_history_tree_executed']
 assert receipt['native_word']['complete_raw_words']==25**5
