from dataclasses import replace
from fractions import Fraction as F
import pytest
from .field import Q
from .source import graph_vectors,population_and_record_check,run
from model.openwave_native_process import source as process

@pytest.fixture(scope='module')
def receipt():return run()

def test_exact_quadratic_field_and_three_graph_modes():
 assert Q(0,1)**2==Q(5)
 assert Q(3,1)/Q(3,1)==Q(1)
 _,rows=graph_vectors();a,b,c=[r[2]for r in rows]
 assert a+b==Q(F(11,6))and a*b==Q(F(151,180))and c==Q(F(9,10))

def test_complete_current_modes_and_host_boundary(receipt):
 assert receipt['complete_basis_current_control']['basis_columns']==66
 assert receipt['event_clock']['complete_raw_record_words']==721**4
 assert receipt['event_clock']['original_prefix_status']=='RESOURCE_STOP'
 assert not receipt['event_clock']['expanded_full_history_tree_executed']
 for row in receipt['mode_records']:
  assert row['event_multiplier_exact']==row['recovered_event_multiplier_exact']
  assert all(sum(map(F,p))==1 for p in row['four_complete_population_distributions_exact'])
 assert not receipt['scope']['physical_mass_phase_or_GeV_assignment']

def test_altered_native_current_is_rejected(monkeypatch):
 rows=list(process.raw_local_branches());i=next(i for i,r in enumerate(rows)if r.weight and r.pairs)
 rows[i]=replace(rows[i],current=rows[i].current+1)
 monkeypatch.setattr(process,'raw_local_branches',lambda law='mixture':tuple(rows))
 with pytest.raises(AssertionError):population_and_record_check(graph_vectors()[1])

def test_sampling_scope_and_actual_record_means(receipt):
 assert F(receipt['sampling_plan']['simultaneous_failure_upper_exact'])<F(1,100)
 for row in receipt['mode_records']:
  means=[Q(*p)for p in row['four_actual_record_means_exact']]
  assert all(means[i+1]/means[i]==Q(*row['event_multiplier_exact'])for i in range(3))
