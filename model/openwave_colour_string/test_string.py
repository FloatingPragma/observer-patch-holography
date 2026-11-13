from fractions import Fraction as F
import numpy as np
from .source import chain,native_transfer
from .certificate import execute,uniform_confinement,unit_elimination

def test_exact_native_transfer_and_coloured_countercontrol():
    d=native_transfer()
    assert d['all_native_baryon_singlets']==120 and d['all_nuclear_magnetic_singlets']==2520
    assert d['nuclear_local_colour_difference_times_exact_complete_projector_nonzeros']==0
    assert d['coloured_complement_identity_fails_and_is_retained']

def test_fixed_source_chain_and_complete_band_records():
    two=execute(2);four=execute(4)
    assert two['parameters_fixed_without_graph_refitting']==four['parameters_fixed_without_graph_refitting']
    assert abs(two['native_gA_diagnostic']-1.2753)<1e-8 and four['native_gA_diagnostic']>1.66
    assert four['lowest_cluster_rank']==2
    assert len(four['all_colour_singlet_eigenvalues_diagnostic'])==816
    r=four['complete_native_records']
    assert abs(sum(r['all_native_occupation_success_probabilities_diagnostic'])-float(F(r['preparation_success_exact'])))<1e-11
    assert abs(sum(r['all_native_occupation_failure_probabilities_diagnostic'])-float(F(r['preparation_failure_exact'])))<1e-11

def test_string_law_uniform_bound_and_shared_length_obstruction():
    d=chain(4);prefix=np.cumsum(d['n'],axis=1)[:,:-1]
    assert np.array_equal(2*np.sum(prefix*(3-prefix),axis=1),4*d['spans'])
    assert F(uniform_confinement()['exact_uniform_ground_mean_span_upper'])>0
    c=unit_elimination()
    assert not c['simultaneous_benchmark_match_possible_with_nonnegative_adjusted_U']
    assert F(c['any_positivity_preserving_kappa_string_tension_squared_upper'])<1
    assert F(c['required_kappa_squared_lower'])>F(c['positivity_preserving_kappa_upper'])**2
