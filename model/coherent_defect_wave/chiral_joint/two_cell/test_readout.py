"""Complete original-attempt probabilities and charged selection rules."""
import sympy as s
from .readout import complete_experiment,unconditional_experiment,reflected_experiments,probe_matrix,sqrt_upper
from .spectrum import sector_certificate


def test_exact_probe_selection_excludes_wrong_charged_ground():
    target,_=probe_matrix(1,0,'remove')
    assert target==((0,1),3,-1)
    ground=sector_certificate(1)
    assert ground['ground_window_candidate_labels']==[{'L':4,'twoSz':0,'eigen_index':0}]
    assert all(row['L']!=target[1] for row in ground['ground_window_candidate_labels'])


def test_complete_readout_has_certifiably_visible_thresholds():
    result=complete_experiment()
    assert result['exact_complete_Kraus_metric_identity']
    assert not result['probabilities_renormalized']
    records={r['kind']:r for r in result['records']}
    assert sum(len(r['complete_energy_windows']) for r in records.values())==460
    for kind,lo,hi in [('add',s.Rational(105,100),s.Rational(106,100)),
                       ('remove',s.Rational(90,100),s.Rational(92,100))]:
        first=records[kind]['complete_energy_windows'][0]
        elo,ehi=map(s.Rational,first['energy_transfer_interval'])
        assert lo<elo<ehi<hi
        assert s.Rational(first['probability_interval'][0])>s.Rational(1,5)
    assert set(records)=={'add','remove','add_fail','remove_fail'}


def test_ground_preparation_failure_stays_in_original_attempt_law():
    result=unconditional_experiment()
    assert result['record_count_including_preparation_failure']==461
    failure=result['records'][0]
    assert failure['probability_interval']==['1052/1053']*2
    low,high=map(s.Rational,result['total_probability_enclosure'])
    assert low<=1<=high
    success_upper=sum(s.Rational(r['probability_interval'][1]) for r in result['records'][1:])
    assert success_upper<s.Rational(1,1000)
    # Dropping the preparation failure cannot masquerade as a normalized law.
    assert success_upper<1


def test_reflected_charge_pairs_use_the_actual_symmetry():
    result=reflected_experiments()
    assert not result['same_endpoint_addition_removal_equality_asserted']
    assert len(result['pairs'])==4
    assert sum(r['compared_energy_windows'] for r in result['pairs'])==460


def test_square_root_bound_is_outward():
    for q in (s.Rational(0),s.Rational(2),s.Rational(1,7)):
        upper=sqrt_upper(q)
        assert upper**2>=q
        if q:assert (upper-s.Rational(1,1<<80))**2<q
