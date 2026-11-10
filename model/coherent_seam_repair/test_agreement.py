"""Independent finite-premise checks, including the wrong-apparatus control."""
import pytest
import sympy as s
from .agreement import ChartedObserver, certificate
from .program import Observer


def test_complete_chart_agreement_and_wrong_detector():
    result = certificate()
    assert result['oriented_local_operator_identities'] == 60
    assert len(result['complete_adaptive_outcomes']) == 4
    assert sum(s.Rational(row['probability']) for row in result['complete_adaptive_outcomes']) == 1
    assert {row['preparation'] for row in result['complete_adaptive_outcomes']} == {'success','failure'}
    assert s.Rational(result['wrong_untransported_detector_probability']) == s.Rational(1,120)
    assert s.Rational(result['correct_transported_detector_probability']) == s.Rational(1,60)


def test_unowned_frame_and_wrong_source_are_rejected():
    with pytest.raises(ValueError):
        ChartedObserver(Observer(),{})
    with pytest.raises(ValueError):
        ChartedObserver(Observer(fibre=5),{12:[]})
    with pytest.raises(ValueError):
        ChartedObserver(Observer(fibre=5),{0:[9]})
