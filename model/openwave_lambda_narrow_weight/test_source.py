from fractions import Fraction as F
import pytest
from .source import window_upper,clusters

def test_projection_bound_retains_an_exact_omitted_component():
    # H=diag(0,2), u=(1,1), x=(0,1/2), E=0. The full residual is (1,0).
    upper=window_upper(F(1),F(1,4),F(1,100),F(2))
    assert F(1,2)<=upper<F(51,100)
    assert window_upper(F(0),F(1,4),F(1,100),F(2))<F(1,2)
    with pytest.raises(ValueError):window_upper(F(-1),F(1),F(1),F(1))

def test_exact_source_resolvent_does_not_imply_a_density():
    # u at energy2, exact inverse x=u/2: any window of width<2 has zero weight.
    assert 0<window_upper(F(0),F(1,4),F(1,10),F(1))<1

def test_native_cluster_trials_keep_complete_ground_certificates():
    p,m,c=clusters()
    assert len(p)>1 and len(m)==18
    assert F(c['proton_prepared_ground_weight_lower_exact'])>F(99,100)
    assert F(c['pion_prepared_ground_weight_lower_exact'])>F(99,100)
