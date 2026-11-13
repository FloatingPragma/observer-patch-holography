from fractions import Fraction as F
import pytest
from .memory import correction,probe_defects,run


def test_full_complement_computed_correction_and_memory_remainder():
    result=run();matrix,rows=correction()
    assert result['actual_computed_self_energy']['positive_nonzero_correction_executed']
    assert result['actual_computed_self_energy']['infinite_complement_Neumann_ratio_upper_exact']=='9/89'
    assert matrix[0,0]==0 and matrix[1,1]>0 and len(rows)>201
    assert result['retained_vacuum_second_step_leakage_norm_squared_exact']=='303/512'


def test_same_clock_projection_improvement_does_not_assert_autonomy():
    result=probe_defects(F(1,1024))
    assert F(result['retained_upper_divided_by_boundary_lower_exact'])<F(1,225)
    assert not result['interacting_vacuum_or_operator_norm_improvement_claimed']
    assert probe_defects(0)['boundary_only_vacuum_semigroup_defect_interval_exact']==['0','0']
    with pytest.raises(ValueError):probe_defects(-1)
