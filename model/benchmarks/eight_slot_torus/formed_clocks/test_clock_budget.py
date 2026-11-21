import math
from fractions import Fraction

import numpy as np
import pytest

from . import clock_budget as budget
from .clock_inference import _hoeffding_lower, decode_elapsed


def test_fixed_budget_and_best_possible_known_phase():
    design = budget.design()
    assert design["difference_arc_width_upper"] < 2*math.pi/4096
    assert design["axis_count_outage_upper"] < budget.ALPHA
    assert budget.necessary_equal_shots([1.]) == 21128268
    assert budget.necessary_equal_shots([0.]) is None
    assert budget.necessary_equal_shots([.5,.5]) > budget.necessary_equal_shots([1.])
    with pytest.raises(ValueError):
        budget.necessary_equal_shots([1.01])


def test_outward_binomial_limits_and_artificial_clock():
    n = 9
    tail = Fraction(1,100)
    for k in range(1,n+1):
        p = Fraction(_hoeffding_lower(k,n,tail))
        assert sum(math.comb(n,j)*p**j*(1-p)**(n-j) for j in range(k,n+1)) <= tail
    def word(phase):
        p = np.array([1+math.cos(phase),1+math.sin(phase),1-math.cos(phase),1-math.sin(phase)])
        counts = np.rint(.07*budget.SHOTS*p/4).astype(np.int64)
        return np.array([[*counts,budget.SHOTS-counts.sum()]])
    first, second = word(-.4), word(-.4-198*math.pi/4096)
    options = dict(rate=-math.pi/4096,interval="hoeffding",alpha=budget.ALPHA,
                   family_size=budget.FAMILY,cell_error=budget.ERROR)
    assert decode_elapsed(first,second,**options)["even_runs_step_two"] == [[198,198]]
    assert decode_elapsed(second,first,**options)["even_runs_step_two"] == [[7994,7994]]


def test_authenticated_reduced_comparison():
    result = budget.experiment()
    assert result["counted_launches"] == 54*2**34
    assert result["early_active_launches"] == result["each_later_active_launches"] == 2*2**34
    assert not result["failures"]
    primary = result["primary_resource_audit"]
    assert primary["all_menu_launches"] == 16*result["counted_launches"]
    assert primary["necessary_shots_unit_incidence"] == 21128268
    assert all(r["decoded_even_runs"] == [[0,8190]] for r in result["records"] if r["menu"]=="dephased")
