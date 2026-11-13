from fractions import Fraction as F
import numpy as np
from .free import symbols, corner_records, refinement_records


def test_shared_speed_and_actual_refinement_bounds():
    records = refinement_records([('electron', F(1)), ('muon', F(207)), ('tau',F(3477))])
    assert records['common_speed_squared_exact'] == '1/128'
    assert len(records['records']) == 36
    coarse = symbols((.1,.2,.3),1,1)
    fine = symbols((.1,.2,.3),F(1,8),1)
    assert abs(fine['positive_energy']-fine['continuum_Dirac_energy']) < abs(coarse['positive_energy']-coarse['continuum_Dirac_energy'])


def test_all_corner_branches_and_unregulated_doublers_retained():
    lifted = corner_records(F(1,4),F(1))
    control = corner_records(F(1,4),F(1),0)
    assert sum(F(r['lift_exact']) == 0 for r in lifted) == 1
    assert sum(F(r['lift_exact']) == 0 for r in control) == 8
    assert max(F(r['lift_exact']) for r in lifted) == F(3,2)


def test_error_bound_random_momenta_does_not_require_small_momentum_fit():
    rng = np.random.default_rng(802)
    for a in (1.,.25):
        for _ in range(20):
            v = symbols(rng.uniform(-2,2,3),a,1)
            assert abs(v['positive_energy']-v['continuum_Dirac_energy']) <= v['Dirac_absolute_error_bound']
            assert abs(v['photon_quadratic_energy']-v['continuum_photon_energy']) <= v['Maxwell_absolute_error_bound']
