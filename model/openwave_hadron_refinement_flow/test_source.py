from fractions import Fraction as F
import sympy as sy
from .certificate import calculate
from .certified_rest import inertia,eigen_intervals
from .positivity import joint_stability
from .exact_native import certificate


def test_exact_native_car_reduction_retains_high_outputs():
    d=certificate()
    assert d['exact_frame_rank']==16
    assert d['complete_fixed_projection_native_dimension']==104
    assert d['exact_coupling_square_identity']


def test_inertia_counts_and_calibrated_full_source_intervals():
    A=sy.Matrix([[2,1],[1,2]]);G=sy.eye(2)
    intervals=eigen_intervals(A,G)
    assert intervals[0][0]<=1<=intervals[0][1]
    assert intervals[1][0]<=3<=intervals[1][1]
    assert inertia(A-2*G)==1
    d=calculate();r=d['selected_bulk_neutron']
    lo,hi=map(F,r['full_source_GeV_interval'])
    assert lo<F('.93956542194')<hi and hi-lo<F('0.000000071')
    assert F(d['fixed_coefficient_bulk_control']['full_source_axial_interval'][0])>F('1.29')
    assert all(F(x['filter_success_on_uniform_physical_input'])+F(x['filter_failure_on_uniform_physical_input'])==1 for x in d['complete_preparation_records'])


def test_full_capacity_count_bound_keeps_nuclear_gate():
    d=joint_stability(F(13))
    assert d['uniform_floor']>F('3.5')
    assert d['cases']['three']>100
