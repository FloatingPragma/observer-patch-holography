"""Independent matrix functions, finite cutoff identities and obstruction controls."""
from fractions import Fraction as F
import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm,sqrtm
import sympy as s
from model.openwave_native_force.force import source_matrices
from model.native_equilibration.vacuum_repair import hamiltonian
from . import geometry,clocks,waves


def test_native_clock_matrices_and_all_ratio_obstructions():
    _,radii,_=geometry.prepared();r,m,_,_=source_matrices()
    for index in (15,23):
        actual=hamiltonian(2,2,radii[index])/np.sqrt(14)
        np.testing.assert_allclose(actual,np.asarray(m.subs(r,s.Rational(radii[index])),float),atol=2e-15,rtol=0)
    result=clocks.certificate()
    assert len(result['rows'])==7
    for row in result['rows']:
        low,high=map(F,row['ratio_minus_metric_lapse_interval'])
        assert low>0 or high<0


def test_interval_sine_against_high_precision_at_all_phase_scales():
    from model.finite_area_response.response import Interval
    with mp.workdps(85):
        for phase in (F(0),F(1,10**10),F(-1,3),F(2),F(5),F(8)):
            actual=mp.sin(mp.mpf(phase.numerator)/phase.denominator)
            bound=clocks.sine(Interval.point(phase))
            assert bound.hi-bound.lo<F(1,10**20)
            low=mp.mpf(bound.lo.numerator)/bound.lo.denominator
            high=mp.mpf(bound.hi.numerator)/bound.hi.denominator
            assert low<=actual<=high
    with pytest.raises(ValueError):clocks.sine(Interval.point(9))


def test_Ramsey_Y_sign_from_direct_two_level_interferometer():
    result=clocks.certificate()
    plus_y=np.array([1,1j])/np.sqrt(2)
    for row in result['rows']:
        for gaps,probability in zip(row['gap_over_sqrt14_intervals'],row['Ramsey_plus_Y_probability_intervals']):
            phase=float(sum(map(F,gaps))/2)
            state=np.array([1,np.exp(-1j*phase)])/np.sqrt(2)
            observed=abs(np.vdot(plus_y,state))**2
            lo,hi=map(F,probability)
            assert abs(observed-float((lo+hi)/2))<1e-15


def test_dense_spatial_action_from_independent_edge_energy():
    radii=[F(7,8),F(9,10),F(1),F(6,5),F(11,10)]
    k,raw,_,_=waves.wave_matrix(radii)
    direct=np.zeros((5,5));values=[1.,*map(float,radii),1.]
    for edge in range(6):
        row=np.zeros(5)
        if edge>0:row[edge-1]=-1
        if edge<5:row[edge]=1
        direct+=(values[edge]**2+values[edge+1]**2)/2*np.outer(row,row)
    np.testing.assert_allclose(raw,direct,atol=1e-15,rtol=0)
    expected=np.diag(np.sqrt(list(map(float,radii))))@direct@np.diag(np.sqrt(list(map(float,radii))))
    np.testing.assert_allclose(k,expected,atol=1e-15,rtol=0)
    assert np.linalg.eigvalsh(k).min()>0


def test_all_modes_matrix_function_and_complete_probability_records():
    radii=[F(9,10),F(1),F(11,10),F(21,20),F(1),F(19,20),F(1)]
    duration=3
    result=waves.execute(radii,duration)
    k=waves.wave_matrix(radii)[0]
    direct=sqrtm(k)
    x=np.arange(7);state=np.exp(-(x-7)**2/8+1j*.7*x);state/=np.linalg.norm(state)
    final=expm(-1j*duration*direct)@state
    np.testing.assert_allclose(result['complete_vacuum_and_position_probabilities'],[.5,*list(.5*np.abs(final)**2)],atol=2e-14,rtol=0)
    assert len(result['normal_frequencies'])==7
    assert abs(sum(result['complete_vacuum_and_position_probabilities'])-1)<1e-13


def test_total_one_quantum_cutoff_keeps_exact_oscillator_equations():
    h=s.diag(0,1,3)
    for index,frequency in ((1,s.Integer(1)),(2,s.Integer(3))):
        annihilator=s.zeros(3);annihilator[0,index]=1
        q=(annihilator+annihilator.T)/s.sqrt(2*frequency)
        p=s.I*s.sqrt(frequency/2)*(annihilator.T-annihilator)
        assert s.simplify(s.I*(h*q-q*h)-p)==s.zeros(3)
        assert s.simplify(s.I*(h*p-p*h)+frequency**2*q)==s.zeros(3)


def test_metric_delay_radial_control_and_Einstein_counterexample():
    result=geometry.metric_certificate()
    assert all(F(x)>0 for x in result['piecewise_metric_transverse_delay_interval'])
    assert result['radial_null_coordinate_speed_squared']=='1'
    value=geometry.einstein_test()
    assert value['smooth_harmonic_q_minus_inverse_radius_G00_at_3_0_0']=='-1/9'
    assert abs(F(value['discrete_linearized_G00']))>1000000*abs(F(value['discrete_laplacian_q']))
    assert not value['same_source_vacuum_Einstein_premise_passed']


def test_Taylor_tail_and_deleted_geometry_control():
    assert waves.operator_tail(3,F(16),16)<F(1,10**21)
    _,radii,_=geometry.prepared()
    curved=waves.execute(radii,12);flat=waves.execute([F(1)]*31,12)
    distance=np.abs(np.array(curved['complete_vacuum_and_position_probabilities'])-flat['complete_vacuum_and_position_probabilities']).sum()/2
    assert distance>1e-6
    assert curved['one_quantum_mean_position']<flat['one_quantum_mean_position']
    with pytest.raises(ValueError):waves.wave_matrix([1,-1,1])
    with pytest.raises(ValueError):waves.operator_tail(3,100,1,degree=3)
