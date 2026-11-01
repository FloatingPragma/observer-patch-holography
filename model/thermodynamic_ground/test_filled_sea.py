from fractions import Fraction as F
from itertools import combinations

import numpy as np

from model.interacting_vacuum import ground
from model.thermodynamic_ground import filled_sea


def test_complete_source_sea_has_exact_nonzero_energy_variance():
    result=filled_sea.certificate()
    assert result['particle_number']==8
    assert result['nonzero_occupation_terms']==3844
    assert result['contact_mean']=='2'
    assert result['contact_variance']==result['interacting_energy_variance']=='2/3'
    assert result['interacting_stationary_pure_state'] is False
    assert result['correlated_interacting_ground_excluded'] is False


def test_independent_small_wick_determinant_formula_and_complete_sum():
    adjacency=ground.band_data()[0]
    p=np.kron((np.eye(8)-adjacency/np.sqrt(3))/2,np.eye(2))
    assert np.linalg.norm(p@p-p)<1e-14
    mean=sum(np.linalg.det(p[np.ix_([2*x,2*x+1],[2*x,2*x+1])]) for x in range(8))
    second=mean
    for x,y in combinations(range(8),2):
        selected=[2*x,2*x+1,2*y,2*y+1]
        second+=2*np.linalg.det(p[np.ix_(selected,selected)])
    assert abs(mean-2)<1e-14
    assert abs(second-mean**2-float(F(2,3)))<1e-14


def test_zero_contact_control_has_no_energy_variance():
    from model.cube_eta_protection import cube
    a,b=filled_sea.reference_arrays()
    d=np.asarray([sum(((mask>>(2*j))&3)==3 for j in range(8))
                  for mask in range(65536)],dtype=np.int64)
    h=cube.operators().hubbard_H.real
    assert np.array_equal(h@a-16*d*a,-24*b)
    assert np.array_equal(h@b-16*d*b,-8*a)
    # Kinetic-only state is stationary; contact creates the variance.
    assert not np.array_equal(d*a,2*a) or not np.array_equal(d*b,2*b)


def test_sea_arrays_are_fresh_and_not_a_retained_infinite_state():
    first=filled_sea.reference_arrays();second=filled_sea.reference_arrays()
    first[0][:]=0
    assert np.count_nonzero(second[0])>0
    assert second[0].shape==(65536,)
