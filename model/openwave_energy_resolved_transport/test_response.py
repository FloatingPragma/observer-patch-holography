"""Independent finite matrices test the full-H projector and time bounds."""
from fractions import Fraction

import numpy as np
import pytest
import scipy.linalg as la
import sympy as s

from .response import (native_pole_residue,tagged_spin_effects,projector_tail_bounds,
                       source_bounds,timeline,record_law)


def test_actual_pole_and_spin_contraction():
    residue=native_pole_residue()
    assert residue['native_removed_spin_target_dimensions']==[17,7]
    assert residue['actual_selected_target_spin']==1
    assert residue['exact_residue_difference_divisible_by_selected_polynomial']
    assert tagged_spin_effects()['all_three_axes_exact']


@pytest.mark.parametrize('seed',range(12))
def test_projector_tail_against_independent_complex_matrix_evolution(seed):
    rng=np.random.default_rng(seed)
    raw=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    coupling=raw+raw.conj().T
    np.fill_diagonal(coupling,0)
    coupling/=np.linalg.norm(coupling,2)
    lam=s.Rational(1,500)
    h0=np.diag([0.,3.,1.,5.]);h=h0+float(lam)*coupling
    energies,vectors=np.linalg.eigh(h)
    selected=(abs(energies)<.1)
    assert selected.sum()==1
    projector=vectors[:,selected]@vectors[:,selected].conj().T
    y=np.array([np.sqrt(.3),np.sqrt(.2),0,0],complex)
    x0=np.array([np.sqrt(.3),0,0,0],complex)
    moment=abs(coupling[2,0])**2*.3
    # Literal float-to-rational moment widened independently for roundoff.
    rational_m=s.Rational(float(moment))
    margin=s.Rational(1,10**14)
    bounds=projector_tail_bounds(coupling=lam,tau=2,edge_norm=2,
        weight_intervals=[(s.Rational(3,10),)*2,(s.Rational(1,5),)*2],
        distances=[s.Rational(1,10),s.Rational(29,10)],
        tagged_energy_separation=s.Rational(9,10),pole_to_other_energy_max=[0,3],
        trial_error=0,moment_interval=(max(0,rational_m-margin),rational_m+margin))
    state=projector@y
    remainder=state[2]+float(lam)*(coupling@x0)[2]
    assert abs(remainder)<=float(bounds['projection_remainder_upper'])+1e-14
    lower,upper=map(float,bounds['initial_tagged_amplitude_interval'])
    assert lower-1e-14<=abs(state[2])<=upper+1e-14
    assert abs(state[2])>0
    assert y[2]==0  # Removing the energy filter removes the t0 signal.
    assert abs(np.vdot(state,h@state))<=.1*np.vdot(state,state).real+1e-14
    assert h[2,2]==1  # The tagged output has a different energy expectation.
    for t in (0,.1,1,20):
        evolved=la.expm(-1j*t*h)@state
        limit=min(float(bounds['tagged_phase_adjusted_amplitude_change_absolute_upper']),
                  t*float(bounds['tagged_phase_adjusted_amplitude_change_per_time_upper']))
        assert abs(evolved[2]-state[2])<=limit+1e-13
        # Complete preparation/pump/bin/arrival law, with q=3/5 and P0 certain.
        q=.6;b=np.vdot(state,state).real;p=abs(evolved[2])**2
        outcomes=[1-q,0,q/2,q*(.5-b),q*p,q*(b-p)]
        assert min(outcomes)>-1e-13 and abs(sum(outcomes)-1)<1e-12


def test_wrong_tail_sign_and_omitted_outside_atom_are_detected():
    lam=.0001
    h0=np.diag([0.,3.,1.])
    perturb=np.array([[0.,1.,1.],[1.,0.,0.],[1.,0.,0.]])
    h=h0+lam*perturb
    values,vectors=np.linalg.eigh(h)
    p=vectors[:,abs(values)<.1]@vectors[:,abs(values)<.1].T
    x0=np.array([1.,0,0])
    leading=lam*(perturb@x0)[2]
    correct=abs((p@x0)[2]+leading)
    wrong=abs((p@x0)[2]-leading)
    assert correct<1e-10 and wrong>lam
    outside=p@np.array([0.,1.,0.])
    assert outside[2]!=0  # This O(lambda^2) component cannot be discarded.
    assert abs(outside[2])<10*lam**2


def test_native_signal_is_positive_before_any_clock_evolution():
    values=source_bounds()
    assert values['bounds']['strict_initial_amplitude_margin']>0
    rows=timeline()
    assert rows[0]['source_time']=='0'
    assert Fraction(rows[0]['per_route_original_input_probability_interval'][0])>0
    assert Fraction(rows[0]['per_route_probability_change_from_t0_absolute_upper'])==0
    assert all(not row['positive_time_change_or_propagation_certified'] for row in rows)
    assert record_law()['coarse_record_count']==17


def test_missing_bin_distance_cannot_supply_the_projector_bound():
    with pytest.raises(ValueError):
        projector_tail_bounds(coupling=s.Rational(1,100),tau=1,edge_norm=1,
            weight_intervals=[(s.Rational(1,2),)*2],distances=[0],
            tagged_energy_separation=1,pole_to_other_energy_max=[0],trial_error=0,
            moment_interval=(s.Rational(1,4),)*2)
