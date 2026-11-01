from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import bounds, cube, dimer


def test_full_source_particle_hole_symmetry():
    result=cube.source_symmetry()
    assert result['full_Fock_dimension']==65536
    assert result['all_matrix_entries_checked'] is True


def test_exact_source_hypotheses_of_imported_Lieb_theorem():
    result=cube.lieb_premises()
    assert result['connected_vertices']==list(range(8))
    assert result['bipartition_sizes']==[4,4]
    assert result['even_half_filling'] is True
    assert result['ground_hidden_spin']==0 and result['ground_multiplicity']==1
    assert result['proof_imported_not_reproved'] is True
    assert result['uniform_hopping_sign_required'] is False


def test_unique_clock_and_old_ground_failure_are_exact():
    U,kappa,V,N=s.symbols('U kappa V N',real=True)
    defect=s.expand(U*(V-N)-kappa*(2*V-N)+kappa*N)
    assert s.expand(defect-(U-2*kappa)*(V-N))==0
    assert s.solve(defect,kappa)==[U/2]
    assert -2-s.sqrt(3)/4 > -4
    row=cube.sector_bounds()
    assert row['centered_global_ground_number']==8
    assert row['old_N4_ground_remains_centered_ground'] is False
    assert s.sqrt(17)/2-s.sqrt(3)/4-s.Rational(3,2)>s.Rational(1,8)


def test_full_dimer_spectrum_and_unique_ground():
    data=dimer.data()
    eig=data['K'].eigenvals()
    ground=data['ground_energy']
    assert sum(eig.values())==16 and eig[ground]==1
    assert all(value>ground for value in eig if value!=ground)
    assert s.simplify(data['N']*data['ground']-2*data['ground'])==s.zeros(16,1)


def test_actual_dimer_charge_measures_and_faithfulness():
    data=dimer.data();rho=data['local_density']
    assert s.trace(rho)==1 and all(x>0 for x in rho.diagonal())
    lines=data['charged_lines']
    assert s.simplify(sum(row['addition_weight'] for row in lines)-s.Rational(1,2))==0
    assert all(row['addition_weight']==row['removal_weight'] for row in lines)
    second=s.simplify(sum(2*row['addition_weight']*row['energy']**2 for row in lines))
    assert second==s.Rational(65,256)==data['signed_second_moment']


def test_zero_coupling_retains_the_null_upper_line():
    data=dimer.data(t=s.Rational(1,16),U=s.Integer(0))
    assert data['charged_lines'][1]['addition_weight']==0
    assert data['charged_lines'][0]['energy']==s.Rational(1,16)
    assert data['local_density']==s.eye(4)/4


def test_canonical_composite_force_norm_premise():
    c=[dimer.annihilator(j) for j in range(4)]
    n=[field.T*field for field in c]
    identity=s.eye(16)
    Z=(identity-2*n[0]*n[1])*(identity-2*n[2]*n[3])
    d=[(2*n[j^1]-identity)*c[j] for j in range(4)]
    assert Z*Z==identity
    for j in range(4):
        assert Z*c[j]*Z==-d[j]
        for k in range(4):
            assert d[j]*d[k]+d[k]*d[j]==s.zeros(16)
            assert d[j]*d[k].T+d[k].T*d[j]==(identity if j==k else s.zeros(16))


def test_dimer_product_upper_bound_and_bond_sign():
    positive=dimer.data()
    # Flipping the orbital phase on site one reverses hopping, preserves N,D.
    phase=s.diag(*[(-1)**((mask>>2).bit_count()) for mask in range(16)])
    negative=phase*positive['K']*phase
    state=phase*positive['ground']
    assert s.simplify(negative*state-positive['ground_energy']*state)==s.zeros(16,1)
    assert 4*positive['ground_energy']==-2-s.sqrt(17)/2
    # Fixed factor number makes each one-field expectation vanish exactly.
    for mode in range(4):
        assert (state.T*dimer.annihilator(mode)*state)[0]==0


def test_pullthrough_constants_and_refinement_rate():
    assert F(22,7)*F(7,4)/2==F(11,4)
    assert bounds.pull_through(F(1,8192))==F(11,536870912)
    assert bounds.original_pull_through(F(1,8192))==F(11,268435456)
    assert bounds.pull_through(F(2))==2
    assert bounds.pull_through(F(0))==0


def test_Ritz_diagnostic_remains_within_proved_bounds():
    data=cube.numerical_diagnostic()
    low=-4-np.sqrt(3)/4;upper=-2-np.sqrt(17)/2
    assert low<data['lowest_two_Ritz_values'][0]<upper
    assert data['balanced_N8_dimension']==4900
    assert max(data['residual_norms'])<1e-8
    assert data['used_for_ground_or_gap_admission'] is False


@pytest.mark.parametrize('value',[-1,0.1,True])
def test_invalid_exact_ratios_rejected(value):
    with pytest.raises(ValueError):bounds.pull_through(value)
    with pytest.raises(ValueError):bounds.original_pull_through(value)
