"""Independent exact source and measure-scaling controls; no physical mass fit."""
from fractions import Fraction as F
import math

import pytest
import sympy as s

from model.source_family_compatibility import mass_scale_selection as m
from model.source_family_compatibility import local_spectral
from model.source_family_compatibility.defect_budget import exact_source


@pytest.fixture(scope='module')
def source_control():
    return m.default_source_scaling(F(3, 2))


@pytest.fixture(scope='module')
def receipt():
    return m.certificate()


@pytest.mark.parametrize('a,v,ratio', [(F(1, 16), 1, 8), (F(3, 7), F(5, 2), F(11, 3))])
def test_literal_source_hopping_and_physical_clock_scale(a, v, ratio):
    row=m.physical_source(a,v,ratio)
    assert row['c']==v/a and row['U']==ratio*v/a
    assert row['literal_hopping_magnitude']==v/(2*a)
    assert row['physical_energy_per_normalized_energy']==row['c']
    assert not row['a_dependent_reference_selected']
    transformed=m.scale_coordinates((F(1,3),F(-4,7),0),F(-3,2),F(5,4))
    assert transformed['momentum']==(F(5,12),F(-5,7),F(0))
    assert transformed['energy']==F(-15,8)


def test_jacobian_and_inverse_spacing_by_independent_exact_integration():
    # Fixed normalized fibre delta_(1+q1); compact tests are quadratic on
    # the bounded energy range in question and vanish outside a larger set.
    # Exact integration over a rectangular momentum support tests the3D
    # Jacobian and distinguishes a/lambda from lambda*a.
    x=s.Symbol('x',real=True)
    a=s.Rational(3,7); lam=s.Rational(5,3)
    sides=(s.Rational(2,3),s.Rational(5,4),s.Rational(7,6))
    lhs=s.integrate((lam*(1/a+x))**2,(x,0,sides[0]/lam))*sides[1]*sides[2]/lam**2
    def integral(spacing):
        return s.integrate((1/spacing+x)**2,(x,0,sides[0]))*sides[1]*sides[2]
    rhs=integral(a/lam)/lam**3
    assert s.simplify(lhs-rhs)==0
    assert s.simplify(lhs-integral(a*lam)/lam**3)!=0
    assert s.simplify(lhs-integral(a/lam)*lam**3)!=0


def test_actual_full_default_source_clock_projectors_and_original_weights(source_control):
    h,rho,a,projectors=local_spectral.dimer_source()
    scaled,_,_=exact_source(2,((0,1,1),),F(3,16),F(3,2))
    simple=lambda x:x.applyfunc(s.simplify)
    assert simple(scaled-s.Rational(3,2)*h)==s.zeros(16)
    assert simple(rho*rho-rho)==s.zeros(16) and s.trace(rho)==1
    records=source_control['complete_charge_energy_records']
    assert len(records)==2*len(projectors)
    assert sum(s.sympify(r['original_probability']) for r in records)==1
    assert sum(r['original_probability']=='0' for r in records)>0
    total=s.zeros(16)
    for label,creator,sign in (('addition',a.H,1),('removal',a,-1)):
        rows=[r for r in records if r['charge']==label]
        assert len(rows)==len(projectors)
        for row,(energy,projection) in zip(rows,projectors,strict=True):
            effect=simple(creator.H*projection*creator)
            total+=effect
            assert s.simplify(s.sympify(row['original_probability'])-s.trace(rho*effect))==0
            assert s.simplify(s.sympify(row['signed_source_energy'])-sign*energy)==0
            assert s.simplify(s.sympify(row['signed_scaled_energy'])-sign*s.Rational(3,2)*energy)==0
        assert s.simplify(sum(s.sympify(r['original_probability'])for r in rows)-s.Rational(1,2))==0
    assert simple(total-s.eye(16))==s.zeros(16)
    assert source_control['same_full_ground_projector']
    assert source_control['zero_probability_records_retained']
    assert not source_control['conditional_probability_renormalization_applied']
    assert not source_control['finite_dimer_is_continuum_momentum_mode']


def test_original_probability_can_escape_without_being_renormalized():
    # A measure implication control: fixed mass z at zero plus two distant
    # atoms. It demonstrates a subprobability vague limit and retained
    # original complement, without claiming a realization of the CAR source.
    z=F(2,5);source=((F(0),z),(F(1),(1-z)/2),(F(-1),(1-z)/2))
    for spacing in (F(1,2),F(1,4),F(1,8)):
        scale=m.physical_source(spacing)['physical_energy_per_normalized_energy']
        physical=tuple((scale*e,w)for e,w in source)
        inside=sum(w for e,w in physical if abs(e)<=1)
        outside=sum(w for e,w in physical if abs(e)>1)
        assert inside==z and outside==1-z and inside+outside==1
        assert inside!=1


@pytest.mark.parametrize('power',[1,2,3,7])
@pytest.mark.parametrize('spacing',[F(1),F(7,100),F(1,256),F(17,10000)])
def test_dense_refinement_exact_bracket_and_inverse_test_map(power,spacing):
    row=m.dense_power_bracket(spacing,power)
    n=row['n']
    assert row['a_n']==F(1,n**power)
    assert row['a_next']==F(1,(n+1)**power)<spacing<=row['a_n']
    assert row['lambda']==spacing/row['a_n']
    assert row['lambda_lower']<row['lambda']<=1
    assert row['inverse_test_dilation']*row['lambda']==1
    assert row['jacobian']==row['lambda']**-3


def test_dense_interpolation_uses_common_compact_uniform_measure_bound():
    # For a_n=n^-2 the relative bracket, rather than absolute mesh width,
    # controls the rescaled compact test. The latter matters at zero scale.
    defects=[]
    for n in (10,100,1000):
        spacing=(F(1,n*n)+F(1,(n+1)**2))/2
        row=m.dense_power_bracket(spacing,2)
        assert row['n']==n
        inverse=row['inverse_test_dilation']
        assert 1<=inverse<4
        # f(p,E)=p1 on a common compact gives an exact nonzero test error.
        error=(inverse-1)*F(3,2)
        assert error<=((F(n+1,n)**2)-1)*F(3,2)
        defects.append(error)
    assert defects[2]<defects[1]<defects[0]


def _dyadic_energy(radius):
    """Continuous exact sparse-limit implication control, not source evidence."""
    radius=F(radius)
    if not radius:return F(0)
    assert radius>0
    scale=F(1)
    while radius<scale:scale/=2
    while radius>=2*scale:scale*=2
    t=radius/scale
    # A bounded-Lipschitz continuous dyadically homogeneous profile.
    profile=2*t-1 if t<=F(3,2) else F(2)
    return scale*profile


def test_sparse_refinement_does_not_imply_one_continuously_dilation_invariant_limit(receipt):
    first=[];second=[]
    for n in (3,7,19):
        a=F(1,2**n);b=F(3,2)*a
        first.append(_dyadic_energy(a)/a)
        second.append(_dyadic_energy(b)/b)
        assert _dyadic_energy(2*a)==2*_dyadic_energy(a)
    assert first==[1,1,1] and second==[F(4,3)]*3
    assert not receipt['arbitrary_sparse_subsequence_excluded_by_theorem']
    assert receipt['dense_refinement_ratio_tending_to_one_suffices']
    assert not receipt['actual_full_continuum_limit_proved']


def test_positive_mass_shell_dilates_are_disjoint_inside_one_compact():
    factors=tuple(F(1)+F(j,100)for j in range(1,101))
    row=m.shell_dilates((1,0,0),2,1,factors)
    assert row['original_invariant_mass_squared']==3
    assert len({r['invariant_mass_squared']for r in row['dilates']})==100
    for r in row['dilates']:
        assert max(abs(x)for x in r['momentum'])<=2
        assert abs(r['energy'])<=4
        assert r['invariant_mass_squared']==3*r['factor']**2
        assert r['joint_measure_weight_factor']==r['factor']**3>=1
    # A compact positive shell piece of weightw would have each scaled
    # copy of weight at leastw. The momentum marginal bounds the union.
    compact_mass_upper=F(4)**3;w=F(3,4)
    assert len(factors)*w>compact_mass_upper
    assert row['hypothetical_shell_geometry_only'] and not row['actual_source_shell_proved']


def test_massless_shell_cannot_use_positive_mass_disjointness_argument():
    with pytest.raises(ValueError,match='positive invariant mass'):
        m.shell_dilates((1,0,0),1,1,(1,F(3,2),2))
    # Dilated massless annuli overlap on the same cone.
    assert max(F(1),F(3,2))<min(F(2),F(3))


def test_continuous_timelike_response_is_allowed_by_homogeneity(receipt):
    # F_p is uniform on [2|p|,3|p|], so its mass is one, its density scales
    # as1/|p|, and it has no atom on any fixed invariant mass shell.
    radius=F(3,5);factor=F(7,4)
    density=1/radius
    assert density*(3*radius-2*radius)==1
    scaled_density=1/(factor*radius)
    assert scaled_density*factor==density
    assert not receipt['continuous_timelike_response_excluded']
    assert not receipt['every_unstable_resonance_excluded']


def test_running_mass_control_breaks_exact_fixed_family_premise(receipt):
    # A supplied running dimensionless coefficient yields a finite physical
    # mass; this countercontrol does not derive that coefficient or mass.
    mass=F(3,7);speed=F(5,4);momentum=F(2,3)
    values=[]
    for spacing in (F(1,8),F(1,16),F(1,32)):
        mu=spacing*mass/speed
        c=m.physical_source(spacing,speed)['c']
        assert c*mu==mass
        q=float(spacing*momentum)
        # A two-generator Clifford slice with a supplied anticommuting mass.
        pauli_x=s.Matrix([[0,1],[1,0]]);pauli_z=s.diag(1,-1)
        lattice_coefficient=s.sin(s.Rational((spacing*momentum).numerator,(spacing*momentum).denominator))
        symbol=lattice_coefficient*pauli_x+mu*pauli_z
        assert symbol*symbol==(lattice_coefficient**2+mu**2)*s.eye(2)
        energy2=float(c)**2*(math.sin(q)**2+float(mu)**2)
        limit=float(speed*speed*momentum*momentum+mass*mass)
        assert abs(energy2-limit)<=float(speed*speed*spacing*spacing*momentum**4/3)
        values.append(mu)
    assert len(set(values))==3
    assert receipt['fixed_reference_and_dimensionless_couplings_required']
    assert not receipt['running_reference_or_observable_dictionary_excluded']
    assert receipt['physical_species_assigned']==0 and receipt['physical_mass_assigned']is None


@pytest.mark.parametrize('bad',[0,-1,True,.5,'1'])
def test_invalid_source_scale_and_clock_inputs_rejected(bad):
    with pytest.raises(ValueError):m.physical_source(bad)
    with pytest.raises(ValueError):m.scale_coordinates((0,0,0),1,bad)


@pytest.mark.parametrize('bad_power',[0,-1,True,F(2),1.5])
def test_dense_sequence_contract_rejects_unproved_refinement_input(bad_power):
    with pytest.raises(ValueError):m.dense_power_bracket(F(1,10),bad_power)


def test_momentum_and_shell_input_contracts():
    with pytest.raises(ValueError):m.scale_coordinates((1,2),1,2)
    with pytest.raises(ValueError):m.shell_dilates((1,0,0),2,1,(1,1))
    with pytest.raises(ValueError):m.threshold_dilation((2,0,0),1,1,F(1,10))
    row=m.threshold_dilation((1,0,0),2,1,F(1,100))
    assert row['invariant_mass_squared']<F(1,10000)
    assert row['channel_support_dilation_invariance_required']
    assert not row['actual_source_causal_cone_proved']


def test_direct_helper_fraction_output_cannot_mutate_future_default():
    first=m.default_source_scaling()
    expected=F(3,2)
    first['positive_clock_factor']._numerator=17
    second=m.default_source_scaling()
    assert second['positive_clock_factor']==expected


def test_certificate_keeps_exact_conditional_scope_and_detached_probabilities(receipt):
    assert receipt['full_joint_vague_limit_required']
    assert receipt['positive_original_weight_positive_mass_shell_excluded_if_premises']
    assert not receipt['actual_nodal_atom_positivity_proved']
    assert not receipt['actual_default_source_homogeneity']['native_large_detector_tape_executed']
    original=F(receipt['physical_source_example']['c'].numerator,receipt['physical_source_example']['c'].denominator)
    receipt['physical_source_example']['c']._numerator=999
    fresh=m.certificate()
    assert fresh['physical_source_example']['c']==original
    assert fresh['actual_default_source_homogeneity']['total_original_probability']==1


@pytest.mark.parametrize('ratio,index', [(F(1,2),0),(F(1,2),12),(F(2,3),7)])
def test_geometric_refinement_retains_only_discrete_covariance(ratio,index):
    row=m.geometric_refinement(ratio,index)
    assert row['a_n']==ratio**index and row['a_next']==ratio*row['a_n']
    assert row['jacobian']==ratio**-3
    assert not row['continuous_dilation_covariance_inferred']
    assert row['discrete_covariance_requires_whole_sequence_limit']


def _dyadic_shell_measure(radius):
    """Normalized continuous sparse-scaling control, not a source prediction."""
    radius=F(radius)
    if not radius:return ((F(0),F(1)),)
    assert radius>0
    scale=F(1)
    while radius<scale:scale/=2
    while radius>=2*scale:scale*=2
    t=radius/scale
    # Entries are (positive energy squared, original weight). Half of each
    # weight can be placed at each energy sign to obtain an even measure.
    return ((radius**2+scale**2,2-t),(radius**2+4*scale**2,t-1))


def test_geometric_shell_tower_is_allowed_without_a_full_continuous_limit():
    for radius in (F(1),F(3,2),F(7,5),F(5,3)):
        original=_dyadic_shell_measure(radius)
        assert sum(w for _,w in original)==1 and all(w>=0 for _,w in original)
        for n in (1,4,12):
            ratio=F(1,2**n)
            actual=_dyadic_shell_measure(ratio*radius)
            assert actual==tuple((ratio*ratio*energy2,weight)for energy2,weight in original)
            assert all(energy2>0 for energy2,_ in actual)
    # Every dyadic shell appears on an open annular momentum region.
    for n in (0,1,7):
        scale=F(1,2**n);p=F(3,2)*scale
        row=_dyadic_shell_measure(p)
        assert row[0][1]==F(1,2) and row[0][0]-p*p==scale*scale
    # The different3/2-scaled refinement has a different tangent measure.
    first=sum(e*w for e,w in _dyadic_shell_measure(1))
    second=sum(e*w for e,w in _dyadic_shell_measure(F(3,2)))/F(3,2)**2
    assert first==2 and second==F(19,9) and first!=second
    # Geometric copies approaching zero can have summable joint weights.
    ratio=F(1,2);piece=F(3,4)
    assert sum(piece*ratio**(3*j)for j in range(100))<piece/(1-ratio**3)
    assert ratio**20<F(1,10**6)


def test_dyadic_shell_control_is_weakly_continuous_at_annular_boundaries():
    for epsilon in (F(1,10),F(1,100),F(1,1000)):
        left=_dyadic_shell_measure(1-epsilon)
        right=_dyadic_shell_measure(1+epsilon)
        # The nonmatching mass branch loses its original weight, whereas
        # the matching energy approaches sqrt2 from either side.
        assert left[0][1]==2*epsilon and right[1][1]==epsilon
        assert left[1][0]==(1-epsilon)**2+1
        assert right[0][0]==(1+epsilon)**2+1
    for n in (5,10,20):
        radius=F(1,2**n)
        assert max(e for e,_ in _dyadic_shell_measure(radius))<=5*radius**2


@pytest.mark.parametrize('ratio,index', [(0,0),(1,0),(F(3,2),0),(F(1,2),-1),(F(1,2),True)])
def test_geometric_refinement_input_contract(ratio,index):
    with pytest.raises(ValueError):m.geometric_refinement(ratio,index)
